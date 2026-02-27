function makeOrientationVideoFromMotionData(data, tsec, opt, EventRecord)
% makeOrientationVideoFromMotionData  Render orientation video + event cue + N/E labels.
%
% Inputs
%   data        struct with fields acc, gyr, magcorr
%   tsec        timestamps (seconds)
%   opt         options struct with FolderProcDataMat
%   EventRecord struct with event info
%
% Jesus 05.02.2026

if nargin < 4 || isempty(EventRecord)
    eventTsec.stimOn1 = [];  % row vector
    eventTsec.stimOn2 = [];  % row vector
    eventTsec.rwd = [];  % row vector
else
    eventTsec.stimOn1 = EventRecord.TimeSecFromMidnight(EventRecord.EventType==1).';
    eventTsec.stimOn2 = EventRecord.TimeSecFromMidnight(EventRecord.EventType==2).';
    eventTsec.rwd = EventRecord.TimeSecFromMidnight(EventRecord.EventType==7).';

    % eventTsec.stimOn1 = eventTsec.stimOn1(:).';  % row vector
    % eventTsec.stimOn2 = eventTsec.stimOn2(:).';  % row vector
    % eventTsec.rwd = eventTsec.rwd(:).';  % row vector
end

% Settings
Fs  = opt.fsmot;
fps.fuse = 100; % calculate fuse at HiRes
fps.vid = 20; % will plot only a subset for video purpouses
hop = max(1, round(Fs/fps.fuse)); % obtain sample skipping

outFile = fullfile(opt.FolderProcDataMat, "FUSE_result.mp4");

% 1) Gather raw signals
acc0 = [data.acc.X(:),     data.acc.Y(:),     data.acc.Z(:)];
gyr0 = [data.gyr.X(:),     data.gyr.Y(:),     data.gyr.Z(:)];
mag0 = [data.magcorr.X(:), data.magcorr.Y(:), data.magcorr.Z(:)];

gyr0 = deg2rad(gyr0);

% 2) Axis remap
accB = [-acc0(:,2), acc0(:,3),  acc0(:,1)];
gyrB = [-gyr0(:,2), gyr0(:,3),  gyr0(:,1)];
magB = [-mag0(:,1), mag0(:,3),  mag0(:,2)];

% Downsample for fusion & render
idx  = 1:hop:numel(tsec);
tVid = tsec(idx);
FsD = 1/median(diff(tVid));

accD = accB(idx,:);
gyrD = gyrB(idx,:);
magD = magB(idx,:);

% --- Magnetometer cleaning (spike gate + fill) ---
magNorm = vecnorm(magD,2,2);
base    = movmedian(magNorm, round(.5*FsD));      % 2 s baseline
ratio   = magNorm ./ max(base,eps);

% Flag strong spikes (tune thresholds)
isBad = ratio > 1.25 | ratio < 0.80;

% Fill per-axis outliers with moving median (robust)
magClean = magD;
for j = 1:3
    magClean(:,j) = filloutliers(magD(:,j), "spline", round(0.5*FsD));
end

% Hard gate: hold last good vector where norm is implausible
for k = 2:size(magClean,1)
    if isBad(k)
        magClean(k,:) = magClean(k-1,:);
    end
end

magD = magClean;

% 3) Orientation estimation (AHRS)
fuse = ahrsfilter("SampleRate", FsD);
fuse.GyroscopeNoise = 0.01;         % (rad/s)²
fuse.GyroscopeDriftNoise = 1e-6;    % (rad/s)²
fuse.MagneticDisturbanceNoise = 300; % (µT)²
fuse.LinearAccelerationDecayFactor = 0.1; 
fuse.MagneticDisturbanceDecayFactor = 0.1; %
fuse.ExpectedMagneticFieldStrength = 25;

n = size(accD,1);
q = quaternion.zeros(n,1);

% Optional: store bias estimates to verify drift
gyroBias = zeros(n,3);

warning off

for i = 1:n
    [q(i), gyroBias(i,:)] = fuse(accD(i,:), gyrD(i,:), magD(i,:));
    % gyroBias(i,:) = fuse.GyroscopeBias;
end

if ~isa(q,"quaternion")
    q = quaternion(q);
end

% Rotation matrices for hgtransform
try
    R = rotmat(q,"point");   % 3x3xN or "frame"
catch
    R = quat2rotm(compact(q));
    R = permute(R,[2 3 1]);
    if size(R,3) ~= numel(tVid)
        R = permute(quat2rotm(compact(q)), [2 3 1]);
    end
end

% 4) Scene setup
figW = 640; figH = 640;
fig = figure( ...
    "Visible","on", ...
    "Color","w", ...
    "Renderer","opengl", ...
    "Units","pixels", ...
    "Position",[100 100 figW figH]);

ax = axes(fig);
hold(ax,"on");
grid(ax,"on");
axis(ax,"equal");
axis(ax,[-1 1 -1 1 -1 1]);
axis(ax,"manual");
view(ax,3);
xlabel(ax,"World X");
ylabel(ax,"World Y");
zlabel(ax,"World Z");

% Fixed world triad at origin
plot3(ax,[0 1],[0 0],[0 0],"k-","LineWidth",1);
plot3(ax,[0 0],[0 1],[0 0],"k-","LineWidth",1);
plot3(ax,[0 0],[0 0],[0 1],"k-","LineWidth",1);

% Reference ring in XY plane
th = linspace(0,2*pi,200);
plot3(ax,0.75*cos(th),0.75*sin(th),0*th,"Color",[0.7 0.7 0.7]);

% Magnetic orientation letters on planes X=+1 and Y=+1
% Assumption: world +X is North, world +Y is East (NEU visualization cue)
text(ax,  1, 0, 0, "N", "FontSize",18, "FontWeight","bold", ...
    "HorizontalAlignment","left", "VerticalAlignment","middle", "Color",[0 0 0]);
text(ax,  0, 1, 0, "E", "FontSize",18, "FontWeight","bold", ...
    "HorizontalAlignment","center", "VerticalAlignment","bottom", "Color",[0 0 0]);
text(ax, -1, 0, 0, "S", "FontSize",18, "FontWeight","bold", ...
    "HorizontalAlignment","left", "VerticalAlignment","middle", "Color",[0 0 0]);
text(ax,  0,-1, 0, "W", "FontSize",18, "FontWeight","bold", ...
    "HorizontalAlignment","center", "VerticalAlignment","bottom", "Color",[0 0 0]);

% 5) Agent geometry (sphere + cone + axes)
Tnode = hgtransform("Parent",ax);

% Sphere
rs = 0.15;
[nx,ny,nz] = sphere(30);
surf(rs*nx, rs*ny, rs*nz, ...
    "Parent", Tnode, ...
    "FaceColor",[0.2 0.2 0.2], ...
    "EdgeColor","none", ...
    "FaceLighting","gouraud", ...
    "AmbientStrength",0.25);

% Cone nose along +Xb
Ln = 0.20;
rb = 0.06;
nCirc = 40;
phi = linspace(0,2*pi,nCirc);

xBase = rs*ones(1,nCirc);
yBase = rb*cos(phi);
zBase = rb*sin(phi);
tip   = [rs+Ln; 0; 0];

V = [[xBase; yBase; zBase], tip].';
F = [(1:nCirc).', [2:nCirc,1].', repmat(nCirc+1,nCirc,1)];

patch( ...
    "Parent", Tnode, ...
    "Vertices", V, ...
    "Faces", F, ...
    "FaceColor",[0.85 0.2 0.2], ...
    "EdgeColor","none", ...
    "FaceLighting","gouraud", ...
    "AmbientStrength",0.25);

% Body axes (green is now both directions and half-length)
L = 1;
L2 = 0.5*(L/2);  % half of current length, split both sides => total length = L/2
% Uncomment red if you want forward axis for debugging:
% line([-L2 L2],[0 0],[0 0], "Parent",Tnode, "Color",[0.85 0.1 0.1], "LineWidth",2); % Xb
line([0 0],[-L2 L2],[0 0], "Parent",Tnode, "Color",[0.1 0.6 0.1], "LineWidth",2);   % Yb (L/R)
line([0 0],[0 0],[0 L/4], "Parent",Tnode, "Color",[0.1 0.3 0.9], "LineWidth",2);    % Zb up (keep as-is)

camlight(ax,"headlight");
material(ax,"dull");

% 5b) Event cue (top-right), lasts 2 seconds
% Use an annotation so it stays in the corner regardless of 3-D view.
hStim1 = annotation(fig, "textbox", [0.73 0.83 0.2 0.05], ...
    "String", "Stim1", ...
    "FitBoxToText","off", ...
    "HorizontalAlignment","center", ...
    "VerticalAlignment","middle", ...
    "FontSize",14, ...
    "FontWeight","bold", ...
    "Color","k", ...
    "BackgroundColor",[0.85 0.2 0.2], ...
    "EdgeColor","none", ...
    "Visible","off");

hStim2 = annotation(fig, "textbox", [0.73 0.83 0.2 0.05], ...
    "String", "Stim2", ...
    "FitBoxToText","off", ...
    "HorizontalAlignment","center", ...
    "VerticalAlignment","middle", ...
    "FontSize",14, ...
    "FontWeight","bold", ...
    "Color","k", ...
    "BackgroundColor",[0.85 0.85 0.2], ...
    "EdgeColor","none", ...
    "Visible","off");

hRwd = annotation(fig, "textbox", [0.73 0.93 0.2 0.05], ...
    "String", "RWD", ...
    "FitBoxToText","off", ...
    "HorizontalAlignment","center", ...
    "VerticalAlignment","middle", ...
    "FontSize",14, ...
    "FontWeight","bold", ...
    "Color","k", ...
    "BackgroundColor",[0.2 0.85 0.2], ...
    "EdgeColor","none", ...
    "Visible","off");

cueDuration = 1;  % seconds

% 6) VideoWriter
vw = VideoWriter(outFile,"MPEG-4");
vw.FrameRate = fps.vid;
vw.Quality   = 50;  % trade off quality vs file size
open(vw);

% 7) Render loop
nF = numel(tVid);
fprintf("Rendering %d frames to %s\n", round(nF/(fps.fuse/fps.vid)), outFile);

for k = 1:(fps.fuse/fps.vid):nF
    % Apply rotation
    M = eye(4);
    M(1:3,1:3) = R(:,:,k);
    Tnode.Matrix = M;

    % Event cue logic: on if within [event, event+2) for any event
    tNow = tVid(k);
    if ~isempty(eventTsec.rwd)
        cueOn1 = any(tNow >= eventTsec.rwd & tNow < (eventTsec.rwd + cueDuration));
        hRwd.Visible = matlab.lang.OnOffSwitchState(cueOn1);
    end
    if ~isempty(eventTsec.stimOn1)
        cueOn2 = any(tNow >= eventTsec.stimOn1 & tNow < (eventTsec.stimOn1 + cueDuration));
        hStim1.Visible = matlab.lang.OnOffSwitchState(cueOn2);
    end
    if ~isempty(eventTsec.stimOn2)
        cueOn3 = any(tNow >= eventTsec.stimOn2 & tNow < (eventTsec.stimOn2 + cueDuration));
        hStim2.Visible = matlab.lang.OnOffSwitchState(cueOn3);
    end
    
    title(ax, sprintf("t = %.3f s", tNow));
    drawnow limitrate

    fr = getframe(fig);
    writeVideo(vw, fr);
end

fprintf("Done.\n");
close(vw);
close(fig);

end