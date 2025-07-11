% function [ev, numSmp, dIn] = readEvents(varargin)
function [ev, dIn] = readEvents(opt)
% Use this function to read event-codes saved in Intan (one file per channel).
%   Current version looks for the first change in a digitla pin, then
%   averages over 2 samples 0.5 ms later to make sure that all pins are set.
%   This offset is later removed and time of the first change is returned.
%   Unit of time is sample index (not seconds). 
%
% INPUTS-OPTIONAL
%  * path           : Data path (if not provided gui will open)
%  * smpRate        : sampling rate (if not provided read from header-file)
%
% OUTPUTS
%  * ev.eventNr     : Vector containing all event-codes.
%  * ev.smpInd      : Index of sample corresponding to each event-code
%  * numSmp         : total number of samples (use to chekc consistency)
%  * dIn            : Raw digital in
%
% EXAMPLES
% 1.
%
% See also: readingExample, readAux, readNeural, readAdc,
%           readHeader, readAmpToKlusta, generatePrm

% VERSION HISTORY:
% Author:        Jonas Rose
% Version:       1.1
% Last Change:
% 08.05.2016, Jonas: Release version
% 17.05.2016, Jonas: sampling rate is picked up from header file or input
% 17.05.2016, Jonas: bugfix, read events as uint16
% 05.05.2022, Aylin: corrected code for new Intan System
% 25.05.2022, Aylin: corrected code for standard and extra event codes
% 11.10.2022, Jesus: Testing old format compatibility (Does not affect new format)
% 19.08.2024, Jesus: Removed numSmp output. Modified input system

%% Defaults
pth     = opt.PathRaw; % folder for reading events
smpRate = 30000; % sampling rate
smpSum  = 2;        % an event code is sums over n samples (to catch instabilities)
skipDur = .5;       % look 0.5 ms after the first change to assume all pins are there.    

%% initializations
% PINS can be obtained from files
INfiles = string(ls(fullfile(pth,"board-DIGITAL-IN*")));
pins = 1:numel(INfiles);

% number of samples to skip - convert from ms to smpInd
% (look for first pin that changes then read the code n samples later)
smpSkip = smpRate*skipDur/1000;

%% get optional inputs
% for downwards-compatibility
% if nargin==1
%     pth = varargin{1};
% else
%     i=1;
%     while i<=length(varargin)
%         switch lower(varargin{i})
%             case 'path'
%                 i = i+1;
%                 pth         = varargin{i};
%             case 'smprate'
%                 i = i+1;
%                 smpRate     = varargin{i};
%             case 'pins'
%                 i = i+1;
%                 pins = varargin{i};
%         end
%         i=i+1;
%     end
% end
% 
% % let the user select the path
% if isempty(pth)
%     pth = uigetdir;
% end
% read sampling rate from header file
% if isnan(smpRate)
%     nfo = readHeader('path',pth,'verbose',0,'noTime');
%     smpRate = nfo.ampSmpRate;
% end

%% read all digital IN
% read the first digital channel
% if isfile('board-DIN-00.dat') % If is old data, 
%     oldsystem = 1;            % Tag it
%     fid = fopen(fullfile(pth,'board-DIN-00.dat')); % Use the old 00 as first ch
% else
fid = fopen(fullfile(pth,'board-DIGITAL-IN-01.dat')); % Corrected for new Intan System, previous: fopen(fullfile(pth,'board-DIN-00.dat'));
% end
tmp = fread(fid,inf,'uint16');
fclose(fid);
% number of samples
numSmp = size(tmp,1);
% if oldsystem
%     pins = 1:15; % Redo pins to 0:1:15
% end
% preallocate for every pin/ channel
% if oldsystem
%    dIn = zeros(numSmp,15,'uint8'); % fix old system to 16
% else
dIn = zeros(numSmp, length(pins), 'uint8');   % Corrected for new Intan System, previous: zeros(numSmp,16,'uint8');
% end

dIn(:,1) = tmp;                            % Corrected for new Intan System
clear tmp;

% read the remaining pins
for i= pins(2:end)
%     if oldsystem % here it will go from 1:15
%         fid = fopen(fullfile(pth,['board-DIN-' sprintf('%02d',i) '.dat']));
%         dIn(:,i+1) = fread(fid,numSmp,'uint16');
%         fclose(fid);
%     else % here it will go from 2:16
    fid = fopen(fullfile(pth,['board-DIGITAL-IN-' sprintf('%02d',i) '.dat'])); % Corrected for new Intan System, previous: fopen(fullfile(pth,['board-DIN-' sprintf('%02d',i) '.dat']));
    dIn(:,i) = fread(fid,numSmp,'uint16');
    fclose(fid);    
%     end
end

%% convert to sample # and event-code
% remove samples during which pin 1 and 2 are up (per default before reset of all pins,
% necessary in order to correctly extract all relevant events) 
% use extra function to set all pins to 0 before recording to make redundant: 'resetPinsIntan'
pinsOff = find(dIn(:,1)==[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],1);  % first sample with pin 1 off
if pinsOff ~= 1 % initially not all pins turned off
    % remove initial samples
    dIn(1:pinsOff-1,:) = [];
end

% find samples during which any pin is up
tmp         = find(any(dIn,2));
% the start of event codes (1ms after the first pin up)
ev.smpInd   = tmp(diff([-2; tmp])>2)+smpSkip;

%%%%%% 1. CONVERT standard events - only read out pin 1 - 4
% convert each event code to an integer
ev.eventStandard = zeros(length(ev.smpInd),1);
for i=1:length(ev.smpInd)
    % convert to decimal, sum over n samples to catch inconsitencies
    ev.eventStandard(i) = binvec2dec(sum(...
        dIn(ev.smpInd(i):ev.smpInd(i)+smpSum,1:4)...
        ));
end

%%%%%% 2. CONVERT extra events - read out 5 - 16 (set pin 1 - 4 to 0)
% convert each event code to an integer
ev.eventExtra = zeros(length(ev.smpInd),1);
dIn(:,1:4) = 0; % set first 4 pins to 0
for i=1:length(ev.smpInd)
    % convert to decimal, sum over n samples to catch inconsitencies
    ev.eventExtra(i) = binvec2dec(sum(...
        dIn(ev.smpInd(i):ev.smpInd(i)+smpSum,1:end)...
        ));
end

% remove the delay caused by skipping samples
ev.smpInd = ev.smpInd-smpSkip;