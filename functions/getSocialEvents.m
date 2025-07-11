function [socialidx] = getSocialEvents(opt)

% Create output variable
socialidx = struct();

% Read excel table and convert to array
T = readtable(fullfile(opt.behavFiles, "VidAssessment.xlsx"));
T = table2array(T(:,3:end));

% Check no missing cells. If any, assume NaNs are 0s.
nanidx = isnan(T);
if any(any(nanidx)), T(nanidx) = 0; end

% Now get indexes into struct, ready to output
socialidx.tuteeLEDview   = logical(T(:,1));
socialidx.tuteeLEDreact  = logical(T(:,2));
socialidx.tuteeSTview    = logical(T(:,3));
socialidx.tuteeSTreact   = logical(T(:,4));
socialidx.tuteeSTpeak    = logical(T(:,5));
socialidx.tuteeRWview    = logical(T(:,6));
socialidx.tuteeRWreact   = logical(T(:,7));
socialidx.tutorPresent   = logical(T(:,8));
socialidx.tutorLEDview   = logical(T(:,9));
socialidx.tutorLEDreact  = logical(T(:,10));
socialidx.tutorSTview    = logical(T(:,11));
socialidx.tutorSTreact   = logical(T(:,12));
socialidx.tutorRWview    = logical(T(:,13));
socialidx.tutorRWreact   = logical(T(:,14));
socialidx.tutee_tutorObs = logical(T(:,15));
socialidx.tutee_tutorInt = logical(T(:,16));

end