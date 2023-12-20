function intan2FieldTrip(sessions, varargin)
% Merged 'intan2mat_wrapper' and 'mat2FieldTrip' functions for a real
% metawrapper.

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

% Default options if not specified
if ~isfield(opt,'lowpass'), opt.lowpass = [  0  400];  end

% Proceed with the main functions
INTANdata = intan2MAT_wrapper(sessions, opt);

MAT2FieldTrip(INTANdata, opt);

end