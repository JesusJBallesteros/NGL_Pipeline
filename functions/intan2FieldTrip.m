function intan2FieldTrip(sessions, varargin)
% Merged 'intan2mat_wrapper' and 'mat2FieldTrip' functions for a real
% metawrapper.

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

% Default options if not specified
if ~isfield(opt,'lowpass'), opt.lowpass = [  0  400];  end

% Proceed with the main fuinctions
INTANdata = intan2mat_wrapper(sessions, opt);
mat2FieldTrip(INTANdata, opt);

end