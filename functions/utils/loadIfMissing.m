function loadIfMissing(varname, filepath, varargin)
% loadIfMissing  Load a variable into the caller's workspace iff missing.
%
% PURPOSE:
%   Collapse the repeated NGL pattern:
%       if ~exist('events','var'), load(fullfile(opt.trialSorted, "events.mat")); end
%   into a single line. Also handles the legacy fallback where a saved
%   file stored the variable under a different name (e.g. 'conditions'
%   in older session caches vs the current 'condition'), and the looser
%   "the file contains exactly one variable, take it" case.
%
% USAGE:
%   loadIfMissing('events',    fullfile(opt.trialSorted, "events.mat"))
%   loadIfMissing('condition', fullfile(opt.trialSorted, "condition.mat"), 'conditions')
%
% INPUTS:
%   varname   - char, the name the variable should have in the caller.
%               Skip the load entirely if a variable of this name
%               already exists in the caller's workspace.
%   filepath  - char, absolute path to the .mat file. If the file does
%               not exist the call is a no-op (no error, no warning).
%               Callers that want hard failure on a missing file should
%               check isfile() themselves before calling.
%   varargin  - optional fallback variable names (char). If the file
%               does not contain `varname`, try each fallback in order
%               and assign the first match to `varname` in the caller.
%               If none match and the file contains exactly one variable,
%               take it.
%
% NOTES:
%   - Uses evalin('caller', ...) and assignin('caller', ...) by design.
%     The original pattern this replaces is already caller-workspace
%     based (~exist('var','var')), so a helper that lives in a function
%     necessarily has to reach into the caller. Use only in scripts
%     and at top-level analysis loops; do NOT call from deeply nested
%     functions where the "caller" is unclear.
%
% SEE ALSO:
%   load, evalin, assignin.
%
% Last modified 02.06.2026 (Jesus) - audit item L

assert(ischar(varname) && ~isempty(varname), ...
    'NGL:loadIfMissing', 'varname must be a non-empty char.');
assert(ischar(filepath) && ~isempty(filepath), ...
    'NGL:loadIfMissing', 'filepath must be a non-empty char.');

%% Already in the caller's workspace -> no-op.
if evalin('caller', sprintf('exist(''%s'',''var'')', varname))
    return
end

%% File missing -> no-op (caller is responsible for hard-error guards).
if ~isfile(filepath)
    return
end

%% Load the .mat and find the right field.
S = load(filepath);
fn = fieldnames(S);

picked = '';
if isfield(S, varname)
    picked = varname;
else
    for k = 1:numel(varargin)
        if isfield(S, varargin{k})
            picked = varargin{k};
            break
        end
    end
    if isempty(picked) && numel(fn) == 1
        picked = fn{1};
    end
end

if isempty(picked)
    warning('NGL:loadIfMissing:noMatch', ...
        ['File %s exists but contains neither ''%s'' nor any of the ', ...
         'requested fallbacks. Leaving the variable undefined in the caller.'], ...
        filepath, varname);
    return
end

assignin('caller', varname, S.(picked));
end
