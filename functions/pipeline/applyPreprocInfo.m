function opt = applyPreprocInfo(opt, preprocInfo, varargin)
% applyPreprocInfo  Overlay NGL01-owned fields from a preprocInfo snapshot
%                   onto a current opt struct.
%
% PURPOSE:
%   When NGL02_postPhy runs after a separate NGL01_Main session (potentially
%   after MATLAB has been restarted), the user's current opt may not match
%   the options that were actually used during preprocessing. This helper
%   overlays the NGL01-owned fields from the saved snapshot onto the live
%   opt, leaving NGL02-owned fields (and other unrelated knobs) intact.
%
%   Any field whose value would change is reported, so the user can spot
%   silent overrides (e.g. opt.numChannels in the SetAndRunMe disagrees
%   with what was used in NGL01).
%
% USAGE:
%   opt = applyPreprocInfo(opt, preprocInfo);
%   opt = applyPreprocInfo(opt, preprocInfo, 'Verbose', false);
%   opt = applyPreprocInfo(opt, preprocInfo, 'Fields', {'numChannels', 'alignto'});
%
% INPUTS:
%   opt          - current options struct (post-set_default).
%   preprocInfo  - struct returned by loadPreprocInfo. Must contain .opt.
%
% NAME-VALUE PAIRS:
%   'Fields'  - cell array of field names to overlay. Default: the NGL01-owned
%               field list defined below.
%   'Verbose' - logical, default true. Print one line per field whose value
%               changes (and a summary line if nothing changes).
%
% OUTPUT:
%   opt  - opt with the requested fields overwritten by preprocInfo.opt.
%          Fields absent from preprocInfo.opt are left unchanged.
%
% NGL01-OWNED FIELDS (default overlay list):
%   Acquisition:   numChannels
%   Events:        RetrieveEvents, alignto, trEvents, addtime, uselog
%   Motion:        GetMotionSensors
%   Preprocessing: noise, lowpass, lowpassFT, highpass, linefilter, CAR,
%                  dwnsmplRate, timebreak
%   Sorting:       kilosort, KSchanMapFile, bombcell, callBcGUI
%
%   NGL02-owned fields are deliberately NOT in this list:
%   doSpikething, doLFPthing, offlineTrack, FLIP, useTrack, trialparsed,
%   artifdet, spectrogram, neurDyn. Also excluded: bin, FieldTrip, doNWB,
%   phy (NGL01-only flags that have no effect at NGL02 time).
%
% SEE ALSO:
%   savePreprocInfo, loadPreprocInfo
%
% Last modified 27.05.2026 (Jesus)

%% Parse optional args
defaultFields = { ...
    'numChannels', ...
    'RetrieveEvents', 'alignto', 'trEvents', 'addtime', 'uselog', ...
    'GetMotionSensors', ...
    'noise', 'lowpass', 'lowpassFT', 'highpass', 'linefilter', 'CAR', ...
    'dwnsmplRate', 'timebreak', ...
    'kilosort', 'KSchanMapFile', 'bombcell', 'callBcGUI' ...
    };

p = inputParser;
p.addParameter('Fields',  defaultFields, @iscell);
p.addParameter('Verbose', true,          @(v) islogical(v) && isscalar(v));
p.parse(varargin{:});
fieldsToOverlay = p.Results.Fields;
verbose         = p.Results.Verbose;

%% Validate inputs
assert(isstruct(opt), 'NGL:applyPreprocInfo', 'opt must be a struct.');
assert(isstruct(preprocInfo) && isfield(preprocInfo,'opt') && isstruct(preprocInfo.opt), ...
    'NGL:applyPreprocInfo', ...
    'preprocInfo must be a struct with a .opt sub-struct (from loadPreprocInfo).');

savedOpt = preprocInfo.opt;
nChanged = 0;

%% Overlay loop
for k = 1:numel(fieldsToOverlay)
    f = fieldsToOverlay{k};
    if ~isfield(savedOpt, f)
        continue   % preprocInfo doesn't carry this field; leave opt alone
    end

    newVal = savedOpt.(f);

    if isfield(opt, f)
        oldVal = opt.(f);
        if ~isequaln(oldVal, newVal)
            nChanged = nChanged + 1;
            if verbose
                fprintf('  preprocInfo override: opt.%s  %s  ->  %s\n', ...
                        f, localShort(oldVal), localShort(newVal));
            end
        end
    else
        nChanged = nChanged + 1;
        if verbose
            fprintf('  preprocInfo added:    opt.%s = %s\n', ...
                    f, localShort(newVal));
        end
    end

    opt.(f) = newVal;
end

if verbose && nChanged == 0
    fprintf('  preprocInfo overlay: no field changes (opt already matched the snapshot).\n');
end
end

% ---------------------------------------------------------------------------
function s = localShort(v)
% Compact one-line stringifier used only for log messages.
    try
        if ischar(v)
            s = ['''' v ''''];
        elseif isstring(v)
            s = ['"' char(v) '"'];
        elseif isnumeric(v) || islogical(v)
            if isempty(v)
                s = '[]';
            elseif isscalar(v)
                s = num2str(v);
            else
                s = ['[' num2str(size(v,1)) 'x' num2str(size(v,2)) ' ' class(v) ']'];
            end
        elseif iscell(v)
            if isempty(v)
                s = '{}';
            elseif all(cellfun(@ischar, v))
                s = ['{' strjoin(cellfun(@(c) ['''' c ''''], v, 'uni', 0), ', ') '}'];
            else
                s = sprintf('{%d-cell}', numel(v));
            end
        else
            s = ['<' class(v) '>'];
        end
    catch
        s = '<unprintable>';
    end
end
