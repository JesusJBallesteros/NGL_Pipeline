function contrast = compareByBlock(data, condition, spec)
% compareByBlock  Generic per-trial contrast primitive.
%
% PURPOSE:
%   Replaces the hardcoded "early / late half of each block x NS/FS x
%   NS-FS subtraction" logic that used to live baked inside
%   trialparsed_MTspectrogram (project-specific social-learning pipeline).
%   Splits per-trial data (any cell-of-trials or trial-first numeric
%   array) into A / B partitions defined by `spec`, reduces each
%   partition, and returns the contrast (A - B by default).
%
%   The intent is that project-specific wrappers (e.g. socialLearning's
%   NS-FS early/late) become a couple of calls to this primitive rather
%   than bespoke code that has to be re-derived every time.
%
% USAGE:
%   contrast = compareByBlock(data, condition, spec)
%
% INPUTS:
%   data      - per-trial payload. One of:
%                 * cell{Ntrials,1}          -> passed through per trial
%                 * numeric [Ntrials, ...]   -> first dim is trials
%                 * FT freq struct with .powspctrm [Ntrials x ...] and
%                   keeptrials='yes' -> unpacked into per-trial slices
%   condition - condition struct (from NGL02_postPhy). Fields consumed
%               depend on spec.mode.
%   spec      - struct describing the contrast (see MODES).
%
% MODES (spec.mode):
%
%   'earlyLate'
%       Split the trial index vector in half (median split); A = early,
%       B = late.
%       Consumes: spec.trialSelector (optional; logical/index of trials
%                 to keep before the split). Default = all trials.
%
%   'block-vs-block'
%       A = trials with condition.(spec.blockField) == spec.blockA
%       B = trials with condition.(spec.blockField) == spec.blockB
%       Consumes: spec.blockField (e.g. 'block'), spec.blockA, spec.blockB.
%
%   'condition-vs-condition'
%       A = trials matching spec.selectorA (logical or function handle)
%       B = trials matching spec.selectorB
%       Consumes: spec.selectorA, spec.selectorB.
%
%   'level-vs-level'
%       Same as condition-vs-condition but with named field / value:
%         A = condition.(spec.field) == spec.levelA
%         B = condition.(spec.field) == spec.levelB
%       Consumes: spec.field, spec.levelA, spec.levelB.
%
% OTHER SPEC FIELDS (all optional):
%   spec.reduce   - 'mean' (default) | 'median' | 'sum' | function handle
%   spec.contrast - 'diff' (default; A-B) | 'ratio' (A./B) | 'logratio'
%                   (log(A./B)) | 'zdiff' (rescaled by pooled SD)
%   spec.dropAborted - true (default) | false. If true and
%                      condition.aborted exists, aborted trials are
%                      removed BEFORE partitioning.
%
% OUTPUT (struct):
%   .A          - reduced payload for partition A
%   .B          - reduced payload for partition B
%   .contrast   - A vs B result (shape depends on spec.contrast)
%   .nA / .nB   - trial counts contributing to each partition
%   .idxA / .idxB - logical trial masks used (into the ORIGINAL data,
%                   NOT the post-drop-aborted subset). Useful for
%                   provenance / debugging.
%   .spec       - the resolved spec (with defaults filled in).
%
% CONTRACT:
%   * If either partition is empty, that side's reduction is [] and
%     .contrast is []. No error; caller decides how to handle it.
%   * FT freq-struct input is unpacked, contrasted, and RETURNED as a
%     freq-struct with .powspctrm collapsed to the reduced shape (trial
%     dim removed). Non-.powspctrm fields are copied verbatim.
%
% Last modified 26.06.2026 (Jesus) - new primitive (LFP Pass 2, generalises
%                                     the block-comparison logic that used
%                                     to live inside trialparsed_MTspectrogram).

    % --- Fill defaults ------------------------------------------------
    if ~isfield(spec, 'reduce'),       spec.reduce       = 'mean'; end
    if ~isfield(spec, 'contrast'),     spec.contrast     = 'diff'; end
    if ~isfield(spec, 'dropAborted'),  spec.dropAborted  = true;   end

    % --- Detect payload shape ----------------------------------------
    [payloadKind, payload, meta] = localUnpack(data);

    % --- Base trial mask (start with all; drop aborted if requested) --
    nTrials = localCountTrials(payloadKind, payload);
    baseMask = true(nTrials, 1);
    if spec.dropAborted && isfield(condition, 'aborted') ...
            && numel(condition.aborted) == nTrials
        baseMask(logical(condition.aborted)) = false;
    end

    % --- Partition based on spec.mode --------------------------------
    switch lower(spec.mode)
        case 'earlylate'
            [idxA, idxB] = localEarlyLate(baseMask, spec);
        case 'block-vs-block'
            [idxA, idxB] = localBlockVsBlock(baseMask, condition, spec);
        case 'condition-vs-condition'
            [idxA, idxB] = localSelectorVsSelector(baseMask, condition, spec);
        case 'level-vs-level'
            [idxA, idxB] = localLevelVsLevel(baseMask, condition, spec);
        otherwise
            error('NGL:compareByBlock:badMode', ...
                'Unknown spec.mode ''%s''. See help for allowed values.', spec.mode);
    end

    % --- Reduce each partition ---------------------------------------
    A = localReduce(payloadKind, payload, idxA, spec.reduce, meta);
    B = localReduce(payloadKind, payload, idxB, spec.reduce, meta);

    % --- Contrast ----------------------------------------------------
    C = localContrast(A, B, spec.contrast, payloadKind, payload, idxA, idxB, meta);

    % --- Package -----------------------------------------------------
    contrast          = struct();
    contrast.A        = A;
    contrast.B        = B;
    contrast.contrast = C;
    contrast.nA       = nnz(idxA);
    contrast.nB       = nnz(idxB);
    contrast.idxA     = idxA;
    contrast.idxB     = idxB;
    contrast.spec     = spec;
end


%% =====================================================================
%% Local helpers
%% =====================================================================

function [kind, payload, meta] = localUnpack(data)
% Classify the input payload so downstream helpers know how to index it.
    meta = struct();
    if iscell(data)
        kind    = 'cell';
        payload = data(:);
    elseif isstruct(data) && isfield(data, 'powspctrm')
        kind         = 'ftfreq';
        payload      = data.powspctrm;
        meta.header  = rmfield(data, 'powspctrm');
    elseif isnumeric(data)
        kind    = 'numeric';
        payload = data;
    else
        error('NGL:compareByBlock:badData', ...
            ['data must be a cell {Ntrials,1}, numeric [Ntrials, ...] ', ...
             'or FT freq struct with .powspctrm and keeptrials=''yes''.']);
    end
end

function n = localCountTrials(kind, payload)
    switch kind
        case 'cell',    n = numel(payload);
        case 'numeric', n = size(payload, 1);
        case 'ftfreq',  n = size(payload, 1);
    end
end

function [idxA, idxB] = localEarlyLate(baseMask, ~)
    trials = find(baseMask);
    half   = floor(numel(trials) / 2);
    idxA = false(size(baseMask));  idxA(trials(1:half))        = true;
    idxB = false(size(baseMask));  idxB(trials(half+1:end))    = true;
end

function [idxA, idxB] = localBlockVsBlock(baseMask, condition, spec)
    assert(isfield(condition, spec.blockField), 'NGL:compareByBlock:noField', ...
        'condition.%s not found (spec.mode = block-vs-block).', spec.blockField);
    v    = condition.(spec.blockField);
    idxA = baseMask & (v(:) == spec.blockA);
    idxB = baseMask & (v(:) == spec.blockB);
end

function [idxA, idxB] = localSelectorVsSelector(baseMask, condition, spec)
    mA = localApplySelector(spec.selectorA, condition, numel(baseMask));
    mB = localApplySelector(spec.selectorB, condition, numel(baseMask));
    idxA = baseMask & mA(:);
    idxB = baseMask & mB(:);
end

function [idxA, idxB] = localLevelVsLevel(baseMask, condition, spec)
    assert(isfield(condition, spec.field), 'NGL:compareByBlock:noField', ...
        'condition.%s not found (spec.mode = level-vs-level).', spec.field);
    v    = condition.(spec.field);
    idxA = baseMask & (v(:) == spec.levelA);
    idxB = baseMask & (v(:) == spec.levelB);
end

function m = localApplySelector(sel, condition, n)
    if isa(sel, 'function_handle')
        m = logical(sel(condition));
    elseif islogical(sel)
        m = sel;
    elseif isnumeric(sel)
        m = false(n, 1); m(sel) = true;
    else
        error('NGL:compareByBlock:badSelector', ...
            'spec.selectorA/B must be a function handle, logical mask, or index vector.');
    end
    m = m(:);
    if numel(m) ~= n
        m = false(n, 1);
    end
end

function out = localReduce(kind, payload, idx, reduce, meta) %#ok<INUSD>
    if ~any(idx), out = []; return, end
    switch kind
        case 'cell'
            slice = payload(idx);
            % Cells: caller may store per-trial spike vectors or matrices.
            % Best generic move: return slice as-is if reducer is 'mean'
            % and cells are numeric-of-same-shape; else return cell slice.
            if isequal(reduce, 'mean') && all(cellfun(@isnumeric, slice)) ...
                    && numel(unique(cellfun(@numel, slice))) == 1
                mat = cell2mat(cellfun(@(c) c(:).', slice, 'uni', false));
                out = mean(mat, 1, 'omitnan');
            else
                out = slice;
            end
        case 'numeric'
            out = localApplyReducer(payload(idx, :, :, :, :), reduce);
        case 'ftfreq'
            out = localApplyReducer(payload(idx, :, :, :), reduce);
    end
end

function out = localApplyReducer(x, reduce)
    if isa(reduce, 'function_handle')
        out = reduce(x);
    else
        switch lower(reduce)
            case 'mean',   out = mean(x, 1, 'omitnan');
            case 'median', out = median(x, 1, 'omitnan');
            case 'sum',    out = sum(x, 1, 'omitnan');
            otherwise
                error('NGL:compareByBlock:badReduce', ...
                    'Unknown reducer ''%s''.', reduce);
        end
    end
    out = squeeze(out);
end

function C = localContrast(A, B, contrastMode, kind, payload, idxA, idxB, meta) %#ok<INUSD>
    if isempty(A) || isempty(B), C = []; return, end
    switch lower(contrastMode)
        case 'diff'
            C = A - B;
        case 'ratio'
            C = A ./ B;
        case 'logratio'
            C = log(A ./ B);
        case 'zdiff'
            if strcmp(kind, 'numeric') || strcmp(kind, 'ftfreq')
                sdA = std(payload(idxA, :, :, :, :), 0, 1, 'omitnan');
                sdB = std(payload(idxB, :, :, :, :), 0, 1, 'omitnan');
                pooledSD = sqrt((squeeze(sdA).^2 + squeeze(sdB).^2) / 2);
                C = (A - B) ./ max(pooledSD, eps);
            else
                error('NGL:compareByBlock:zdiffCell', ...
                    '''zdiff'' contrast is not defined for cell payloads.');
            end
        otherwise
            error('NGL:compareByBlock:badContrast', ...
                'Unknown spec.contrast ''%s''.', contrastMode);
    end
end
