function aggregated = loadAggregatedSpikes(input, area)
% loadAggregatedSpikes  Load one per-area NGL03 aggregated output.
%
% PURPOSE:
%   Returns ONE per-area aggregated struct used by every cross-subject
%   script (NGL04_fireRate, NGL04_PCA, ...). Each call loads exactly
%   one area's file:
%       <input.analysis>/aggregated_<area>.mat                 (study-level)
%     fallback (stitched from per-subject files):
%       <input.analysis>/<subj>/<subj>_aggregated_<area>.mat
%
%   The returned struct is FLAT — no per-area sub-structs. NGL04 calls
%   this once per area in its area loop.
%
% USAGE:
%   aggregated = loadAggregatedSpikes(input);            % single-area / 'main'
%   aggregated = loadAggregatedSpikes(input, 'NCL');     % per-area
%   aggregated = loadAggregatedSpikes(input, 'STR');
%
% INPUTS:
%   input - NGL pipeline input struct, fields used:
%             .analysis     <studyName>/analysisCode/<analysis>
%             .nsubjects
%             .subjects(x).name
%   area  - char, area tag. Empty / omitted -> 'main' (single-area
%           sentinel; matches NGL03_aggregate naming).
%
% OUTPUT:
%   aggregated - struct with cell arrays sized (nSubj x maxSess):
%             .allspike      per-(subj,sess) FLAT spike struct
%             .allneurons    per-(subj,sess) FLAT neurons struct
%             .allcondition  per-(subj,sess) condition struct
%             (other fields preserved verbatim from the source files)
%
% CONTRACT:
%   Empty cells are valid for sessions a given subject doesn't have.
%   Callers should guard with isempty/isstruct before indexing.
%
% LEGACY HANDLING:
%   If the requested per-area file is missing AND a legacy
%   aggregated.mat (or <subj>_aggregated.mat) sits in the analysis tree,
%   raise NGL:loadAggregatedSpikes:legacyFile with a pointer to
%   migrate_aggregated_to_perArea.m. Silently falling back would
%   re-introduce the nested-shape bug that this refactor exists to
%   eliminate.
%
% Last modified 26.06.2026 (Jesus) - per-area file layout; legacy guard.

    if nargin < 2 || isempty(area), area = 'main'; end
    area = char(area);

    studyFile  = fullfile(input.analysis, ['aggregated_' area '.mat']);
    legacyFile = fullfile(input.analysis, 'aggregated.mat');

    if isfile(studyFile)
        fprintf('loadAggregatedSpikes: study-level %s\n', studyFile);
        aggregated = load(studyFile);
        return
    end

    if isfile(legacyFile)
        error('NGL:loadAggregatedSpikes:legacyFile', ...
            ['Found legacy nested aggregated.mat at\n   %s\n', ...
             'but the per-area file\n   %s\n', ...
             'does not exist. NGL04 now requires per-area files. ', ...
             'Run migrate_aggregated_to_perArea.m once to split the ', ...
             'legacy file into per-area copies, then re-run NGL04.'], ...
            legacyFile, studyFile);
    end

    fprintf(['loadAggregatedSpikes: study-level aggregated_%s.mat not ', ...
             'found; stitching per-subject files.\n'], area);
    aggregated = struct();
    anyHit     = false;
    for x = 1:input.nsubjects
        sname      = input.subjects(x).name;
        subjFile   = fullfile(input.analysis, sname, ...
                              [sname '_aggregated_' area '.mat']);
        legacySubj = fullfile(input.analysis, sname, ...
                              [sname '_aggregated.mat']);
        if ~isfile(subjFile)
            if isfile(legacySubj)
                error('NGL:loadAggregatedSpikes:legacyFile', ...
                    ['Found legacy per-subject\n   %s\n', ...
                     'but per-area file\n   %s\n', ...
                     'does not exist. Run migrate_aggregated_to_perArea.m.'], ...
                    legacySubj, subjFile);
            end
            warning('NGL:loadAggregatedSpikes:missingSubj', ...
                'Per-subject aggregated file missing for %s / %s: %s', ...
                sname, area, subjFile);
            continue
        end
        anyHit = true;
        S  = load(subjFile);
        fn = fieldnames(S);
        for k = 1:numel(fn)
            field = fn{k};
            row   = S.(field);     % {1 x nSess_for_this_subj}
            if ~iscell(row), continue, end
            if ~isfield(aggregated, field), aggregated.(field) = {}; end
            for y = 1:numel(row)
                aggregated.(field){x, y} = row{y};
            end
        end
    end

    if ~anyHit
        warning('NGL:loadAggregatedSpikes:nothing', ...
            ['No per-subject aggregated_%s.mat files found anywhere ', ...
             'under %s; returned aggregated struct is empty.'], ...
            area, input.analysis);
    end
end
