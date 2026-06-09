function aggregated = loadAggregatedSpikes(input)
% loadAggregatedSpikes  Load NGL03_acrossSession output.
%
% PURPOSE:
%   Returns the aggregated struct used by every cross-subject script
%   (NGL04_fireRate, NGL04_PCA, ...). Prefers the study-level
%   aggregated.mat (built when opt.aggregateSubjects=true); falls back
%   to stitching per-subject *_aggregated.mat files into a
%   (nSubjects x maxSessions) struct of cells.
%
% USAGE:
%   aggregated = loadAggregatedSpikes(input);
%
% INPUT:
%   input - NGL pipeline input struct, fields used:
%             .analysis     <studyName>/analysisCode
%             .nsubjects
%             .subjects(x).name
%
% OUTPUT:
%   aggregated - struct with cell arrays sized (nSubj x maxSess):
%             .allspike      per-(subj,sess) spike struct (loadSpikes output)
%             .allneurons    per-(subj,sess) neurons struct (sort2trials)
%             .allcondition  per-(subj,sess) condition struct
%             (other fields preserved verbatim from the source files)
%
% CONTRACT:
%   Empty cells are valid for sessions a given subject doesn't have.
%   Callers should guard with isempty/isstruct before indexing.
%
% Last modified 09.06.2026 (Jesus)

    studyFile = fullfile(input.analysis, 'aggregated.mat');
    if isfile(studyFile)
        fprintf('loadAggregatedSpikes: study-level %s\n', studyFile);
        aggregated = load(studyFile);
        return
    end
    fprintf('loadAggregatedSpikes: study-level aggregated.mat not found; stitching per-subject files.\n');
    aggregated = struct();
    for x = 1:input.nsubjects
        sname    = input.subjects(x).name;
        subjFile = fullfile(input.analysis, sname, [sname '_aggregated.mat']);
        if ~isfile(subjFile)
            warning('NGL:loadAggregatedSpikes:missingSubj', ...
                'Per-subject aggregated file missing for %s: %s', sname, subjFile);
            continue
        end
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
end
