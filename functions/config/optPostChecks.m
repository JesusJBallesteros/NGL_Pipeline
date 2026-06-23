function optPostChecks(opt)
% optPostChecks  Cross-field validation rules that the schema can't
%                express (each rule depends on more than one option).
%
% PURPOSE:
%   The schema's per-entry validators see one leaf at a time. Anything
%   that depends on TWO fields lives here. Throws NGL:invalidOption on
%   the first violation (same identifier as the schema validator and
%   the legacy inline asserts).
%
% RULES:
%   - opt.stepSz_ms <= opt.binSize_ms
%     (canonical FR bins must overlap or be contiguous)
%   - opt.fireRatePlot.stepSz_ms <= opt.fireRatePlot.binSize_ms
%   - opt.pcaPlot.stepSz_ms      <= opt.pcaPlot.binSize_ms
%   - opt.aggregateSubjects=true requires opt.aggregateSessions=true
%     (subjects only aggregate after sessions).
%
% NOT ENFORCED (intentionally preserved as legacy behaviour):
%   - opt.popDyn.alignIdx <= numel(opt.alignto). The legacy set_default
%     commented this out so the new auto-iterate PCA path is free to
%     ignore alignIdx. If you ever revive jPCA/GPFA placeholders that
%     read alignIdx, add the check back here.
%
% Last modified 09.06.2026 (Jesus)

    if opt.stepSz_ms > opt.binSize_ms
        error('NGL:invalidOption', ...
            'opt.stepSz_ms (%g) must be <= opt.binSize_ms (%g); otherwise bins do not overlap.', ...
            opt.stepSz_ms, opt.binSize_ms);
    end

    if opt.fireRatePlot.stepSz_ms > opt.fireRatePlot.binSize_ms
        error('NGL:invalidOption', ...
            'opt.fireRatePlot.stepSz_ms (%g) must be <= opt.fireRatePlot.binSize_ms (%g).', ...
            opt.fireRatePlot.stepSz_ms, opt.fireRatePlot.binSize_ms);
    end

    if opt.pcaPlot.stepSz_ms > opt.pcaPlot.binSize_ms
        error('NGL:invalidOption', ...
            'opt.pcaPlot.stepSz_ms (%g) must be <= opt.pcaPlot.binSize_ms (%g).', ...
            opt.pcaPlot.stepSz_ms, opt.pcaPlot.binSize_ms);
    end

    if opt.aggregateSubjects && ~opt.aggregateSessions
        error('NGL:invalidOption', ...
            ['opt.aggregateSubjects=true requires opt.aggregateSessions=true ', ...
             '(subjects only aggregate after sessions).']);
    end

    if isfield(opt,'regenFrom') && isfield(opt.regenFrom,'preproc') && opt.regenFrom.preproc
        if isempty(opt.regenFrom.nChannels)
            assert(isnumeric(opt.numChannels) && isscalar(opt.numChannels) && opt.numChannels > 0, ...
                'NGL:invalidOption', ...
                ['opt.regenFrom.preproc=true requires opt.numChannels to be a positive scalar ', ...
                 '(used as the fallback channel count when opt.regenFrom.nChannels is empty).']);
        end
    end
end
