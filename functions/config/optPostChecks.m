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
%   - opt.gaze.downsampleStep XOR opt.gaze.targetFps (only one gets set)
%   - opt.gaze.maxFrames      XOR opt.gaze.maxSeconds
%   - opt.gaze.startTime < opt.gaze.endTime (when both set)
%   - opt.gaze.previewFrame set -> crf/preset ignored (info warning)
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

    % --- Gaze (NGL06_videoAnalysis) cross-field guards ---------------
    if isfield(opt, 'gaze')
        g = opt.gaze;
        if isfield(g,'targetFps') && ~isempty(g.targetFps) ...
                && isfield(g,'downsampleStep') && ~isempty(g.downsampleStep) ...
                && g.downsampleStep ~= 2   % 2 is the schema default; treat as unset
            error('NGL:invalidOption', ...
                ['opt.gaze.downsampleStep and opt.gaze.targetFps both set ', ...
                 '(%g and %g). Set exactly one: downsampleStep for fixed integer ', ...
                 'decimation, targetFps to derive it from fps/step.'], ...
                g.downsampleStep, g.targetFps);
        end
        if isfield(g,'maxFrames') && ~isempty(g.maxFrames) ...
                && isfield(g,'maxSeconds') && ~isempty(g.maxSeconds)
            error('NGL:invalidOption', ...
                ['opt.gaze.maxFrames and opt.gaze.maxSeconds are mutually ', ...
                 'exclusive (both were set to %g and %g). Set at most one.'], ...
                g.maxFrames, g.maxSeconds);
        end
        if isfield(g,'startTime') && ~isempty(g.startTime) ...
                && isfield(g,'endTime') && ~isempty(g.endTime) ...
                && g.startTime >= g.endTime
            error('NGL:invalidOption', ...
                ['opt.gaze.startTime (%g s) must be < opt.gaze.endTime (%g s).'], ...
                g.startTime, g.endTime);
        end
        if isfield(g,'previewFrame') && ~isempty(g.previewFrame)
            warning('NGL:gazePreviewIgnored', ...
                ['opt.gaze.previewFrame is set (%g); output will be a single ', ...
                 'PNG. Video-only options (crf=%g, preset=''%s'') are ignored.'], ...
                g.previewFrame, ...
                getfieldOr(g, 'crf', NaN), getfieldOr(g, 'preset', '<unset>'));
        end
    end
end

function v = getfieldOr(s, f, dflt)
    if isfield(s, f), v = s.(f); else, v = dflt; end
end
