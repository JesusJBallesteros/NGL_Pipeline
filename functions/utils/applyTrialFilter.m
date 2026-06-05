function mask = applyTrialFilter(condition, trial2plot)
% applyTrialFilter  Build a logical row mask for trial selection from a
%                   `condition` struct and a `param.trial2plot` rule.
%
% PURPOSE:
%   Centralises the trial-selection rule that used to live inside
%   calculate_fireRate_general's main loop. Now that the firing-rate
%   functions preserve the FULL trial axis, each consumer that wants
%   a subset calls this helper to get the same row mask used historically.
%
% USAGE:
%   mask = applyTrialFilter(condition, 'allInitiated')   % drop aborted only
%   mask = applyTrialFilter(condition, 'correct')         % keep correct only
%   mask = applyTrialFilter(condition, 'incorrect')       % etc.
%
% INPUTS:
%   condition  - per-trial condition struct. Must contain at minimum a
%                .aborted field (used for the 'allInitiated' rule). For
%                any other trial2plot value, the named field must exist.
%   trial2plot - char selecting the rule:
%                  'allInitiated' -> keep trials where .aborted is false.
%                  any other      -> keep trials where condition.<trial2plot>
%                                    is logical-true.
%
% OUTPUT:
%   mask - logical row vector of length Ntotal (= numel(condition.aborted)
%          or numel(condition.<trial2plot>) for other rules).
%
% NOTES:
%   - Returns a ROW vector regardless of how the source field is shaped.
%     Use mask(:) on the caller side if you need a column.
%
% SEE ALSO:
%   calculate_fireRate_general, calculate_fireRate_byBlock,
%   plot_fireRate_session, calculate_neural_pca,
%   calculate_neural_trialEmbedding.
%
% Last modified 02.06.2026 (Jesus)

assert(isstruct(condition), 'NGL:applyTrialFilter', ...
    'condition must be a struct.');
assert(ischar(trial2plot) && ~isempty(trial2plot), 'NGL:applyTrialFilter', ...
    'trial2plot must be a non-empty char.');

if strcmpi(trial2plot, 'allInitiated')
    assert(isfield(condition,'aborted'), 'NGL:applyTrialFilter', ...
        ['param.trial2plot=''allInitiated'' but condition has no ', ...
         '.aborted field.']);
    mask = ~logical(condition.aborted(:))';
else
    assert(isfield(condition, trial2plot), 'NGL:applyTrialFilter', ...
        'condition has no field ''%s'' (param.trial2plot).', trial2plot);
    mask = logical(condition.(trial2plot)(:))';
end
end
