function ok = hasCuratedClusters(spike_area)
% hasCuratedClusters  Return true iff at least one cluster carries a
%                     Phy human-curation label of 'good' or 'mua'.
%
% PURPOSE:
%   loadSpikes can return a non-empty spike struct (cluster_info.tsv
%   exists, .label cell has entries) even when the user has not finished
%   Phy curation — e.g. only the 'noise' clusters were tagged and the
%   rest were left blank. Downstream fire-rate / population-dynamics
%   processing expects at least one 'good' or 'mua' cluster to be there;
%   running with everything-unlabelled produces empty plots at best and
%   crashes at worst.
%
%   This helper makes that "is curation effectively complete?" check
%   explicit so NGL02_postPhy can mark areas as skipped with a clear
%   reason text file (see noteSkippedArea) instead of failing deep in
%   the analysis stack.
%
% USAGE:
%   ok = hasCuratedClusters(spike.NCL)        % multi-area
%   ok = hasCuratedClusters(spike)            % single-area
%
% INPUT:
%   spike_area - struct from loadSpikes. The function inspects
%                .HumanLabel (the Phy 'group' column). If that field is
%                absent or every entry is empty / 'noise', the function
%                returns false.
%
% OUTPUT:
%   ok - logical. true iff at least one cluster in spike_area has
%        HumanLabel == 'good' or 'mua' (case-insensitive).
%
% NOTES:
%   - Project-specific custom labels are not recognised. If your study
%     uses additional categories ('su', 'fs', ...), extend the accept
%     list here. We deliberately do NOT accept HumanLabel == ''
%     (unlabelled) because that's the symptom this check is meant to
%     catch.
%
% SEE ALSO:
%   loadSpikes, noteSkippedArea, NGL02_postPhy.
%
% Last modified 23.06.2026 (Jesus)

    ok = false;
    if ~isstruct(spike_area) || isempty(spike_area), return; end
    if ~isfield(spike_area, 'HumanLabel') || isempty(spike_area.HumanLabel), return; end

    accept = {'good', 'mua'};
    hl = spike_area.HumanLabel;
    for k = 1:numel(hl)
        v = hl{k};
        if isempty(v), continue; end
        if any(strcmpi(char(string(v)), accept))
            ok = true;
            return
        end
    end
end
