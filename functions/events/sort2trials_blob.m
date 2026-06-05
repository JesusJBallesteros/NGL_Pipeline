function neurons = sort2trials_blob(spike, blobMerges, opt) %#ok<INUSD>
% sort2trials_blob  Bin per-cluster spike times around social-arena
%                    blob-interaction events.
%
% PURPOSE:
%   Companion to sort2trials.m for the social-arena paradigm. Where
%   sort2trials operates on a (2 x nAlignments) cell trialdef from
%   trialdefGen, this function operates on a numeric [nInteractions x 2+]
%   matrix produced by the offline video tracker (typically blob.Merges).
%   It uses a fixed +/- 5 s padding window around each interaction's
%   start/end times.
%
%   Output shape mirrors the "raw" social branch that used to live as
%   the else-branch of sort2trials.m: a flat cell-of-cells, with no
%   alignment-name keying.
%
% USAGE:
%   neurons = sort2trials_blob(spike, blob.Merges, opt)
%
% INPUTS:
%   spike      - NGL spike struct from loadSpikes. Required: .timestamp
%                (cell-per-cluster, SECONDS).
%   blobMerges - [nInteractions x >=2] numeric matrix:
%                  column 1 = interaction start (seconds since session start)
%                  column 2 = interaction end   (seconds since session start)
%   opt        - resolved options struct (unused today; kept on the
%                signature for symmetry with sort2trials and to leave
%                room for opt.socialPadSec when it becomes a knob).
%
% OUTPUT:
%   neurons    - {nClust x 1} cell-of-cells:
%                  neurons{c}{i} = vector of ms timestamps for cluster c,
%                                  interaction i, relativised to the
%                                  interaction start.
%
% KNOWN ISSUES:
%   - The +/-5 s padding window is hardcoded. Should be moved to an opt
%     field (e.g. opt.socialPadSec) when this branch gets exercised
%     more seriously. See audit item R.
%   - No ROI propagation. spike.roi is not copied onto a corresponding
%     neurons field today; downstream callers that need it should read
%     from spike directly.
%
% SEE ALSO:
%   sort2trials (standard alignment-keyed path).
%
% Last modified 02.06.2026 (Jesus) - split from sort2trials (#10 R)

assert(isnumeric(blobMerges) && size(blobMerges,2) >= 2, ...
    'NGL:sort2trials_blob:badBlobMerges', ...
    'blobMerges must be a numeric matrix with at least 2 columns [start end] in seconds.');

nclus  = length(spike.label);
ntrial = size(blobMerges, 1);

neurons = cell(nclus, 1);

PAD_SEC = 5;  % fixed +/- padding around each interaction; see KNOWN ISSUES.

for c = 1:nclus
    neurons{c, 1} = cell(ntrial, 1);
    st = spike.timestamp{1,c} * 1000;  % seconds -> ms

    for i = 1:ntrial
        startMs = (blobMerges(i,1) - PAD_SEC) * 1000;
        endMs   = (blobMerges(i,2) + PAD_SEC) * 1000;
        idx     = (st >= startMs) & (st < endMs);

        % Relativise to interaction start (ms).
        neurons{c,1}{i,1} = st(idx) - blobMerges(i,1) * 1000;
    end
end
end
