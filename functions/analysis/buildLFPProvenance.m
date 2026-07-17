function prov = buildLFPProvenance(sourceFT, opt, extra)
% buildLFPProvenance  Small provenance struct for every LFP output file.
%
% PURPOSE:
%   Standardises the "what produced this .mat" breadcrumb that every LFP
%   compute (computeContinuousTFR, computeTrialparsedTFR, planned
%   spike-field / regression analyses) attaches to its output file.
%   Cheap to write, invaluable when a downstream plot looks suspicious
%   and you need to know which FT file + opt snapshot produced it.
%
% USAGE:
%   prov = buildLFPProvenance(sourceFT, opt);
%   prov = buildLFPProvenance(sourceFT, opt, struct('area','NCL','align','stim2'));
%
% INPUTS:
%   sourceFT - char, absolute path to the FT file this output derived
%              from ('' when unknown; the sourceMtime field is NaN then).
%   opt      - resolved options struct (snapshotted verbatim).
%   extra    - (optional) struct of extra fields merged in verbatim.
%              Use for area, alignment, band, contrast spec, cluster id,
%              or anything else useful for diagnosis.
%
% OUTPUT (struct):
%   .source_FTfile   input sourceFT
%   .source_mtime    datenum of the FT file at write time (NaN if missing)
%   .opt_snapshot    the opt struct passed in (deep-copied at save time)
%   .matlabVer       'R20xxx' version string
%   .ftVer           FieldTrip version marker if discoverable, else ''
%   .savedAt         yyyy-mm-dd HH:MM:SS timestamp
%   .host            char(getenv('COMPUTERNAME')) on Windows / hostname
%                    equivalent
%   .callingFunction top of the MATLAB call stack (immediate caller name)
%   + every field in `extra` merged verbatim.
%
% Last modified 26.06.2026 (Jesus) - new helper (LFP Pass 2).

    if nargin < 3, extra = struct(); end
    prov = struct();
    prov.source_FTfile = char(sourceFT);
    prov.source_mtime  = NaN;
    if ~isempty(sourceFT) && isfile(sourceFT)
        d = dir(sourceFT);
        if ~isempty(d), prov.source_mtime = d.datenum; end
    end
    prov.opt_snapshot   = opt;
    prov.matlabVer      = version('-release');
    prov.ftVer          = localFieldTripVersion();
    prov.savedAt        = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    prov.host           = localHostname();
    prov.callingFunction = localCaller();

    % Merge extras verbatim.
    fns = fieldnames(extra);
    for k = 1:numel(fns)
        prov.(fns{k}) = extra.(fns{k});
    end
end

function v = localFieldTripVersion()
    v = '';
    try
        if exist('ft_version', 'file') == 2
            v = ft_version;
        end
    catch
        % ft_version prints to stdout instead of returning on some versions;
        % swallow.
    end
    if isempty(v), v = 'unknown'; end
end

function h = localHostname()
    h = getenv('COMPUTERNAME');            % Windows
    if isempty(h), h = getenv('HOSTNAME'); end   % *nix
    if isempty(h), h = 'unknown'; end
end

function name = localCaller()
    st = dbstack('-completenames');
    if numel(st) >= 3
        % st(1) = this helper, st(2) = buildLFPProvenance caller, st(3)+ = deeper
        name = st(3).name;
    elseif numel(st) >= 2
        name = st(2).name;
    else
        name = 'command line';
    end
end
