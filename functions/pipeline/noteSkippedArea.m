function txtPath = noteSkippedArea(opt, areaName, reason)
% noteSkippedArea  Record on disk that a session/area was skipped during NGL02.
%
% PURPOSE:
%   Called from NGL02_postPhy when a per-area pipeline branch can't run
%   (typically: loadSpikes errored, or returned zero clusters). Writes a
%   small human-readable text file next to the session outputs so the
%   user can see in `dir <subj>/<sess>` exactly which area was skipped
%   and why, without having to re-open MATLAB.
%
% USAGE:
%   noteSkippedArea(opt, areaName, reason)
%   txtPath = noteSkippedArea(opt, areaName, reason)
%
% INPUTS:
%   opt      - resolved opt struct. Required:
%                .analysis           per-session analysis output folder
%                .SavFileName        session name (for the message body)
%   areaName - char, e.g. 'NCL' / 'STR'. Use 'all' for single-area runs
%              when the entire session is being skipped.
%   reason   - char, free-form explanation. Multi-line OK.
%
% OUTPUT (on disk):
%   <opt.analysis>/<area>_skipped.txt
%
% Last modified 23.06.2026 (Jesus)

    assert(isfield(opt,'analysis') && ~isempty(opt.analysis), ...
        'NGL:noteSkippedArea', 'opt.analysis must be set (run after prepforsession).');
    if ~isfolder(opt.analysis), mkdir(opt.analysis); end

    safeArea = regexprep(char(string(areaName)), '[^\w\-]', '_');
    if isempty(safeArea), safeArea = 'all'; end
    txtPath = fullfile(opt.analysis, [safeArea '_skipped.txt']);

    sessName = '';
    if isfield(opt,'SavFileName'), sessName = opt.SavFileName; end

    fid = fopen(txtPath, 'w');
    if fid < 0
        warning('NGL:noteSkippedArea', 'Could not write %s', txtPath);
        return
    end
    fprintf(fid, 'NGL02 skipped area\n');
    fprintf(fid, '==================\n');
    fprintf(fid, 'session   : %s\n', sessName);
    fprintf(fid, 'area      : %s\n', areaName);
    fprintf(fid, 'timestamp : %s\n', char(datetime('now')));
    fprintf(fid, '\nReason:\n%s\n', reason);
    fclose(fid);
end
