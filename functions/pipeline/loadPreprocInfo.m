function preprocInfo = loadPreprocInfo(folderPath, mode)
% loadPreprocInfo  Read a preprocInfo snapshot written by savePreprocInfo.
%
% PURPOSE:
%   Companion reader for savePreprocInfo. Used by NGL02_postPhy (and any
%   downstream stage) to recover the options that NGL01_Main actually used
%   for a given session or study, even after MATLAB has been restarted.
%
% USAGE:
%   info = loadPreprocInfo(input.analysisCode,    'master');
%   info = loadPreprocInfo(opt.FolderProcDataMat, 'session');
%
% INPUTS:
%   folderPath - char, directory to read from.
%                 'master'  -> the project's analysisCode folder
%                 'session' -> a session's preprocessing folder (the same
%                              folder savePreprocInfo wrote to:
%                              opt.FolderProcDataMat per current convention)
%   mode       - 'master'  reads <folderPath>\preprocInfo_lastRun.mat
%                'session' reads <folderPath>\preprocInfo.mat
%
% OUTPUT:
%   preprocInfo - struct with fields written by savePreprocInfo
%                 (.opt, .Areas, .runDate, .toolboxVersion, .MATLABversion,
%                 and in 'session' mode also .subject, .session).
%
% ERRORS:
%   Throws NGL:loadPreprocInfo:notFound if the expected file is missing,
%   and NGL:loadPreprocInfo:badFile if the file exists but does not contain
%   a 'preprocInfo' variable of the expected shape. Callers that want
%   silent best-effort behaviour should wrap in try/catch.
%
% SEE ALSO:
%   savePreprocInfo, applyPreprocInfo
%
% Last modified 27.05.2026 (Jesus)

%% Validate inputs
assert(ischar(folderPath) && ~isempty(folderPath), ...
    'NGL:loadPreprocInfo', 'folderPath must be a non-empty char.');
assert(ischar(mode) && ismember(mode, {'master','session'}), ...
    'NGL:loadPreprocInfo', 'mode must be ''master'' or ''session''.');

%% Resolve target file name from mode
switch mode
    case 'master',  fname = 'preprocInfo_lastRun.mat';
    case 'session', fname = 'preprocInfo.mat';
end
target = fullfile(folderPath, fname);

%% Check existence
if ~isfile(target)
    error('NGL:loadPreprocInfo:notFound', ...
          'preprocInfo snapshot not found at: %s', target);
end

%% Load and shape-check
S = load(target, 'preprocInfo');
if ~isfield(S, 'preprocInfo') || ~isstruct(S.preprocInfo) || ~isfield(S.preprocInfo,'opt')
    error('NGL:loadPreprocInfo:badFile', ...
          ['File %s does not contain a valid preprocInfo struct ', ...
           '(expected a struct with at least an .opt field).'], target);
end

preprocInfo = S.preprocInfo;
end
