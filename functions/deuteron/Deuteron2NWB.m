function Deuteron2NWB(input, opt)
% Deuteron2NWB  Convert a Deuteron session to NWB format using matNWB.
%
% PURPOSE:
%   Loads the raw data matrix written by Deuteron2Kilosort (_raw.mat),
%   merges project-level metadata from analysisCode/nwb_metadata.yaml,
%   parses the session start time from the raw-data folder name (DDMMYYYY),
%   builds an NwbFile with an ElectricalSeries, and writes a .nwb file to
%   the preprocessing folder.  The temporary _raw.mat is deleted after a
%   successful write.
%
%   Design mirrors intan2NWB_neuroconv.m as closely as possible, but uses
%   the bundled matNWB MATLAB toolbox instead of a Python NeuroConv call,
%   because Deuteron sessions do not have a single-file header that
%   NeuroConv's IntanRecordingInterface can consume.
%
% USAGE:
%   Deuteron2NWB(input, opt)
%   Called from Deuteron_PipelineWrapper when opt.doNWB = true,
%   immediately after Deuteron2Kilosort.
%
% INPUTS:
%   input  - struct from set_default; relevant fields:
%              .analysisCode   path to project analysisCode\ folder
%   opt    - options struct; relevant fields:
%              .FolderProcDataMat  preprocessing output folder
%              .SavFileName        session name (= raw-data folder name)
%              .sampleRate         recording sample rate (Hz)
%              .numChannels        electrode count
%
% OUTPUT:
%   <SavFileName>.nwb written to opt.FolderProcDataMat.
%   Requires matNWB toolbox (bundled under toolboxes/matnwb/).
%   YAML metadata is read by readyaml() (functions/utils/readyaml.m),
%   a pure-MATLAB parser with no version dependency.
%
% CALLS:
%   NwbFile, types.core.*, types.hdmf_common.*, util.table2nwb, nwbExport
%
% SEE ALSO:
%   intan2NWB_neuroconv.m, nwb_metadata_template.yaml
%
% Last modified 15.05.2026 (Jesus)

%% Skip guard — same 1 MB threshold as INTAN path
nwbPath = fullfile(opt.FolderProcDataMat, [char(opt.SavFileName), '.nwb']);
if isfile(nwbPath)
    d = dir(nwbPath);
    if d.bytes > 1e6
        disp('- NWB file already exists for this session. Skipping.');
        return
    end
end

%% Load raw data matrix produced by Deuteron2Kilosort
rawMatPath = fullfile(opt.FolderProcDataMat, [char(opt.SavFileName), '_raw.mat']);
if ~isfile(rawMatPath)
    error('NGL:missingRawMat', ...
        ['Raw data matrix not found: %s\n', ...
         'Ensure opt.doNWB=true when running Deuteron2Kilosort.'], rawMatPath);
end
disp('- Loading raw data matrix...')
raw        = load(rawMatPath, 'data_mat', 'sampleRate_raw');
data_mat   = raw.data_mat;      % int16, [nCh × nSamples], µV
sampleRate = raw.sampleRate_raw;
nCh        = size(data_mat, 1);

%% Read project-level metadata from nwb_metadata.yaml
yamlPath = fullfile(input.analysisCode, 'nwb_metadata.yaml');
meta = struct();
if isfile(yamlPath)
    try
        meta = readyaml(yamlPath);   % requires external tool
        fprintf('- Metadata YAML loaded: %s\n', yamlPath);
    catch ME
        warning('NGL:yamlReadFailed', ...
            'readyaml() could not parse %s: %s. Proceeding with minimal metadata.', ...
            yamlPath, ME.message);
    end
else
    warning('NGL:noMetaYaml', ...
        'nwb_metadata.yaml not found in analysisCode. Proceeding with minimal metadata.');
end

% Safe nested field accessor (returns default if any level is absent)
getF = @(fields, def) safeField(meta, fields, def);

%% Load channel map from analysisCode (optional — enriches electrode table)
% The file is the same .mat used by Kilosort (opt.KSchanMapFile), expected
% at fullfile(input.analysisCode, opt.KSchanMapFile).
chanMapPath = '';
if isfield(opt, 'KSchanMapFile') && ~isempty(opt.KSchanMapFile)
    chanMapPath = fullfile(input.analysisCode, opt.KSchanMapFile);
end
useChanMap = false;
if ~isempty(chanMapPath) && isfile(chanMapPath)
    cm = load(chanMapPath, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords');
    useChanMap = true;
    fprintf('- Channel map loaded: %s\n', chanMapPath);
else
    if ~isempty(chanMapPath)
        warning('NGL:noChanMap', ...
            'Channel map not found at: %s. Using minimal electrode table.', chanMapPath);
    end
    % Minimal fallback so downstream code is path-uniform
    cm.chanMap   = (1:nCh)';
    cm.connected = true(nCh, 1);
    cm.xcoords   = zeros(nCh, 1);
    cm.ycoords   = ((0:nCh-1) * 50)';   % 50 µm linear placeholder
    cm.kcoords   = ones(nCh, 1);
end
shankIDs = unique(cm.kcoords(:), 'sorted');
nShanks  = numel(shankIDs);

sessionDesc = getF({'NWBFile','session_description'}, ...
                   'Deuteron extracellular electrophysiology recording.');
labName     = getF({'NWBFile','lab'},         'Neural Basis of Learning Lab');
experimenter= getF({'NWBFile','experimenter'}, {'JDOE'});
expDesc     = getF({'NWBFile','experiment_description'}, 'Free-moving animal');
species     = getF({'Subject','species'}, 'Columba livia');
sex         = getF({'Subject','sex'},     'U');
age         = getF({'Subject','age'},     'P5Y');
egLocation  = getF({'Ecephys','ElectrodeGroup','location'}, 'unknown');

%% Parse session start time from folder name (DDMMYYYY convention)
[~, sessionFolder] = fileparts(opt.PathRaw);
sessionStart = parseSessionDate(sessionFolder);

%% Build NWB file
disp('- Building NWB file...')
identifier = sprintf('%s_%s', char(opt.SavFileName), char(input.studyName));

nwb = NwbFile( ...
    'session_description',    sessionDesc, ...
    'identifier',             identifier, ...
    'session_start_time',     sessionStart, ...
    'file_create_date',       {datetime('now', 'TimeZone', 'local')}, ...
    'general_lab',                    labName, ...
    'general_experiment_description', expDesc);

if ~isempty(experimenter)
    if ischar(experimenter), experimenter = {experimenter}; end
    nwb.general_experimenter = experimenter;
end

%% Subject (optional — only if species is provided in YAML)
if ~isempty(species)
    nwb.general_subject = types.core.Subject( ...
        'species',     species, ...
        'sex',         sex, ...
        'age',         age, ...
        'description', 'See nwb_metadata.yaml for project-level subject details.');
end

%% Device
nwb.general_devices.set('Deuteron', types.core.Device( ...
    'description', ...
    'Deuteron Technologies Neurolog miniature wireless logger. Recorded at 32 kHz.'));

%% ElectrodeGroups — one per shank (kcoords group) in the channel map
% Single-shank probes produce one group ('shank0'); multi-shank probes
% produce shank0, shank1, … matching unique kcoords values.
for s = 1:nShanks
    egName = sprintf('shank%d', s - 1);
    nwb.general_extracellular_ephys.set(egName, types.core.ElectrodeGroup( ...
        'description', sprintf('Channels on shank %d.', shankIDs(s)), ...
        'location',    egLocation, ...
        'device',      types.untyped.SoftLink('/general/devices/Deuteron')));
end

%% Electrode table — one row per channel with geometry from chanMap
% Columns id, location, group, group_name are required by NWB.
% rel_x / rel_y are probe-relative coordinates in µm (from xcoords/ycoords).
% shank_id mirrors kcoords; connected flags non-noisy channels.
% egOVs must be an object array (not a cell array) — matNWB's io.mapData2H5
% rejects cell arrays containing non-char content.  Indexed assignment into
% an uninitialized variable builds the ObjectView array element by element.
egNames = cell(nCh, 1);
for c = 1:nCh
    sIdx          = find(shankIDs == cm.kcoords(c), 1);
    egN           = sprintf('shank%d', sIdx - 1);
    egOVs(c, 1)   = types.untyped.ObjectView( ...           
                        sprintf('/general/extracellular_ephys/%s', egN));
    egNames{c}    = egN;
end

tbl = table( ...
    double(cm.chanMap(:) - 1), ...       % 0-indexed physical channel id
    repmat({egLocation}, nCh, 1), ...    % brain region label
    egOVs, ...                           % ObjectView to ElectrodeGroup
    egNames, ...                         % ElectrodeGroup name string
    double(cm.xcoords(:)), ...           % lateral position, µm (probe-relative)
    double(cm.ycoords(:)), ...           % depth position,  µm (probe-relative)
    double(cm.kcoords(:)), ...           % shank index
    logical(cm.connected(:)), ...        % false = noisy / dead channel
    'VariableNames', ...
        {'id','location','group','group_name','rel_x','rel_y','shank_id','connected'});

nwb.general_extracellular_ephys_electrodes = util.table2nwb(tbl, ...
    ['Electrode table for Deuteron recording. ' ...
     'rel_x/rel_y in µm (probe-relative). ' ...
     'connected=false marks noisy or dead channels.']);

%% ElectrodeTableRegion (all channels)
etr = types.hdmf_common.DynamicTableRegion( ...
    'table',       types.untyped.ObjectView('/general/extracellular_ephys/electrodes'), ...
    'description', 'All electrode channels.', ...
    'data',        (0:nCh-1)');

%% ElectricalSeries
% data_mat is int16 [nCh × nSamples] in µV.
% NWB convention: first dimension is time, so we transpose.
% data_conversion = 1e-6 scales µV → V (the NWB standard unit).
es = types.core.ElectricalSeries( ...
    'description',        'Raw wideband recording (ADC-µV only; no filtering applied).', ...
    'data',               data_mat', ...        % [nSamples × nCh]
    'data_conversion',    single(1e-6), ...     % µV -> V
    'data_unit',          'volts', ...
    'filtering',          'none', ...
    'starting_time',      0.0, ...
    'starting_time_rate', single(sampleRate), ...
    'electrodes',         etr);
nwb.acquisition.set('ElectricalSeries', es);

%% Write
disp('- Writing NWB file...')
nwbExport(nwb, nwbPath);
fprintf('- NWB file written: %s\n', nwbPath);

%% Clean up temporary raw mat
delete(rawMatPath);
disp('- Temporary _raw.mat removed.')

end % main function

% HELPER FUNCTIONS
function val = safeField(s, fields, default)
% SAFEFIELD  Return nested struct field; return default if any level absent.
    val = default;
    try
        tmp = s;
        for k = 1:numel(fields)
            tmp = tmp.(fields{k});
        end
        if ~isempty(tmp)
            val = tmp;
        end
    catch
    end
end

function t = parseSessionDate(folderName)
% PARSESESSIONDATE  Parse YYYYMMDD from raw-data folder name.
%   Different from INTAN, but INTAN is the one that will need to be check
%   Returns a datetime with Europe/Berlin timezone.
%   Falls back to current time with a warning if parsing fails.
    try
        day   = str2double(folderName(7:8));
        month = str2double(folderName(5:6));
        year  = str2double(folderName(1:4));
        t = datetime(year, month, day, 0, 0, 0, 'TimeZone', 'Europe/Berlin');
        fprintf('- Session start time parsed from folder ''%s'': %s\n', ...
            folderName, char(t));
    catch
        warning('NGL:sessionDateParseFailed', ...
            ['Could not parse session date from folder name ''%s'' ', ...
             '(expected DDMMYYYY). Using current time as fallback.'], folderName);
        t = datetime('now', 'TimeZone', 'local');
    end
end
