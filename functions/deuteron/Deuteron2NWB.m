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
%   YAML parsing requires MATLAB R2023b+; older versions proceed with
%   minimal metadata and emit a warning.
%
% CALLS:
%   NwbFile, types.core.*, types.hdmf_common.*, util.table2nwb, nwbExport
%
% SEE ALSO:
%   intan2NWB_neuroconv.m, nwb_metadata_template.yaml
%
% Last modified 13.05.2026 (Jesus)

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
        meta = yamlread(yamlPath);   % requires R2023b+
        fprintf('- Metadata YAML loaded: %s\n', yamlPath);
    catch
        warning('NGL:yamlReadFailed', ...
            ['yamlread() failed — requires MATLAB R2023b+. ', ...
             'Proceeding with minimal NWB metadata.']);
    end
else
    warning('NGL:noMetaYaml', ...
        'nwb_metadata.yaml not found in analysisCode. Proceeding with minimal metadata.');
end

% Safe nested field accessor (returns default if any level is absent)
getF = @(fields, def) safeField(meta, fields, def);

sessionDesc = getF({'NWBFile','session_description'}, ...
                   'Deuteron extracellular electrophysiology recording.');
labName     = getF({'NWBFile','lab'},         'Neural Basis of Learning Lab');
institution = getF({'NWBFile','institution'}, 'Ruhr-Universität Bochum');
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
    'lab',                    labName, ...
    'institution',            institution, ...
    'experiment_description', expDesc);

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

%% ElectrodeGroup
nwb.general_extracellular_ephys.set('ElectrodeGroup0', types.core.ElectrodeGroup( ...
    'description', 'All recorded channels.', ...
    'location',    egLocation, ...
    'device',      types.untyped.SoftLink('/general/devices/Deuteron')));

%% Electrode table (one row per channel, minimal columns)
egOV = types.untyped.ObjectView('/general/extracellular_ephys/ElectrodeGroup0');
tbl  = table( ...
    (0:nCh-1)', ...
    repmat({egLocation},       nCh, 1), ...
    repmat({egOV},             nCh, 1), ...
    repmat({'ElectrodeGroup0'}, nCh, 1), ...
    'VariableNames', {'id', 'location', 'group', 'group_name'});
nwb.general_extracellular_ephys_electrodes = util.table2nwb(tbl, ...
    'Electrode table for Deuteron recording.');

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
% PARSESESSIONDATE  Parse DDMMYYYY from raw-data folder name.
%   Returns a datetime with Europe/Berlin timezone.
%   Falls back to current time with a warning if parsing fails.
    try
        day   = str2double(folderName(1:2));
        month = str2double(folderName(3:4));
        year  = str2double(folderName(5:8));
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
