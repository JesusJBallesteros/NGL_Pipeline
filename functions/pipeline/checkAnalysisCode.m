function checkAnalysisCode(input, opt)
% checkAnalysisCode  Validate that all required config files exist in analysisCode\.
%
% PURPOSE:
%   Called from set_default (Section 7) after input.analysisCode is resolved.
%   Checks for required and recommended config files based on the active
%   pipeline options (opt). Missing required files abort with an error that
%   lists every gap at once so the user can fix them all in one go. Missing
%   recommended files produce warnings only.
%
% REQUIRED FILES (always):
%   NGL_machineConfig.m      machine-specific paths (Python envs, toolbox root)
%   chanMap*.mat             at least one Kilosort channel map
%   master_kilosort4.py
%   parameters.py | parameters_<area>.py | def_parameters.py
%   bombcellConfig.m
%   master_neuroconv.py
%   nwb_metadata.yaml
%   eventDefinitions.m
%
% RECOMMENDED FILES (warning only):
%   conditions_script.m    trial condition definitions
%   postPhy_param.m        NGL02 post-Phy parameters
%
% TEMPLATE SOURCE:
%   All config files except NGL_machineConfig.m and chanMap*.mat have
%   templates in the pipeline repo under configfiles\. Copy the relevant
%   templates into your project's analysisCode\ and customise them.
%
% USAGE:
%   checkAnalysisCode(input, opt)
%   Called automatically from set_default — do not call directly.
%
% Last modified 12.05.2026 (Jesus)

ac = input.analysisCode;   % short alias

missing_required    = {};
missing_recommended = {};

%% Always required
chanMaps = dir(fullfile(ac, 'chanMap*.mat'));
if isempty(chanMaps)
    missing_required{end+1} = 'chanMap*.mat  (at least one Kilosort channel map)';
end

%% Kilosort
if opt.kilosort
    if ~isfile(fullfile(ac, 'master_kilosort4.py'))
        missing_required{end+1} = 'master_kilosort4.py  (Kilosort 4 Python wrapper)';
    end

    hasParams = isfile(fullfile(ac, 'parameters.py')) || ...
                isfile(fullfile(ac, 'def_parameters.py')) || ...
                ~isempty(dir(fullfile(ac, 'parameters_*.py')));
    if ~hasParams
        missing_required{end+1} = ...
            'parameters.py / parameters_<Area>.py / def_parameters.py  (Kilosort 4 run parameters)';
    end
end

%% Bombcell
if opt.bombcell
    if ~isfile(fullfile(ac, 'bombcellConfig.m'))
        missing_required{end+1} = 'bombcellConfig.m  (Bombcell quality metric thresholds)';
    end
end

%% NWB
if opt.doNWB
    if ~isfile(fullfile(ac, 'master_neuroconv.py'))
        missing_required{end+1} = 'master_neuroconv.py  (NeuroConv NWB conversion wrapper)';
    end
    if ~isfile(fullfile(ac, 'nwb_metadata.yaml'))
        missing_required{end+1} = 'nwb_metadata.yaml  (project-level NWB metadata)';
    end
end

%% Event retrieval
if opt.RetrieveEvents
    if ~isfile(fullfile(ac, 'eventDefinitions.m'))
        missing_required{end+1} = 'eventDefinitions.m  (event code definitions)';
    end
end

%% Recommended (always)
if ~isfile(fullfile(ac, 'conditions_script.m'))
    missing_recommended{end+1} = 'conditions_script.m  (trial condition labels)';
end

if ~isfile(fullfile(ac, 'postPhy_param.m'))
    missing_recommended{end+1} = 'postPhy_param.m  (NGL02 post-Phy parameters)';
end

%% Report

if ~isempty(missing_recommended)
    fprintf('\n[checkAnalysisCode] NOTICE: recommended files absent from %s\n', ac);
    for i = 1:numel(missing_recommended)
        fprintf('    missing: %s\n', missing_recommended{i});
    end
    fprintf('  Templates are in configfiles\\ of the pipeline repository.\n\n');
end

if ~isempty(missing_required)
    fprintf('\n[checkAnalysisCode] ERROR: required config files missing from:\n  %s\n\n', ac);
    for i = 1:numel(missing_required)
        fprintf('    MISSING: %s\n', missing_required{i});
    end
    fprintf('\n  Copy the relevant templates from configfiles\\ in the pipeline\n');
    fprintf('  repository into your analysisCode\\ folder and customise them.\n\n');
    error('NGL:missingConfigFiles', ...
        '%d required file(s) missing from analysisCode\\. See list above.', ...
        numel(missing_required));
end

fprintf('[checkAnalysisCode] analysisCode\\ OK (%s)\n', ac);

end
