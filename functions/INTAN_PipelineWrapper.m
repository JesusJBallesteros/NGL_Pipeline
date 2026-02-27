function input = INTAN_PipelineWrapper(input, varargin)
% Version 07.06.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults 
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FieldTrip'),       opt.FieldTrip           = false;        end
if ~isfield(opt,'doNWB'),           opt.doNWB               = false;        end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end
if ~isfield(opt,'noise'),           opt.noise               = [];           end
if ~isfield(opt,'lowpass'),         opt.lowpass             = 10000;        end
if ~isfield(opt,'lowpassFT'),       opt.lowpassFT           = 250;          end
if ~isfield(opt,'highpass'),        opt.highpass            = 0;            end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end
if ~isfield(opt,'CAR'),             opt.CAR                 = false;        end

%% 01. Find out INTAN settings and header file. Extract info.
%  Uses a modified Intan function, to make the basic information
%  available at 'info{ss}' and a more detailed info at
%  the '.INTAN_hdr' sub-structure.
input.sessions(input.run(1)) = findSetting(input.sessions(input.run(1)));

%% 02. Event data retrieval and trial definition. INTAN version
if isfile(fullfile(opt.trialSorted, "trialdef.mat"))
    load(fullfile(opt.trialSorted, "trialdef.mat"))
    if ~exist("trialdef","var") && exist("trialDefinition","var")
        trialdef = trialDefinition.trl; clear trialDefinition
    end
else
    % 'trialdef' outputted for later feed into fieldtrip transf.
    [~, trialdef, ~, opt] = EventProcess(input, opt);
end

%% 03. Create NWB file
if opt.doNWB % We want a .NWB file.
    % Run NeuroConv python app for the INTAN to NWB conversion:               
    % This NEEDS A PYTHON installation in the corresponding Conda Enviroment!
    % Detailed explanation:
    % WHAT IT IS: function to convert data from INTAN to .NWB format.
    % WHAT IT DOES: Checks for Python engine in computer. Checks given input 
    %  path to .rhs or .rhd file. Converts to output path .nwb file.
    %
    % Requires Python installed in the machine. Requires Neuroconv installed (in proper py enviroment).
    % To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
    intan2NWB_neuroconv(input, opt);  
end 

%% 04. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
if opt.bin && ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']))
    % Based on Sara, Aylin and Lukas' scripts.
    Intan2Kilosort_wrapper(input.sessions(input.run(1)), opt);
end

%% 05. Run functions to convert INTAN dat to FIELDTRIP structure.
if opt.FieldTrip && ~isfile(fullfile(opt.analysis,[opt.SavFileName '_FTcont.mat']))
    % Includes a mix of INTAN funtions. CREATES and GIVES proper
    % FieldTrip format without trial-parsing.
    INTANdata = intan2MAT_wrapper(input.sessions(input.run(1)), opt);
    MAT2FieldTrip(INTANdata, opt, trialdef, 1); %(data, options, trialdefinitions, do continuous)
end

%% Get Motion Data into Matlab
if opt.GetMotionSensors
    disp('Extracting Motion Sensor data ...')
    GetMotionSensors(opt, input);
end
end