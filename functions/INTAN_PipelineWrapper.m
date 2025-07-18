function input = INTAN_PipelineWrapper(input, varargin)
%
%
% Version 07.06.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults 
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FieldTrip'),       opt.FieldTrip           = true;         end
if ~isfield(opt,'doNWB'),           opt.doNWB               = false;        end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = true;         end
if ~isfield(opt,'lowpass'),         opt.lowpass             = 150;          end
if ~isfield(opt,'noise'),           opt.noise               = [];           end

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
    [~, trialdef, ~] = EventProcess(input, opt);
end

%% 03. Create NWB file
if opt.doNWB % We want a .NWB file.
% % DEPRECATE 
% % Run wrapper for the INTAN to NWB functionality:               
% % This NEEDS A PYTHON installation and the tooldbox inside!
% % Detailed explanation:
% % WHAT IT IS: function to convert data from INTAN to .NWB format.
% % WHAT IT DOES: Checks for Python engine in computer. Adds the necessary
% %  dependences. Locates input session, copies ALL files to the IntanToNWB
% %  folder and merges them into a new 'info.nwb' file. This file 
% %  is renamed to 'session_name.nwb'. Moves this new file back to 
% %  the original session folder. Removes the copied data from the 
% %  IntanToNWB folder.
% %
% % Requires Python installed in the machine. 
% %  To date, MATLAB 2021b accepts up to Python 3.9. Install the
% %  64 bits version:
% % (https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
% %  To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
% intan2NWB_wrapper(input, opt); TO DEPRECATE?

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

end