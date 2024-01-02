function sessions = INTAN_PipelineWrapper(input, varargin)
%
%
% Version 02.01.2024 (Jesus)

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults 
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FTfile'),          opt.FTfile              = true;         end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end

%% 01. Find out INTAN settings and header file. Extract info.
%  Uses a modified Intan function, to make the basic information
%  available at 'info{ss}' and a more detailed info at
%  the '.INTAN_hdr' sub-structure.
input.sessions(input.run(1)) = findSetting(input.sessions(input.run(1)));

%% 02. Create NWB file
if input.useNWB % We want a .NWB file.

  % Run wrapper for the INTAN to NWB functionality:               
    % This NEEDS A PYTHON installation and the tooldbox inside!
    % Detailed explanation:
    % WHAT IT IS: function to convert data from INTAN to .NWB format.
    % WHAT IT DOES: Checks for Python engine in computer. Adds the necessary
    %  dependences. Locates input session, copies ALL files to the IntanToNWB
    %  folder and merges them into a new 'info.nwb' file. This file 
    %  is renamed to 'session_name.nwb'. Moves this new file back to 
    %  the original session folder. Removes the copied data from the 
    %  IntanToNWB folder.
    %
    % Requires Python installed in the machine. 
    %  To date, MATLAB 2021b accepts up to Python 3.9. Install the
    %  64 bits version:
    % (https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
    %  To check access to Python Modules from MATLAB, look that 'pe' is correctly populated when running the script.
  intan2NWB_wrapper(input, opt);
end 

%% 03. Run wrapper for the INTAN to Kilosort. Creates .bin and .h5 files
if input.ExtractData 
    if opt.bin && ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']))
        % Based on Sara, Aylin and Lukas' scripts.
        Intan2Kilosort_wrapper(input.sessions(input.run(1)), opt);
    end
end

%% 04. Run wrapper for the INTAN to FIELDTRIP.
if opt.FTfile && ~isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName '_continous_FT.mat']))
    % Includes a mix of INTAN funtions. CREATES and GIVES proper
    % FieldTrip format without trial-parsing. 
    intan2FieldTrip(input.sessions(input.run(1)), opt)

    % 04.1 Plotting. Uses Chronux Multitaper approach to generate fast
    % single-tappered Spectrograms on a subset of channels for a small chunck
    % of time. Just to have a preview of how the signal looks like in
    % the LFP range.
%     if isfield(input, 'test_ch') && ~isempty(input.test_ch)
%         plot_testsignal(FT_data, input.test_ch, opt)
%     end
end

end