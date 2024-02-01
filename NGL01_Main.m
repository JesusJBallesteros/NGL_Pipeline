%% Pipeline process INTAN and Deuteron continous data.
% Will read and process INTAN, DEUTERON (or ALLEGO) data, from selected sessions for a given animal.
% The main pipeline will be: INTAN/DEUTERON raw formats to be located, then
% converted to .bin files (spike sorting), and Fieldtrip .mat structures
% (for LFP). Once sorted, spike data will be attached to the FieldTrip
% structure. For Arena experiments, motion sensor data will be extracted and
% interpreted. Data will be trial-parsed using EventCodes.
%
% The hard disk data structure SHOULD fit the IKN standard published at:
% gitlab.ruhr-uni-bochum.de/ikn/howto/-/wikis/Neurophysiology/hard-disk-data-structure
%
% DEPENDENCIES:
% Requires that all pipeline dependencies are properly located. 
% I suggest to include the 'mainfolder' in Matlab's permanent path system.
% The function 'set_default' will take care of the rest of folders on each run.
%
% INPUTS:
%       input.datadrive, char array with the drive where data is located. As 'D:\'
%       input.studyName, char array with the project name, matching the
%                          folder name where all data will be stored. As 'studyName'
%       input.subjects,  char array with either 'all' OR a single subject name e.g. 'DOE'
%       input.dates,     char array with either 'all' OR a cell array of dates 
%                          for a SINGLE subject e.g. {'YYYYMMDD' 'yyyymmdd' ...)
%
% OPTIONS: is a struct with many possible fields. All should have a
% corresponding default inside whatever function is being called. Main ones
% are:     
%     opt.bin,              Creation of .bin file, input to Kilosort 2/4.
%     opt.FTfile,           Creation of .mat file with FieldTrip format.
%     opt.RetrieveEvents,   Retrieve event log from Deuteron system.
%     opt.GetMotionSensors, Retrieve data from motion sensors in Deuteron.
%     opt.kilosort,         Asks to proceed with KS processing and waits to retrieve its results.
%     opt.set_filter,       If Deuteron data was adquired with a wideband.
%     opt.lowpass,          Lowpass band to extract LFP from wideband.
%     opt.highpass,         Highpass band to extract spike activity.
%
% OUTPUTS:
% For one single session or for a batch of sessions, from one single animal:
%       Fieldtrip (.mat), binary (.bin), HDF5 (.h5) and/or .nwb files from
%           1. Deuteron .DT2 or .DF1 data.
%           2. INTAN file-per-type and file-per-channel format data.
%           3. (ALLEGO data?)
%       EventRecord.mat file, from Deuteron session.
%       MotionData.mat file, From Deuteron sensors.
%       Plots snippets of time- and frequency-domain data, from FieldTrip
%       
% Last modified 26.04.2023 (Jesus)

% TODO LIST
% Order channels as incremental ordinals. (in 'Deuteron_ExtractEvents' ~136)
% If Deuteron2Kilosort(opt) filter for DF1 format works, set filter out of format cases (generalize)
% Continue with 'Deuteron_GetDigInEvents' when we get a recording with EVENTS
% Check for Deuteron_GetDigInEvents(EventRecord) status.
% Check for FT trial-parsing using EventRecord with MAT2FieldTrip(data, opt, varargin)
%    Create a 'trial-parsed' stream in 'mat2FieldTrip' VS. add post-hoc parsing
% Extract nChannels from EventsRecord. Find first 'File started' then use
%    'strsplit(EventRecord(50).Details,{';','='})' and find the 6th cell
% There seems to be an ERROR on 2nd and following runs of the NWB functionalities.
%    Figure out what's going on with the NWB/H5 DLLs that block either when the other has been performed...

% Version 01.02.2024 (Jesus)

%% 00. Check current inputs.
% Check if input variable exist already. Parse values.
if ~exist("input","var")
    input = struct( 'datadrive' , datadrive , ...   % force char array
                    'studyName' , studyname , ...   % force char array
                    'toolbox'   , toolbox   , ...   % force char array
                    'subjects'  , [], ...           % do NOT force char array
                    'dates'     , []        );      % do NOT force char array
    input.dates     = dates;    % place as it comes
    input.subjects  = subjects; % place as it comes
end

% Set default inputs and dependencies.
input = set_default(input);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
        input.run = [x y]; % Current run, to pass to functions.
        %% 02. Prepare to proceed with a single session.
        [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);
    
        %% 03. Proceed to appropiated pipeline.
        switch input.sessions(input.run(1)).info.fileformat
            case {'DT2', 'DF1'} 
               % 03.1 Deuteron Pipeline
               opt = Deuteron_PipelineWrapper(input, opt);
    
            case {'fileperch', 'filepertype'}
               % 03.2 INTAN Pipeline
               INTAN_PipelineWrapper(input, opt);
    
            otherwise
               warning('Something went wrong during format verification. Skipping Session');
               continue
        end 
        %% 04. Kilosort
        if opt.kilosort
            % Kilosort will run without GUI.
            master_kilosort(input, opt)
        end
    
        %% 05. Bombcell
        if opt.bombcell
            % Kilosort will run without GUI.
            Bombcell_Main(opt) 
        end

        %% 06. Open Phy to manual curation or just inspection
        if opt.phy
            % Will change to current session directory and open phy.
            % ! Keeps MATLAB busy until interface is closed.
            cd(opt.FolderProcDataMat)
            system('phy template-gui params.py');
        end

        %% 07. Read sorted and curated clusters to MATLAB. (TO BE WRAPPED)
        if opt.cooking
            % To obtain waveforms
            opt.getwF          = false;
                gwfparams.dataType = 'int16';  % Data type of .dat file
                gwfparams.nCh = 32;            % Number of channels that were streamed in .dat file
                gwfparams.wfWin = [-20 41];    % Number of samples around spiketime to include in waveform
                gwfparams.nWf = 1;             % Proportion of total waveforms per unit to extract
    
            % Proceed with reading data from preprocessed files
            % Load information as stored post-Phy curation.
            disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.')
            [extract] = loadKSdir(opt.PathRaw); 
    
            % Prepare spike structure
            clusters        = sort(unique(extract.cids)); % get unique clusters and sort
            nclust          = numel(clusters); % number of clusters
            spike.label     = cell(1,nclust); % for clusters IDs
            spike.timestamp = cell(1,nclust); % spikes timestamps 
    
            % Loop clusters
            for cl = 1:nclust
                % Proceed to extract details
                spike.label{cl}           = num2str(clusters(cl));
                spike.timestamp{cl}       = extract.st(extract.clu==clusters(cl))*1000; % from ms to sec
    
                if opt.getwF
                    % a few more params for 'getWaveForms' dep on cluster
                    gwfparams.spikeTimes = ceil(extract.st(extract.clu==clusters(cl))*32000); % Vector of cluster spike times (in samples) same length as .spikeClusters
                    gwfparams.spikeClusters = extract.clu(extract.clu==clusters(cl));
                    
                    % Get waveform
                    wF = getWaveForms(gwfparams);
        
                    % Refine
                    wF.waveForms = squeeze(wF.waveForms);
                    wF.waveFormsMean = squeeze(wF.waveFormsMean);
    
                    % Find averaged max amplitude channel
                    [wF.maxamplch, ~, ~] = find(wF.waveFormsMean==min(min(wF.waveFormsMean)));
                    spike.waveform{cl}   = permute(wF.waveForms(:,wF.maxamplch,:),[2,3,1]);
                    
                    gwfparams.spikeClusters = [];
                    gwfparams.spikeTimes = [];
                    wF = [];
                end
            end
            clear clusters nclust cl
    
            %% 08. Use Events to parse spike data into trials (TO BE WRAPPED)
            % We probably can access the LFP FT file
            FTfile = dir(fullfile(opt.PathRaw, '*tparsed_FT.mat')); % UPDATE 
            load(FTfile.name); clear FTfile
    
            % If the data has been organised to reflect the temporal structure of the
            % experiment (i.e. the trials), the SPIKE structure should contain a
            % cell-array with the spike times relative to an experimental trigger. The
            % FT_SPIKE_MAKETRIALS function can be used to reorganise the SPIKE
            % structure such that the spike times are expressed relative to a trigger
            % instead of relative to the acquisition devices internal timestamp clock.
            % The time field then contains only those spikes that occurred within one of
            % the trials . The spike times are now expressed on seconds relative to the
            % trigger.
            cfg = [];
                cfg.trl = FT_data.cfg.trl; % Trial times come in samples
                cfg.timestampspersecond = 3200; % timestamps per second to synchr spike timestamps to trial samples.
                % For some reason this has to be downsampled? (vs what I expected, 32000, as sampled from Deuteron)
    
            % Transform manual 'spike' into a FT ready structure. 
            spike = ft_spike_maketrials(cfg, spike);
    
            % Test plot of absolute spike timestamps (x-axis, in min) vs spike
            % times relative to the reference event for trials (y-axis, in sec).
            plot(spike.timestamp{1}/60000, spike.time{1}, '.');
            axis([-.2 31 -.5 6])
            xlabel('Absolute time (min)'), ylabel('Trial time (s)');
    
            % Combine FT_data (LFP) and spike (timestamps)
            cfg = [];
            % needs    
                FT_data.hdr.FirstTimeStamp     = FT_data.cfg.trl(1,1);
                FT_data.hdr.TimeStampPerSample = 1; % CHECK how this works
    
            Comb_data = ft_appendspike([], FT_data, spike);
        end

        %% Clean up to move on to next session
        clear FT_data INTANdata txt

    end % sessions loop
end % subjects loop