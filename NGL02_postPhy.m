%% NGLXX_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.

%% Input storage drive, project name and toolbox folder:
input.datadrive     = 'F:\';
input.studyName     = 'Pilot_SocialLearning'; 
input.toolbox       = 'C:\Code\ephys-data-pipeline'; % Default: 'C:\Code\ephys-data-pipeline'

% To run the script on all subjects and sessions included in your project,
% just leave both as 'all'. For a session-to-session process, explicit the
% subject and session/s to process. 
input.subjects       = {'485'}; % 'all';  % char array 'all', or a cell with a single subject denomination e.g. {'DOE'} or {'042'}
input.dates          = {'20231110'}; % char array 'all', or cell array of dates for a single subject e.g. {'YYYYMMDD' ...}

%% General Options. What you want to obtain:
% Those used for all sessions. Specific options can be set below or defaulted in the functions.
opt = struct();
    opt.plotgeneral    = false; % plot in NGL03 script. General data about #clusters, units per animal/session, etc.
                                % Useful only if extracting results from several animals and sessions
                                 
    opt.plotrasters    = true;  % plot in NGL03 script. Rasters of activity aligned to a selected t0
    opt.plotfrs        = false; % plot in NGL03 script. Firing rates calculated from the above.

    % To obtain waveforms
    opt.getwF          = false;
        gwfparams.dataType = 'int16';  % Data type of .dat file
        gwfparams.nCh = 32;            % Number of channels that were streamed in .dat file
        gwfparams.wfWin = [-20 41];    % Number of samples around spiketime to include in waveform
        gwfparams.nWf = 1;             % Proportion of total waveforms per unit to extract
    
%% 00. Check inputs, set defaults and dependencies.
input = set_default(input);
opt.postPhy = true; % Always true, to differentiate from NGL01

for s = 1:input.nsubjects
    % Read requested sessions from specified animal folder.
    sessions = findSessions(input, opt);

    %% 01. Loop subjects and sessions to process.
    for ss = 1:sessions(s).nsessions
        extract = [];

        % Navigate to session's raw data folder.
        cd(fullfile(sessions(s).folder,sessions(s).list{ss}));
                
        % Progress report.
        fprintf(sprintf('\n --> Subject %s, session %d out of %d: %s \n', ...
                     input.subjects(s).name, ss, sessions(s).nsessions, sessions(s).list{ss}));

        % Determine where processed data will be saved, done for every session.
        opt.PathRaw           = pwd;
        opt.FolderProcDataMat = fullfile(input.sorted, input.subjects(s).name, sessions(s).list{ss});
        opt.SavFileName       = sessions(s).list{ss}; 
        
        % Report and create folder.
        disp(strcat('Processed data will be saved to: >', opt.FolderProcDataMat));
        if exist(opt.FolderProcDataMat,"dir") == 0
            mkdir(opt.FolderProcDataMat);
        end

        %% 02. Proceed with reading data from preprocessed files
        % Load information as stored post-Phy curation.
        try     [extract] = loadKSdir(opt.PathRaw); 
        catch,  [extract.spikeTemplates] = [];
        end

        %% 03. Add extracted spikes to a FT compatible structure
        % The output spike structure could contain the following (for now)
        %   spike.label     = 1xNchans cell-array, with unit labels
        %   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
        %   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nleads x Nsamples X Nspikes)
        % Described in more detail in FT_DATATYPE_SPIKE
        if opt.getwF
            % To get the waveforms of the spikes (not kilosort's templates), use 'getWaveForms'.
            % It gets a fraction of the total waveforms, from permutated data (representative 
            % from all long session, vs. simply all the spikes from the first third of the time.)
            % If that requested fraction is under 1000 waveforms, and the
            % cluster has at least 1000 waveforms, that will be the number
            % of waveforms extracted.
            gwfparams.dataDir = opt.PathRaw;         % .dat file containing the raw
            gwfparams.apD = dir(fullfile(opt.PathRaw, '*.bin')); 
            gwfparams.fileName = gwfparams.apD(1).name;         
        end

        % Prepare spike structure
        clusters        = sort(unique(extract.cids)); % get unique clusters and sort
        nclust          = numel(clusters); % number of clusters
        spike.label     = cell(1,nclust); % for clusters IDs
        spike.timestamp = cell(1,nclust); % spikes timestamps 

        disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.');

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

        %% 04. Use Events to parse spike data into trials
        % We probably can access the LFP FT file
        FTfile = dir(fullfile(opt.PathRaw, '*tparsed_FT.mat')); 
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


        %% TEST plots
        cfg             = [];
            cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
            cfg.outputunit  = 'rate'; % give as an output the firing rate
            cfg.latency     = [0 6]; % between -1 and 3 sec.
            cfg.vartriallen = 'yes'; % variable trial lengths are accepted
            cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
        
        psth = ft_spike_psth(cfg, spike);
            
        cfg              = [];
            cfg.topplotfunc  = 'line'; % plot as a line
            cfg.spikechannel = spike.label([1]);
            cfg.latency      = [0 4.5];
            cfg.errorbars    = 'std'; % plot with the standard deviation
            cfg.interactive  = 'no'; % toggle off interactive mode

        figure, ft_spike_plot_raster(cfg, spike, psth)

%         % TODO: trials needs to be readed from an event mat file created
%         % after a proper experiment. Not valid for testing recordings. 
%         
%         % If 'trials' does not exists, defaults here to 1 (single, long trial)
%           % smth like... if exists trials then use it, otherwise, trials = 1 
%           trials = 1;
% 
%         % TODO: create conditions
%         % This variable will be saved under '...\data\spikeSorted\...' for further access
% 
%         %% 03. Here we can create the 'neurons' cell variable, according to the IKN standard.
%         % This variable will be saved under '...\data\spikeSorted\...' for further access
%         getneurons(results{ss,s});

    end
end

%% Plotting script
% Do not clear the workspace after running NGL02
% use NGL03_plot
