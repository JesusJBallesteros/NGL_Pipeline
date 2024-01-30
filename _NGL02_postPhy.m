%% NGLXX_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 03.01.20024

opt.plotgeneral    = false; % plot in NGL03 script. General data about #clusters, units per animal/session, etc.
                            % Useful only if extracting results from several animals and sessions
                             
opt.plotrasters    = true;  % plot in NGL03 script. Rasters of activity aligned to a selected t0
opt.plotfrs        = false; % plot in NGL03 script. Firing rates calculated from the above.

% % To obtain waveforms
% opt.getwF          = false;
%     gwfparams.dataType = 'int16';  % Data type of .dat file
%     gwfparams.nCh = 32;            % Number of channels that were streamed in .dat file
%     gwfparams.wfWin = [-20 41];    % Number of samples around spiketime to include in waveform
%     gwfparams.nWf = 1;             % Proportion of total waveforms per unit to extract
    
%% 00. Check inputs, set defaults and dependencies.
% input = set_default(input);
% opt.postPhy = true; % Always true, to differentiate from NGL01

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

%         %% 02. Proceed with reading data from preprocessed files
%         % Load information as stored post-Phy curation.
%         try     [extract] = loadKSdir(opt.PathRaw); 
%         catch,  [extract.spikeTemplates] = [];
%         end
% 
%         %% 03. Add extracted spikes to a FT compatible structure
%         % The output spike structure could contain the following (for now)
%         %   spike.label     = 1xNchans cell-array, with unit labels
%         %   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%         %   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nleads x Nsamples X Nspikes)
%         % Described in more detail in FT_DATATYPE_SPIKE
%         if opt.getwF
%             % To get the waveforms of the spikes (not kilosort's templates), use 'getWaveForms'.
%             % It gets a fraction of the total waveforms, from permutated data (representative 
%             % from all long session, vs. simply all the spikes from the first third of the time.)
%             % If that requested fraction is under 1000 waveforms, and the
%             % cluster has at least 1000 waveforms, that will be the number
%             % of waveforms extracted.
%             gwfparams.dataDir = opt.PathRaw;         % .dat file containing the raw
%             gwfparams.apD = dir(fullfile(opt.PathRaw, '*.bin')); 
%             gwfparams.fileName = gwfparams.apD(1).name;         
%         end
% 
%         % Prepare spike structure
%         clusters        = sort(unique(extract.cids)); % get unique clusters and sort
%         nclust          = numel(clusters); % number of clusters
%         spike.label     = cell(1,nclust); % for clusters IDs
%         spike.timestamp = cell(1,nclust); % spikes timestamps 
% 
%         disp('Extracting curated clusters from Phy files. If many waveforms are requested, it may take a while.');
% 
%         % Loop clusters
%         for cl = 1:nclust
%             % Proceed to extract details
%             spike.label{cl}           = num2str(clusters(cl));
%             spike.timestamp{cl}       = extract.st(extract.clu==clusters(cl))*1000; % from ms to sec
% 
%             if opt.getwF
%                 % a few more params for 'getWaveForms' dep on cluster
%                 gwfparams.spikeTimes = ceil(extract.st(extract.clu==clusters(cl))*32000); % Vector of cluster spike times (in samples) same length as .spikeClusters
%                 gwfparams.spikeClusters = extract.clu(extract.clu==clusters(cl));
%                 
%                 % Get waveform
%                 wF = getWaveForms(gwfparams);
%     
%                 % Refine
%                 wF.waveForms = squeeze(wF.waveForms);
%                 wF.waveFormsMean = squeeze(wF.waveFormsMean);
% 
%                 % Find averaged max amplitude channel
%                 [wF.maxamplch, ~, ~] = find(wF.waveFormsMean==min(min(wF.waveFormsMean)));
%                 spike.waveform{cl}   = permute(wF.waveForms(:,wF.maxamplch,:),[2,3,1]);
%                 
%                 gwfparams.spikeClusters = [];
%                 gwfparams.spikeTimes = [];
%                 wF = [];
%             end
%         end
% 
%         clear clusters nclust cl
% 

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