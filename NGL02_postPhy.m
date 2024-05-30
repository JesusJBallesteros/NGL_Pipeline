%% NGL02_postPhy (in progress)
% To run after manual curation of desired sessions is completed. Will read
% the resulting KS results after manual curation.
%
% Jesus 28.05.2024

%% Additional/default options
% opt.plotgeneral    = false; % plot in NGL03 script. General data about #clusters, units per animal/session, etc.
%                             % Useful only if extracting results from several animals and sessions
%                              
% opt.plotrasters    = false;  % plot in NGL03 script. Rasters of activity aligned to a selected t0
% opt.plotfrs        = false; % plot in NGL03 script. Firing rates calculated from the above.
% 
% % % To obtain waveforms
% opt.getwF          = false;
%     gwfparams.dataType = 'int16';  % Data type of .dat file
%     gwfparams.nCh = 32;            % Number of channels that were streamed in .dat file
%     gwfparams.wfWin = [-20 41];    % Number of samples around spiketime to include in waveform
%     gwfparams.nWf = 1;             % Proportion of total waveforms per unit to extract

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
else
    disp('Using INPUTS from NGL01_MAIN.')
end

% Set default inputs and dependencies. In case NGL01 did not before.
input = set_default(input, opt);

%% 01. Find and list requested sessions and subjects.
input.sessions = findSessions(input);

for x = 1:input.nsubjects % Subjects.
    for y = 1:input.sessions(x).nsessions % Sessions.
            input.run = [x y]; % Current run, to pass to functions.
            
            %% 02. Prepare to proceed with a single session.
            [input.sessions(input.run(1)).info, opt] = prepforsession(input, opt);           
            mkdir(opt.spikeSorted)
            mkdir(opt.trialSorted)

            %% 03. Extract preprocessed data (spikes and events)
            spike = loadSpikes(opt);
            
            %% 04. Use Events to parse spike data into trials
            % We probably can access the LFP FT file
%             FTfile = dir([opt.FolderProcDataMat, '\*tparsed_FT.mat']); % UPDATE 
%             load(FTfile.name); clear FTfile
%     
%             % If the data has been organised to reflect the temporal structure of the
%             % experiment (i.e. the trials), the SPIKE structure should contain a
%             % cell-array with the spike times relative to an experimental trigger. The
%             % FT_SPIKE_MAKETRIALS function can be used to reorganise the SPIKE
%             % structure such that the spike times are expressed relative to a trigger
%             % instead of relative to the acquisition devices internal timestamp clock.
%             % The time field then contains only those spikes that occurred within one of
%             % the trials . The spike times are now expressed on seconds relative to the
%             % trigger.
%             cfg = [];
%                 cfg.trl = FT_data.cfg.trl; % Trial times come in samples
%                 cfg.timestampspersecond = 3200; % timestamps per second to synchr spike timestamps to trial samples.
%                 % For some reason this has to be divided? (vs what I expected, 32000, as sampled from Deuteron)
%     
%             % Transform manual 'spike' into a FT ready structure. 
%             spike_trials = ft_spike_maketrials(cfg, spike);
%     
%             % Test plot of absolute spike timestamps (x-axis, in min) vs spike
%             % times relative to the reference event for trials (y-axis, in sec).
%             plot(spike.timestamp{1}/60000, spike.time{1}, '.');
%             axis([-.2 31 -.5 6])
%             xlabel('Absolute time (min)'), ylabel('Trial time (s)');
%     
%             % Combine FT_data (LFP) and spike (timestamps)
%             cfg = [];
%             % needs    
%                 FT_data.hdr.FirstTimeStamp     = FT_data.cfg.trl(1,1);
%                 FT_data.hdr.TimeStampPerSample = 1; % CHECK how this works
%     
%             Comb_data = ft_appendspike([], FT_data, spike);
            spike_trials = trialParse_spikes(opt);

%%         % FT RASTER
%         cfg              = [];
%             cfg.topplotfunc  = 'line'; % plot as a line
%             cfg.spikechannel = spike.label([1]);
%             cfg.latency      = [0 4.5];
%             cfg.errorbars    = 'std'; % plot with the standard deviation
%             cfg.interactive  = 'no'; % toggle off interactive mode
%         figure, ft_spike_plot_raster(cfg, spike, psth)

%         % FT PSTH
%         cfg             = [];
%             cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
%             cfg.outputunit  = 'rate'; % give as an output the firing rate
%             cfg.latency     = [0 6]; % between -1 and 3 sec.
%             cfg.vartriallen = 'yes'; % variable trial lengths are accepted
%             cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
%         psth = ft_spike_psth(cfg, spike);
%             

    end
end