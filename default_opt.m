function opts = default_opt()
% default_opt: Returns the default options for the NGL toolbox.
% Every recognized option must appear here. Downstream functions should
% NEVER hard-code a default value themselves, they should rely on this
% struct already being complete and validated.

    % Data format
    opts.numChannels     = 32;       % Expected channel count (override for 32-ch Deuteron)
    opts.bin             = true;     % Normally, we always check if the .bin file exists
    opts.FieldTrip       = false;    % Produce a FieldTrip-ready .mat file
    opts.doNWB           = false;    % INTAN-NeuroConv NWB export (testing)

    % Events
    opts.RetrieveEvents  = true;     % Extract event log from session
    opts.alignto         = {'itiOn'};% Alignment events; cell array of char vectors
    opts.trEvents        = {};       % 'Special' ITI events (treatments, tutors, etc.)
    opts.addtime         = 0;        % Padding around trial start/end in ms

    % Motion sensors
    opts.GetMotionSensors = false;   % Extract Deuteron head-direction sensor data

    % Data preprocessing
    opts.noise           = [];       % Reserved for noise-rejection parameters
    opts.lowpass         = 0;        % High boundary for low-pass (Hz). 0 = off.
    opts.lowpassFT       = 250;      % Low-pass for FieldTrip LFP stream (Hz)
    opts.highpass        = 0;        % Low boundary for high-pass (Hz). 0 = off.
    opts.linefilter      = 0;        % Line-noise notch centre frequency. 0 = off.
    opts.CAR             = 0;        % Common-average re-referencing. 0 = off.

    % Sorting & curation
    opts.kilosort        = 4;        % Default to Kilosort 4
    opts.KSchanMapFile   = '';       % Empty = linear array; set to 'chanMapXXX.mat' for custom
    opts.bombcell        = true;     % Run Bombcell QC on Kilosort output
    opts.phy             = false;    % Open Phy after sorting (blocks MATLAB)

    % ── NGL02_postPhy options ─────────────────────────────────────────────────
    opts.doSpikething    = true;     % Process single-unit/spike data
    opts.doLFPthing      = true;     % Process LFP data
    opts.offlineTrack    = false;    % Run offline video blob detection
    opts.FLIP            = false;    % Run vFLIP laminar power analysis
    opts.useTrack        = false;    % Index spiking against social-tracking events
    opts.trialparsed     = false;    % Load trial-parsed FT file (vs continuous)
    opts.artifdet        = false;    % Run LFP artifact detection and rejection
    opts.spectrogram     = false;    % Run multitaper time-frequency analysis
    opts.neurDyn.do      = false;    % Run neural-dynamics analysis

end
