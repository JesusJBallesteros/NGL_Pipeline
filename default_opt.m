function opts = default_opt()
% default_opt  Return the canonical default options struct for the NGL toolbox.
%
% PURPOSE:
%   Single source of truth for every configurable option. Every recognised
%   option name MUST appear here with a safe default value. Downstream
%   functions must never hard-code their own defaults — they rely on this
%   struct already being complete and validated before they are called.
%
% USAGE:
%   opts = default_opt();
%   Called internally by set_default. Users should not call this directly;
%   instead, set options in the opt struct inside NGL_SetAndRunMe.m.
%
% OUTPUT:
%   opts - struct with all NGL option fields pre-filled to safe defaults.
%          Field names are the canonical names recognised by set_default.
%
% ADDING A NEW OPTION:
%   1. Add the field here with its default value and a comment.
%   2. Add validation logic in set_default (Section 2) if needed.
%   3. Document it in wiki_NGL01_pipeline.md (Section 6).
%
% Last modified 07.05.2026 (Jesus)

    % Data format
    opts.numChannels     = 32;      % Expected channel count (override for 32-ch Deuteron)
    opts.bin             = true;    % Normally, we always check if the .bin file exists
    opts.FieldTrip       = true;    % Produce a FieldTrip-ready .mat file
    opts.doNWB           = false;    % INTAN-NeuroConv NWB export (testing)

    % Events
    opts.RetrieveEvents  = true;     % Extract event log from session
    opts.alignto         = {'itiOn'};% Alignment events; cell array of char vectors
    opts.trEvents        = {};       % 'Special' ITI events (treatments, tutors, etc.)
    opts.addtime         = 0;        % Padding around trial start/end in ms
    opts.uselog          = false;    % By default, use Deuteron data files to extract events. 
                                      % When true, uses the text log. For cases when the events 
                                      % were not properly transmitted to the system but logged.

    % Motion sensors
    opts.GetMotionSensors = false;   % Extract head-direction sensor data

    % Data preprocessing
    opts.noise           = [];       % Reserved for noise-rejection parameters
    opts.lowpass         = 9000;     % High boundary for low-pass (Hz). [] = off.
    opts.lowpassFT       = 200;      % Low-pass for FieldTrip LFP stream (Hz)
    opts.highpass        = [];       % Low boundary for high-pass (Hz). [] = off.
    opts.linefilter      = 0;        % Line-noise notch centre frequency. 0 = off.
    opts.CAR             = 0;        % Common-average re-referencing. 0 = off.
    opts.dwnsmplRate     = [];       % LFP downsample target (Hz). [] = auto (937.5 Hz).

    % Sorting & curation
    opts.kilosort        = 1;        % Default to Kilosort 4
    opts.KSchanMapFile   = '';       % Empty = linear array; set to 'chanMapXXX.mat' for custom
    opts.bombcell        = true;     % Run Bombcell QC on Kilosort output
    opts.phy             = false;    % Open Phy after sorting (blocks MATLAB)

    % NGL02_postPhy options
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
