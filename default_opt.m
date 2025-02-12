%%  Create the structure 'opt' with default values. IN PROGRESS
% DO NOT MODIFY here, use input values in script,
% whose will run after this is executed.
opt = struct();

%% General
opt.numChannels      = 32;
opt.CAR              = 1;
opt.linefilter       = 0;
opt.bin              = true;
opt.highpass         = 400;
opt.FieldTrip        = false;
opt.lowpass          = 250;
opt.GetMotionSensors = false;
opt.RetrieveEvents   = true;
opt.alignto          = {'itiOn', 'stimOn1', 'rwd'};
opt.trEvents         = {'na1'};
opt.addtime          = 1500;
opt.kilosort         = 4;
opt.KSchanMapFile    = 'chanMapE32-S2_DeutSN11.mat';
opt.spkTh            = -6;
opt.bombcell         = false;
opt.rerun            = false;
opt.phy              = false;
%% Deuteron pipeline
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FieldTrip'),       opt.FieldTrip           = true;         end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = false;        end
if ~isfield(opt,'lowpass'),         opt.lowpass             = 250;          end
if ~isfield(opt,'highpass'),        opt.highpass            = 400;          end
if ~isfield(opt,'StpSz'),           opt.StpSz               = 1000000;      end
if ~isfield(opt,'parsetrial'),      opt.parsetrial          = false;        end
if ~isfield(opt,'CAR'),             opt.CAR                 = true;         end
%% INTAN pipeline 
if ~isfield(opt,'bin'),             opt.bin                 = true;         end
if ~isfield(opt,'FieldTrip'),       opt.FieldTrip           = true;         end
if ~isfield(opt,'useNWB'),          opt.useNWB              = false;        end
if ~isfield(opt,'RetrieveEvents'),  opt.RetrieveEvents      = true;         end
if ~isfield(opt,'GetMotionSensors'),opt.GetMotionSensors    = true;         end
if ~isfield(opt,'lowpass'),         opt.lowpass             = 150;          end
%% NGL02_postPhy
if ~isfield(opt, 'doSpikething') || isempty(opt.doSpikething),  opt.doSpikething = true;    end
if ~isfield(opt, 'doLFPthing') || isempty(opt.doLFPthing),      opt.doLFPthing   = true;    end
if ~isfield(opt, 'offlineTrack') || isempty(opt.offlineTrack),  opt.offlineTrack = false;   end
if ~isfield(opt, 'FLIP') || isempty(opt.FLIP),                  opt.FLIP         = false;   end
% if exist('regions','var'),                                      opt.multregion   = true;    end

%% README.TXT 
% It contains details about the project. File can also be modified later.
readmecontent = ["Study name: Preparation for social learning paradigms", ...
                 "Readme date: 22/02/2024"                          , ...
                 "Person (1) responsible for data repository: Jesus", ...
                 "Person(s) responsible for study: Juan, Jesus "    , ...
                 "Hardware used: Deuteron. E32-S2"                  , ...
                 "Related Publication(s): None."                    , ...
                 "Short description of study: Autoshaping with social learning." ];
% end