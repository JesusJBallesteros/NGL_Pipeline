%% NGL_RunFixture. Run the pipeline over the synthetic test subject MRX.
%
% PURPOSE:
%   The fixture is analysed exactly like a real study: same stages, same opt,
%   same functions. This is the NGL_SetAndRunMe a user would write for it -
%   copy it next to your own, run it whenever you want to know whether the
%   pipeline still returns the planted answers.
%
%   It is a SEPARATE STUDY, not a subject added to yours. Two reasons:
%     * study-level outputs (aggregates, group plots) are named by analysis,
%       not by subject - a fixture run inside your study would write over them;
%     * synthetic data must never reach a real result. set_default enforces
%       this (it skips synthetic subjects when subjects = 'all', and refuses
%       to run them beside real ones), but the enforcement is a backstop. The
%       separation is the actual protection.
%
% USAGE:
%   1. Edit `datadrive` and `studyname` below to a location of your choice.
%      Use a studyname you will recognise as not-real, e.g. 'PIPELINE_TEST'.
%   2. Run section by section, exactly as you would NGL_SetAndRunMe.
%   3. `checkFixture` at the end says whether the answers still come back.
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

%% 1) PREPARE. A study of its own, holding nothing but the fixture.
clear all %#ok<CLALL>

datadrive = 'D';                 % drive letter for the test study
studyname = 'PIPELINE_TEST';     % NOT your real study name

readmecontent = ["Study name: PIPELINE_TEST", ...
                 "Synthetic data only. Subject MRX is generated, not recorded.", ...
                 "Built by tests/fixture/makeFakeStudy.m; see docs/test_fixture.md.", ...
                 "Nothing here is a measurement. Do not publish, do not aggregate", ...
                 "with real subjects."];

NGL00_Prep

%% 1b) BUILD THE DATA (first run only; skip it afterwards).
% The repository carries the generator, not the data: a session is ~270 MB and
% the seed reproduces it byte for byte. This writes the same tree NGL01 would
% have written, so every later stage has what it expects.
addpath(fullfile(input.toolbox, 'tests', 'fixture'));
makeFakeStudy(fullfile([datadrive ':'], studyname));

%% 2) SET. The fixture, and options that suit it.
subjects = {'MRX'};                                    % never 'all' here
dates    = {'19850214', '19860704', '19870321'};       % pre-1990 on purpose

areas = {'NCL','STR'};      % the two shanks of the fixture's channel map

opt = struct();
    % B.1 ACQUISITION. The data is already preprocessed, so NGL01's own work
    % is off: there is no raw file to convert and no sorter to run.
    opt.numChannels   = 32;
    opt.KSchanMapFile = 'chanMap_ATLAS-E32+R-50-S2-L10-200NT_INTAN.mat';
    opt.kilosort      = false;
    opt.bombcell      = false;
    opt.phy           = false;
    opt.doNWB         = false;
    opt.bin           = false;

    % B.2 EVENTS / TRIALS. Already written by the generator.
    opt.RetrieveEvents = false;
    opt.alignto        = {'itiOn', 'stimOn1', 'stimOn2'};

    % B.3 LFP. One band covering what is planted (6 Hz theta, 20 Hz beta),
    % and a window a second wider than the effects on either side - which the
    % fixture can afford because its trial-end event sits at a fixed maximum.
    opt.TFRmethod    = 'wavelet';
    opt.freqInterest = {4:1:40};
    opt.toi          = [-2 1.8];
    opt.timeResol    = 0.05;

    opt.lfp.session.contrast = true;
    opt.lfp.contrast.pairs   = {'correct vs incorrect'};
    opt.lfp.contrast.areas   = {'NCL', 'STR'};   % STR is the negative control
    opt.lfp.session.csd      = true;

    % B.4 FREQUENCY TAGGING. opt.nft.blocks is one set of windows for the whole
    % run, but each session's tagging block starts at a slightly different
    % second - the trials before it are drawn per session. fixtureNFTblocks
    % returns the window common to every session requested, so the same setting
    % is correct in all of them.
    opt.nft.do     = true;
    opt.nft.blocks = fixtureNFTblocks(fullfile([datadrive ':'], studyname), ...
                                      'sessions', dates);
    opt.nft.areas  = {'NCL', 'STR'};      % STR is the negative control again

param = struct();
    param.binSize  = 50;        % ms, PSTH bin
    param.stepSz   = 10;        % ms, PSTH step
    param.interval = [-500 500];

%% 3) RUN. The stages that have something to do on pre-made data.
% NGL01_Main         % not needed: the generator wrote what it produces
% NGL02_postPhy      % not needed: spike.mat is already in its post-Phy shape
NGL02_LFP            % LFP quick-look
NGL07_LFPanalysis    % contrasts + CSD
NGL08_NFT            % frequency tagging (the nft block)

%% 4) CHECK. Did the analyses return what was planted?
% This is the part a real study has no equivalent of, and the reason the
% fixture exists. It errors if any check fails.
checkFixture(fullfile([datadrive ':'], studyname));
