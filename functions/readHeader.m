function [nfo] = readHeader(varargin)

% function [nfo] = readHeader(pth)
%
% Use this function to read Intan header file (one file per channel).
%   Code of core functionality is taken from Intan documentation.
%
% INPUTS-OPTIONAL
%  * path           : Data path.
%  * verbose        : write details to command window
%  * warnMe         : warn if parameters are not met (pass code as string)
%  * noTime         : don't collect time-stamps of samples or numSmp (fast)
%
% OUTPUTS
%  * nfo            : Header information
%  * tms            : time stamps of all samples
%
% EXAMPLES
% 1.
% [nfo] = readHeader('path',datPth, 'verbose',0, 'warnme', ...
%     ['round(nfo.ampSmpRate,2) == 30000 && ' ...
%     'round(nfo.desired_lower_bandwidth,2) == 0.1 &&' ...
%     'round(nfo.desired_upper_bandwidth,2) == 5000']);
%
% See also: readingExample, readAux, readEvents, readNeural, readAdc,
%           readAmpToKlusta, generatePrm

% VERSION HISTORY:
% Author:        Jonas Rose
% Version:       1.2

% BUGS/ TODO:
% - uses Intan's function omitting some info, maybe expand

% 05.17.2016, Jonas: Release version
% 02.08.2016, Jonas: File not found error
% 29.01.2019, Jonas: fixed some minor bugs


%% initializations
verbose     = true;
pth         = [];
checkThis   = '';
warnMe      = false;
noTime      = false;

%% get optional inputs
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case 'path'
            i           = i+1;
            pth         = varargin{i};
        case 'verbose'
            i           = i+1;
            verbose     = varargin{i};
        case 'warnme'
            i           = i+1;
            checkThis   = varargin{i};
            warnMe      = true;
        case 'notime'
            i           = i+1;
            noTime      = true;
    end
    i=i+1;
end

% let the user select the path
if isempty(pth)
    pth = uigetdir;
end

if verbose
    disp(['Reading from directory: ' pth]);
end

%% read the header file
fid = fopen(fullfile(pth,'info.rhd'),'r');
if fid==-1
    error(['File not found, inorrect path? Check: ' fullfile(pth,'info.rhd')]);
end

%% taken from Intan's "read_Intan_RHD2000_file"

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
if fread(fid, 1, 'uint32') ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number [primary secondary].
nfo.version = [fread(fid, 1, 'int16') fread(fid, 1, 'int16')];

% Read information of sampling rate and amplifier frequency settings.
nfo.ampSmpRate                  = fread(fid, 1, 'single');
nfo.auxSmpRate                  = nfo.ampSmpRate / 4;
nfo.dsp_enabled                 = fread(fid, 1, 'int16');
nfo.actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
nfo.actual_lower_bandwidth      = fread(fid, 1, 'single');
nfo.actual_upper_bandwidth      = fread(fid, 1, 'single');
nfo.desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
nfo.desired_lower_bandwidth     = fread(fid, 1, 'single');
nfo.desired_upper_bandwidth     = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
nfo.notch_filter_mode = fread(fid, 1, 'int16');
nfo.notch_filter_frequency = 0;
if (nfo.notch_filter_mode == 1)
    nfo.notch_filter_frequency = 50;
elseif (nfo.notch_filter_mode == 2)
    nfo.notch_filter_frequency = 60;
end

% Close data file.
fclose(fid);

%% read the timepoints of each sample
if ~noTime
    fid          = fopen(fullfile(pth,'time.dat'),'r');
    nfo.times    = fread(fid, inf,'int32');
    fclose(fid);
    nfo.times    = nfo.times/ nfo.ampSmpRate;
    nfo.numSmp   = length(nfo.times);
end

%% show some results on screen
if verbose
    disp(['Sampling rate: ' num2str(nfo.ampSmpRate) ' Hz']);
    disp(['Bandpass: ' ...
        num2str(nfo.desired_lower_bandwidth) ' (' num2str(nfo.actual_lower_bandwidth) ') - ' ...
        num2str(nfo.desired_upper_bandwidth) ' (' num2str(nfo.actual_upper_bandwidth) ')']);
end

%% trigger a warning if specified condition is not met
if warnMe && ~eval(checkThis)
    warning(['Encountered unexpected parameters: ' ...
        'Sampling rate: ' num2str(nfo.ampSmpRate) ' Hz ' ...
        'Bandpass: ' ...
        num2str(nfo.desired_lower_bandwidth) ' (' num2str(nfo.actual_lower_bandwidth) ') - ' ...
        num2str(nfo.desired_upper_bandwidth) ' (' num2str(nfo.actual_upper_bandwidth) ')']);
end