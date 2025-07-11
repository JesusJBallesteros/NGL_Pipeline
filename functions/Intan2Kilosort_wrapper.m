function Intan2Kilosort_wrapper(sessions, varargin)
% Adaptation from the common pipeline for Intan. Prepares recorded data in 
% the high pass for spike sorting with Kilosort. Uses the high-pass files 
% from INTAN to create .h5 and .bin files. It reads the INTAN file, either 
% a file per channel or a file for the whole bunch. Then, converts the ADC 
% step values to microvolts by multiplying by 0.195. The data comes out as 
% ch x samples in int16 format, ready for Kilosort. 
% This data is saved channel by channel and in chunks to a .h5 file. 
% This data is saved as a whole into a .bin file. 
%
% DEPENDENCIES: 
%   Intan2Kilosort_filepertype
%   Intan2Kilosort_fileperchannel
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    opt:      struct. optional inputs to override the defaults:
%               opt.StpSz          = 1000000;  int that determines the chunk size to writo into the .h5 file
%               opt.RetrieveEvents = false;    Logic that determines if we want to retrieve events.
%               opt.highpass       = [450 5000]; Array of [lowest highest] ends for the band-pass filter, in Hz 
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author: Aylin, Lukas & Sara
% Jesus 13.04.2023

if nargin < 2, opt = struct();
elseif nargin == 2, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'highpass'),       opt.highpass      = 400; end
if ~isfield(opt,'StpSz'),          opt.StpSz          = 1000000;    end

%% Collect parameters that not need to necessarily defaulted to a given value. 
% To proceed, list all files (multiple or single, depending on filetype).
% For highpass files, no further filtering is necessary (In priciple! make
% sure you are recording a proper, useful, highpass within INTAN).
opt.myFiles = dir('high*.dat');
opt.set_filter = 0;

% If no highpass files are found, it will use the raw 'amp' data and filtering will be applied.
if isempty(opt.myFiles)
    opt.myFiles = dir('amp*.dat');
    opt.set_filter = 1;
end

% How many channels, from Intan_hdr.
opt.numChannels = sessions.info.nChannels;

% Check that numChannels and the number of files coincide. Fix if necessary
if sessions.info.nfiles < opt.numChannels
    % Report
    disp('Found less channel files than expected by header. Using number of files as truth.')
    opt.numChannels = sessions.info.nfiles;
    sessions.info.nChannels = opt.numChannels;
end

% Sample rate, from Intan_hdr.
opt.sampleRate  = sessions.info.amplifier_sample_rate;

% ChunkSize of HDF5 file (e.g. 5 minutes = 300 s @30000 Hz = 9600000 samples).
opt.HDF5chunkSize = 300*opt.sampleRate;

% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
% Therefore, the matrix should be compiled with the channels in an increasing order.
% opt.channelOrder = 1:1:opt.numChannels; % Deprecating

% To obtain the number of samples per file, first read file info. Can be a 
% file per channel or only one for all, it does not matter. Read the first.
fileinfo = dir(opt.myFiles(1).name);

%% Main call
if strcmp(sessions.info.fileformat,'filepertype')
    
    % To get the number of samples, divide the file size by number of 
    % channels, times bytes that each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/(opt.numChannels * 2); 
    
    Intan2Kilosort_filepertype(opt);

elseif strcmp(sessions.info.fileformat,'fileperch')

    % To get the number of samples, divide the file size by the bytes 
    % each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/2;

    Intan2Kilosort_fileperch(opt);
    
end

end