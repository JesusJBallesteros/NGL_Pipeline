function Intan2Kilosort_wrapper(sessions, ss, varargin)
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
%   Intan2Kilosort_filepertypeV2
%   Intan2Kilosort_fileperchannelV2
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    opt:      struct. optional inputs to override the defaults:
%               opt.StpSz          = 1000000;  int that determines the chunk size to writo into the .h5 file
%               opt.h5             = true;     Logic that determines if we want to create the .h5 file.
%               opt.bin            = true;     Logic that determines if we want to create the .bin file.
%               opt.RetrieveEvents = false;    Logic that determines if we want to retrieve events.
%               opt.highpass       = [500 7500]; Array of [lowest highest] ends for the band-pass filter, in Hz 
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
%
% Version 07.03.2023 Jesus

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'bin'),            opt.bin            = true;       end
if ~isfield(opt,'RetrieveEvents'), opt.RetrieveEvents = false;      end
if ~isfield(opt,'highpass'),       opt.highpass       = [500 7500]; end
if ~isfield(opt,'StpSz'),          opt.StpSz          = 1000000;    end

if ~isfield(opt,'h5'),             opt.h5             = true;       end

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
opt.numChannels = sessions.info(ss).nchannels; 

% Sample rate, from Intan_hdr.
opt.sampleRate  = sessions.info(ss).amplifier_sample_rate;

% ChunkSize of HDF5 file (e.g. 5 minutes = 300 s @30000 Hz = 9600000 samples).
opt.HDF5chunkSize = 300*opt.sampleRate;

% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
% Therefore, the matrix should be compiled with the channels in an increasing order.
% opt.channelOrder = 1:1:opt.numChannels; % Deprecating

% To obtain the number of samples per file, first read file info. Can be a 
% file per channel or only one for all, it does not matter. Read the first.
fileinfo = dir(opt.myFiles(1).name);

%% Main call
if strcmp(sessions.info(ss).fileformat,'filepertype')
    
    % To get the number of samples, divide the file size by number of 
    % channels, times bytes that each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/(opt.numChannels * 2); 
    
    Intan2Kilosort_filepertype(opt);

elseif strcmp(sessions.info(ss).fileformat,'fileperch')

    % To get the number of samples, divide the file size by the bytes 
    % each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/2;

    Intan2Kilosort_fileperch(opt);
    
end

end