function Intan2Kilosort_wrapperV2(sessions, ss, varargin)
% Adaptation from the common pipeline for Intan. Prepares recorded data in 
% the high pass for spike sorting with Kilosort. Uses the high-pass files 
% from INTAN to create .h5 and .bin files. It reads the INTAN file, either 
% a file per channel or a file for the whole bunch. Then, converts the ADC 
% step values to microvolts by multiplying by 0.195. The data comes out as 
% ch x samples in int16 format, ready for Kilosort. 
% This data is saved channel by channel and in chunks to a .h5 file. 
% This data is save as a whole into a .bin file. 
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
%               opt.FolderProcDataMat = sessions.info{ss}.savefolder; % Path to processed folder
%               opt.SavFileName       = sessions.list(ss).name; % File name
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    28.02.2023 (Jesus)

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'StpSz'),          opt.StpSz          = 1000000;  end
if ~isfield(opt,'h5'),             opt.h5             = true;     end
if ~isfield(opt,'bin'),            opt.bin            = true;     end
if ~isfield(opt,'RetrieveEvents'), opt.RetrieveEvents = false;    end

% Paths and naming
if ~isfield(opt,'FolderProcDataMat'), opt.FolderProcDataMat = sessions.info{ss}.savefolder;        end
if ~isfield(opt,'SavFileName'),       opt.SavFileName       = sessions.list(ss).name;              end

% Create proicessed folder.
mkdir(opt.FolderProcDataMat);

%% Collect parameters to proceed with file creation
% List all files (multiple or single depending on type)
opt.myFiles = dir('high*.dat');

% To obtain the number of samples per file, first read file info. Can be a 
% file per channel or only one for all, it does not matter. Read the first.
fileinfo = dir(opt.myFiles(1).name);

% How many channels
opt.numChannels = sessions.info{ss}.nchannels; 

% Sample rate
opt.sampleRate  = sessions.info{ss}.amplifier_sample_rate;

% ChunkSize of HDF5 file (e.g., 5 minutes is, 300s at 30000Hz = 9600000 samples)
%  this chunk size works well. optimal? Once it is, this variable no longer requires user input.
opt.HDF5chunkSize = 300*opt.sampleRate; 

% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
% Therefore, the matrix should be compiled with the channels in an increasing order.
opt.channelOrder      = 1:1:opt.numChannels; 

%% Main call
disp('Generating single channel .h5 files from INTAN...')
if strcmp(sessions.info{ss}.fileformat,'filepertype')

    % To get the number of samples, divide the file size by number of 
    % channels, times bytes that each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/(opt.numChannels * 2); 

    Intan2Kilosort_filepertype(opt);

elseif strcmp(sessions.info{ss}.fileformat,'fileperch')

    % To get the number of samples, divide the file size by the bytes 
    % each int16 word takes (int16 = 2 bytes).
    opt.num_samples = fileinfo.bytes/2;

    Intan2Kilosort_fileperchV2(opt);
end

end