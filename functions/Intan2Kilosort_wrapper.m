function I2K = Intan2Kilosort_wrapper(sessions, ss, varargin)
% Adaptation from the common pipeline for Intan. 
% Prepares recorded data for spike sorting with Kilosort.
% Uses the high-pass files to create the .h5 and .bin files
%
% DEPENDENCIES
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    input:    struct. optional inputs to override the defaults:
%              retrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                       (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                        experiment, removing paradigm irrelevant periods) 
%              stpSz: int. relative to HDF5file: chunks in which ... 
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    27.02.2023 (Jesus)

if nargin < 3, opt = struct();
elseif nargin == 3, opt = varargin{1};
end

%% Defaults
if ~isfield(opt,'RetrieveEvents'), opt.RetrieveEvents = false;    end
if ~isfield(opt,'StpSz'),          opt.StpSz          = 1000000;  end
if ~isfield(opt,'h5'),             opt.h5             = true;     end
if ~isfield(opt,'bin'),            opt.bin            = true;     end

% Paths and naming
if ~isfield(opt,'PathRaw'),              opt.PathRaw              = pwd;                               end
if ~isfield(opt,'FolderSingleChannels'), opt.FolderSingleChannels = fullfile(pwd,'oneFilePerChannel'); end
if ~isfield(opt,'FolderProcDataMat'),    opt.FolderProcDataMat    = fullfile(pwd,'processed');         end
if ~isfield(opt,'SavFileName'),          opt.SavFileName          = [sessions.list(ss).name];          end

% Create folders in case they don't exist.
mkdir(opt.FolderSingleChannels);	    % create folder for single channel files
mkdir(opt.FolderProcDataMat);           % create folder for data matrix

%% Prepare files for matrix compilation. set final parameters 
% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
%  Therefore, the matrix should be compiled with the channels in an increasing order.
% ChunkSize of HDF5 file (e.g., 5 minutes is 300s *numChannels000Hz = 9600000)
%  this chunk size works well. There is still a question of whether this
%  value is optimal. Once it is, this variable no longer requires user input.
opt.myFiles           = dir(fullfile(opt.PathRaw, 'high*.dat'));
opt.numChannels       = sessions.info{ss}.nchannels; 
opt.sampleRate        = sessions.info{ss}.amplifier_sample_rate;
opt.HDF5chunkSize     = 300*opt.sampleRate; 
opt.channelOrder      = 1:1:opt.numChannels; 

%% MAIN CALL
disp('Generating single channel .h5 files from INTAN...')
if strcmp(sessions.info{ss}.fileformat,'filepertype')
    % converts Deuteron DT2 files into the 'oneFilePerChannel' format.
    % Also filters the data and retrieve event codes if requested. 
    out = Intan2Kilosort_filepertype(opt);

elseif strcmp(sessions.info{ss}.fileformat,'fileperch')
    % TODO: implement the new format conversion from Sara
    % converts Deuteron DF1 files into the 'oneFilePerChannel' format.
    % Also filters the data and retrieve event codes if requested. 
    out = Intan2Kilosort_fileperch(opt);
end

%% Pre-define data matrix and average subtracted matrix
% create HDF5 file to compile data matrix
h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
        '/allChnMat', [opt.numChannels Inf], ...
        'ChunkSize', [1 opt.HDF5chunkSize], ...
        'Datatype', 'int16')

% total number of samples per channel 
% out.maxSz = in.myFiles(in.channelOrder(1,1),1).bytes/(in.numFilesPerChannel*2); 
out.h5info = h5info(fullfile(out.dataDirectory, out.myFiles(1,:)),'/channel_1');
out.maxSz   = out.h5info.Dataspace.Size(2);

%% Fill both matrixes
disp('Filling full matrices from single channel files...')

% binary files process faster if the data is appended to it in chunks  
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 

% Write each channel stepwise into a matrix (hdf5) and into a binary file
out.ChunkStart = 1:opt.StpSz:out.maxSz-mod(out.maxSz,opt.StpSz);

for j = 1:opt.StpSz:out.ChunkStart(end)
    sngChn = cell(opt.numChannels,1);

    for i = 1:opt.numChannels
        sngChn{i,1} = h5read(fullfile(opt.FolderSingleChannels, out.myFiles(i,:)), ...
                         ['/channel_' num2str(i)], [1 j], [1 opt.StpSz]);
        
        if opt.h5
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
                     '/allChnMat', sngChn{i,1}, [i j-(out.ChunkStart(1)-1)], [1 opt.StpSz]);
        end

    end

    fwrite(fidDataMat, cell2mat(sngChn), 'int16');
end

fclose(fidDataMat);

if ~opt.h5
    delete(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']));
end

if ~opt.bin
    delete(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.bin']));
end

% Book keeping
I2K.in = opt;
I2K.out = out;
save('I2K.mat','I2K','-mat');

end