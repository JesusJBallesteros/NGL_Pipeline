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
%               createAvrgDatMat: logic. if true, another matrix (and respective binary file) are created 
%                                         with the average of all channels subtracted from every channel 
%               retrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                       (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                        experiment, removing paradigm irrelevant periods) 
%               ApplyNotchFilter: logic. Use Line-noise filter
%               ApplyHighPassFilter: logic. Use High-pass filter
%               ApplyBandPassFilter: logic. Use Band-pass filter
%               stpSz: int. relative to HDF5file: chunks in which ... 
%               input.bandfiles: str. general name for files to use. Normally 'high*.dat' but possibility for 'amp*.dat'
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    05.01.2023 (Jesus)

if nargin < 3, in = struct(); end

%% Defaults
if ~isfield(in,'createAvrgDatMat'),     in.createAvrgDatMat     = false;    end
if ~isfield(in,'retrieveEvents'),       in.retrieveEvents       = false;    end
if ~isfield(in,'stpSz'),                in.stpSz                = 1000000;  end
if ~isfield(in,'bandfiles'),            in.bandfiles            = 'high*.dat'; end
if ~isfield(in,'keeph5'),               in.keeph5               = true;     end
if ~isfield(in,'keepbin'),              in.keepbin              = true;     end

% Paths and naming
if ~isfield(in,'pathRaw'),              in.pathRaw                  = pwd;                                       end
if ~isfield(in,'folderSingleChannels'), in.folderSingleChannels     = fullfile(pwd,'oneFilePerChannel');         end
if ~isfield(in,'folderProcDataMat'),    in.folderProcDataMat        = in.pathRaw;                                end
if ~isfield(in,'folderProcDataMatAveraged'), in.folderProcDataMatAveraged = in.pathRaw;                          end 
if ~isfield(in,'savFileName'),          in.savFileName              = [sessions.list(ss).name];              end
if ~isfield(in,'savFileNameAvrg'),      in.savFileNameAvrg          = ['Averaged_',sessions.list(ss).name];  end

% Create folders in case they don't exist.
mkdir(in.folderSingleChannels);	    % create folder for single channel files
mkdir(in.folderProcDataMat);        % create folder for data matrix
if in.createAvrgDatMat
    mkdir(in.folderProcDataMatAveraged);% create folder for avereged data matrix
end

%% Prepare files for matrix compilation. set final parameters 
% Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
%  Therefore, the matrix should be compiled with the channels in an increasing order.
% ChunkSize of HDF5 file (e.g., 5 minutes is 300s *numChannels000Hz = 9600000)
%  this chunk size works well. There is still a question of whether this
%  value is optimal. Once it is, this variable no longer requires user input.
in.myFiles           = dir(fullfile(in.pathRaw, in.bandfiles));
in.numChannels       = sessions.info{ss}.nchannels; 
in.sampleRate        = sessions.info{ss}.amplifier_sample_rate;
in.HDF5chunkSize     = 300*in.sampleRate; 
in.channelOrder      = 1:1:in.numChannels; 

%% MAIN CALL
disp('Generating single channel .h5 files from INTAN...')
if strcmp(sessions.info{ss}.fileformat,'filepertype')
    % converts Deuteron DT2 files into the 'oneFilePerChannel' format.
    % Also filters the data and retrieve event codes if requested. 
    out = Intan2Kilosort_filepertype(in);

elseif strcmp(sessions.info{ss}.fileformat,'fileperch')
    % TODO: implement the new format conversion from Sara
    % converts Deuteron DF1 files into the 'oneFilePerChannel' format.
    % Also filters the data and retrieve event codes if requested. 
    out = Intan2Kilosort_fileperch(in);
end

%% Pre-define data matrix and average subtracted matrix
% create HDF5 file to compile data matrix
h5create(fullfile(in.folderProcDataMat,[in.savFileName '.h5']), '/allChnMat', [in.numChannels Inf], 'ChunkSize', [1 in.HDF5chunkSize], 'Datatype', 'int16')

if in.createAvrgDatMat
    % create HDF5 file to compile averaged data matrix
    h5create(fullfile(in.folderProcDataMatAveraged,[in.savFileNameAvrg '.h5']), '/avgSubtracted', [in.numChannels Inf], 'ChunkSize', [1 in.stpSz], 'Datatype', 'int16') 
end

% total number of samples per channel 
% out.maxSz = in.myFiles(in.channelOrder(1,1),1).bytes/(in.numFilesPerChannel*2); 
out.h5info = h5info(fullfile(out.dataDirectory, out.myFiles(1,:)),'/channel_1');
out.maxSz   = out.h5info.Dataspace.Size(2);

%% Fill both matrixes
disp('Filling full matrices from single channel files...')

% binary files process faster if the data is appended to it in chunks  
fidDataMat = fopen(fullfile(in.folderProcDataMat,[in.savFileName '.bin']), 'a'); 

if in.createAvrgDatMat
    fidDataMatAvg = fopen(fullfile(in.folderProcDataMatAveraged,[in.savFileNameAvrg '.bin']), 'a');
end

% Write each channel stepwise into a matrix (hdf5) and into a binary file
out.ChunkStart = 1:in.stpSz:out.maxSz-mod(out.maxSz,in.stpSz);

for j = 1:in.stpSz:out.ChunkStart(end)
    sngChn = cell(in.numChannels,1);

    if in.createAvrgDatMat
        % mean, to calculate averaged matrix is computed over one chunk of data at a time
        meanForSub = int16(zeros(in.numChannels,in.stpSz)); 
    end

    for i = 1:in.numChannels
        sngChn{i,1} = h5read(fullfile(in.folderSingleChannels, out.myFiles(i,:)), ['/channel_' num2str(i)], [1 j], [1 in.stpSz]);
        
        if in.keeph5
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(in.folderProcDataMat, [in.savFileName '.h5']), '/allChnMat', sngChn{i,1}, [i j-(out.ChunkStart(1)-1)], [1 in.stpSz]);
        end

        if in.createAvrgDatMat
            meanForSub(i,:) = sngChn{i,1}/in.numChannels;   
        end
    end

    if in.createAvrgDatMat
        % write into a new matrix and a new binary file the average subtracted channel data 
        subtrAverage(meanForSub,numChannels, sngChn, in.folderProcDataMatAveraged, in.savFileNameAvrg, j, in.stpSz, fidDataMatAvg, out.ChunkStart);
    end

    fwrite(fidDataMat, cell2mat(sngChn), 'int16');
end

fclose(fidDataMat);

if ~in.keeph5
    delete(fullfile(in.folderProcDataMat, [in.savFileName '.h5']));
end

if ~in.keepbin
    delete(fullfile(in.folderProcDataMat, [in.savFileName '.bin']));
end

if in.createAvrgDatMat
    fclose(fidDataMatAvg);
end

% Book keeping
I2K.in = in;
I2K.out = out;
save('I2K.mat','I2K','-mat');

end