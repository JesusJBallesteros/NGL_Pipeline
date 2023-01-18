function D2K = Deuteron2Kilosort_wrapper(sessions, varargin)
% Adaptation from the common pipeline for Deuteron. Prepares recorded data for spike sorting with Kilosort.
% For now, uses *.DT2 files and creates .h5 and .bin files.
%
% DEPENDENCIES
%    Deuteron2Kilosort: function to compile recorded data in a single file per channel.
%                      Can also filter the data and retrieve event codes. 
%
% INPUTS:
%    sessions: struct. Variable containing info about sessions in process
%    ss:       int. Current session ordinal in the pipeline
%    input:    struct. optional inputs to override the defaults:
%               createAvrgDatMat: logic. if true, another matrix (and respective binary file) are created 
%                                     with the average of all channels subtracted from every channel 
%               retrieveEvents: logic. possibility to load event codes to build a restriced matrix
%                                   (e.g., the matrix starts at the first 'itiON' and ends at 'end'
%                                    experiment, removing paradigm irrelevant periods) 
%               ApplyHighPassFilter: logic. Use High-pass filter
%               stpSz: int. relative to HDF5file: chunks in which ...  
%
% OUTPUT:
%    Binary file, channels(rows) per sample (columns), with channels
%    in increasing order as required for Kilosort.
%    h5 file, channels(rows) per sample (columns), with channels
%    in increasing order as required for Kilosort.
%    Also, a file named 'D2K.mat' containing used inputs and outputs.
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    18.01.2023 (Jesus)

global ss

if nargin < 2, in = struct(); end

%% Defaults
if ~isfield(in,'createAvrgDatMat'),     in.createAvrgDatMat     = false;    end
if ~isfield(in,'retrieveEvents'),       in.retrieveEvents       = false;    end
if ~isfield(in,'ApplyHighPassFilter'),  in.ApplyHighPassFilter  = false;    end
if ~isfield(in,'stpSz'),                in.stpSz                = 1000000;  end
if ~isfield(in,'keeph5'),               in.keeph5               = true;     end
if ~isfield(in,'keepbin'),              in.keepbin              = true;     end

% Paths and naming
if ~isfield(in,'pathRaw'),              in.pathRaw                  = pwd;                                                  end
if ~isfield(in,'folderSingleChannels'), in.folderSingleChannels     = fullfile(pwd,'oneFilePerChannel');                    end
if ~isfield(in,'folderProcDataMat'),    in.folderProcDataMat        = in.pathRaw;                                           end %fullfile(pwd,'dataMatrix');
if ~isfield(in,'folderProcDataMatAveraged'), in.folderProcDataMatAveraged = in.pathRaw;                                     end %fullfile(pwd,'dataMatrixAveraged');
if ~isfield(in,'savFileName'),          in.savFileName              = sessions.list(ss).name;                               end
if ~isfield(in,'savFileNameAvrg'),      in.savFileNameAvrg          = ['Averaged_',sessions.list(ss).name];                 end
if ~isfield(in,'dllFolder'),            in.dllFolder                = 'C:\Code\Scripts\ephys-data-pipeline\function\dlls';  end

% Create folders in case they don't exist.
mkdir(in.folderSingleChannels);	    % create folder for single channel files
mkdir(in.folderProcDataMat);        % create folder for data matrix
if in.createAvrgDatMat
    mkdir(in.folderProcDataMatAveraged);% create folder for avereged data matrix
end

%% MAIN CALL
% Converts Deuteron DT2 and DF1 files into the 'oneFilePerChannel' format.
% Also filters the data and retrieve event codes if requested. 
disp('Generating single channel files from Deuteron...')
out = Deuteron2Kilosort(in, sessions);

%% Get already existing Parameters
numChannels     = sessions.info{ss}.numChannels;
HDF5chunkSize   = sessions.info{ss}.HDF5chunkSize;

%% Pre-define data matrix and average subtracted matrix.
% Create HDF5 file to compile data matrix
h5create(fullfile(in.folderProcDataMat, [in.savFileName '.h5']), ...
         '/allChnMat', [numChannels Inf], ...
         'ChunkSize', [1 HDF5chunkSize], ...
         'Datatype', out.datatype)

% create HDF5 file to compile averaged data matrix
if in.createAvrgDatMat
    h5create(fullfile(in.folderProcDataMatAveraged,[in.savFileNameAvrg '.h5']), ...
             '/avgSubtracted', [numChannels Inf], ...
             'ChunkSize', [1 in.stpSz], ...
             'Datatype', out.datatype) 
end

% Find out total number of samples per channel by reading .h5 file metadata
out.h5info = h5info(fullfile(out.dataDirectory, out.myFiles(1,:)),'/channel_1');
out.maxSz   = out.h5info.Dataspace.Size(2);

%% Adjustments to the writing of the matrix
if in.retrieveEvents == true
    disp('Retrieving Events from Deuteron...')

% TODO: to test    
%     input.startEvent   = events{1,2}; % first trial
%     input.endEvent     = events{end,2};
%     output.ChunkStart          = input.startEvent:input.stpSz:output.maxSz-mod(output.maxSz,input.stpSz);
%     output.ChunkStart(output.ChunkStart>input.endEvent) = [];
     disp('To retrieve events is still not functional here ...');
     return
else
    in.startEvent = 1;
    out.ChunkStart = in.startEvent:in.stpSz:out.maxSz-mod(out.maxSz,in.stpSz);
end
    
%% Fill both matrixes
disp('Filling full matrices from single channel files...')

% binary files process faster if the data is appended to it in chunks  
fidDataMat = fopen(fullfile(in.folderProcDataMat,[in.savFileName '.bin']), 'a'); 

if in.createAvrgDatMat
    fidDataMatAvg = fopen(fullfile(in.folderProcDataMatAveraged,[in.savFileNameAvrg '.bin']), 'a');
end

% Write each channel stepwise into a matrix (hdf5) and into a binary file
for j = in.startEvent:in.stpSz:out.ChunkStart(end)
    sngChn = cell(numChannels,1);

    if in.createAvrgDatMat
        % mean, to calculate averaged matrix is computed over one chunk of data at a time
        meanForSub = int16(zeros(numChannels,in.stpSz)); 
    end

    for i = 1:numChannels
        sngChn{i,1} = h5read(fullfile(in.folderSingleChannels, out.myFiles(i,:)), ...
                            ['/channel_' num2str(i)], [1 j], [1 in.stpSz]);
        
        if in.keeph5
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(in.folderProcDataMat, [in.savFileName '.h5']), ...
                '/allChnMat', sngChn{i,1}, [i j-(out.ChunkStart(1)-1)], [1 in.stpSz]);
        end

        if in.createAvrgDatMat
            meanForSub(i,:) = sngChn{i,1}/numChannels;   
        end
    end

    if in.createAvrgDatMat
        % write into a new matrix and a new binary file the average subtracted channel data 
        subtrAverage(meanForSub,numChannels, ...
                    sngChn, in.folderProcDataMatAveraged, ...
                    in.savFileNameAvrg, j, in.stpSz, ...
                    fidDataMatAvg, out.ChunkStart);
    end
    
    % Write bin file
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

D2K.in = in;
D2K.out = out;
save('D2K.mat','D2K','-mat');

end