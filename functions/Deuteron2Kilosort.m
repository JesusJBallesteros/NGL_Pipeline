function out = Deuteron2Kilosort(in, sessions)
% This function is a dependency of the script Deuteron2Kilosort_wrapper,
% only necessary if the recording system in use is Deuteron
% It compiles the data save in DT2 files in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
% Filtering and the elimination of experimentally irrelevant recording 
% periods only seem to be necessary when recording with Deuteron.
%
% DEPENDENCIES:
%  importEventsDeuteronWithDLL: interaction between Matlab and Deuteron's .NET
%                 app for event reading. Retrieves the event codes in 
%                 order, with corresponding time points (as minutes from midnight),
%                 sample number, and pins with detected rising edge
%
%  highPassFilter: These functions filter the data of each channel. The filtered 
%                   channels are also saved, without deleting the original ones.
%                   If filtered channels exist, the final matrix is compiled with them. 
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    Ene 2023 (Jesus)
 
global ss

%% Get already existing Parameters
numChannels     = sessions.info{ss}.numChannels;
HDF5chunkSize   = sessions.info{ss}.HDF5chunkSize;
Files           = sessions.info{ss}.files;
numFiles        = length(sessions.info{ss}.files);
sampleRate      = sessions.info{ss}.sampleRate;
ext             = sessions.info{ss}.fileformat;

%% Total session Data, divided per channels
% Create an empty file per channel:
% space is pre-allocated to save every sample of neural data (HDF5 files with infinite slots)
for i = 1:numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(in.folderSingleChannels,IndChnl);

    if ~isfile([fileName,'.h5'])
        % Create appropiate h5 file per channel
        if ~strcmp(ext, 'DF1'), out.datatype = 'int16';
        else,                   out.datatype = 'int16';
        end

        h5create([fileName '.h5'], ...
                 ['/channel_' num2str(i)], ...
                 [1 Inf], 'ChunkSize', ...
                 [1 HDF5chunkSize], ...
                 'Datatype', out.datatype);
    end
end

% list with all pre-allocated files
allFileNames = ls([in.folderSingleChannels,'\Ch*']); 

if ~strcmp(ext, 'DF1')
    %% Allocate data to its respective single-channel file
    % open each neural data file, resize data for detection with Kilosort,
    % allocate data to its respective single-channel file.
    indexPos = 0;
    for i = 1:numFiles
        fid = fopen(fullfile(in.pathRaw, Files(i).name));
        data = fread(fid, 'uint16'); % each data point of neural data is a 16 bit word
        fclose(fid);
        
        dataMatrix = reshape(data', numChannels, []); % data are now in form of channels x samples
%         KSRawdata = int16(int32(dataMatrix) - int32(intmax('uint16')/2)); % resize for detection with Kilosort
        KSRawdata = dataMatrix;

        % distribute each row of data to its respective single-channel file
        for b = 1:numChannels
            h5write(fullfile(in.folderSingleChannels,allFileNames(b,:)), ...
                ['/channel_' num2str(b)], ...	 % dataset name
                KSRawdata(b,:), ...              % dataset: data of a channel stored in the DT2 file (already scaled)
                [1 indexPos+1], ...              
                [1 size(KSRawdata,2)])
        end
        
        indexPos = indexPos+size(KSRawdata,2);
        clear dataMatrix
    end
    
    %% Set and apply filters. 
    % TODO: This section will be improved in a near future
    if in.ApplyHighPassFilter
        HighPassfilterSettings = designfilt('highpassfir','FilterOrder',3,'CutoffFrequency',350,'SampleRate',sampleRate);
    end
    
    if in.ApplyHighPassFilter == true
        % path of the final single channel files
        filePathProssFilt = fullfile(in.folderSingleChannels,'Filtered');
        mkdir(filePathProssFilt); % create folder for filtered channels
    
        Deuteron_highPassFilter(numChannels,in.folderSingleChannels,filePathProssFilt,HDF5chunkSize,HighPassfilterSettings);
    end
    
    %% Files for matrix compilation
    if in.ApplyHighPassFilter == true 
        out.myFiles       = ls(fullfile(filePathProssFilt,'*.h5')); % files that will be used to compile the final matrix
        out.dataDirectory = fullfile(filePathProssFilt);
    else
        out.myFiles       = ls(fullfile(in.folderSingleChannels,'*.h5'));
        out.dataDirectory = fullfile(in.folderSingleChannels);
    end
    
    %% If events are used to crop the matrix  
    if in.retrieveEvents == true
        out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
    else 
        out.events = false; 
    end

else % ext = DF1
    %% Allocate data to its respective single-channel file
    % open each neural data file, resize data for detection with Kilosort,
    % allocate data to its respective single-channel file.
    indexPos = 0;

    % Get metadata
%     numberOfAdcBits = sessions.info{ss}.numADCBits;
%     offset = 2 ^ (numberOfAdcBits - 1);
%     voltageResolution = sessions.info{ss}.voltageRes;

    for i = 2:numFiles % Skips EVENTS file (do not contains neural data)
        fid = fopen(fullfile(in.pathRaw, Files(i).name), 'r');
            data = fread(fid, Inf, 'uint8=>uint8');
        fclose(fid);
        
        % Extract metadata from block header
        tmp.constId = (hex2num(HeaderConstants.HexConstId));
        tmp.constIdBytes = typecast(tmp.constId, 'uint8');
        tmp.blockStartIndices = FindDataBlockStart(data, tmp.constIdBytes);
        tmp.numberOfBlocks = length(tmp.blockStartIndices);
        tmp.startOfFirstHeader = tmp.blockStartIndices(1);
        tmp.endOfFirstHeader = tmp.startOfFirstHeader + HeaderConstants.HeaderTotalBytes;
        tmp.firstHeader = data(tmp.startOfFirstHeader:tmp.endOfFirstHeader);
        HeaderStruct = ExtractHeaderData(tmp.firstHeader);
            
        % Extract neural data from blocks. Check where each type data is in partition info
        neuralIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.NeuralData), HeaderStruct.PartitionInfo, 'un', 0)));
        dataAsBytes = ExtractDataByType(data, HeaderStruct, neuralIndex, tmp.blockStartIndices, tmp.numberOfBlocks);
        clear data

        % Convert neural data to physical units (V)
         % cast bytes to unsigned 16 bit integers and store as float
         neuralData= typecast(dataAsBytes, 'uint16');
         clear dataAsBytes

%          % Convert uint16 values to voltage values
%          for dataPointIndex = 1:length(neuralData)
%             neuralData(dataPointIndex) = voltageResolution * (neuralData(dataPointIndex) - offset); 
%          end
    
        % sort by channels
        neuralData = reshape(neuralData, numChannels, []); % data are now in form of channels x samples
        KSRawdata = int16(int32(neuralData) - int32(intmax('uint16')/2)); % resize for detection with Kilosort
%         clear neuralData

%         % get timestamps of neural data
%         timestampsNeural = GetTimestamps(HeaderStruct.Timestamp, sampleRate, size(neuralData, 2));

        % distribute each row of data to its respective single-channel file
        for b = 1:numChannels
            h5write(fullfile(in.folderSingleChannels,allFileNames(b,:)), ...
                ['/channel_' num2str(b)], ...	 % dataset name
                KSRawdata(b,:), ...             % dataset: data of a channel stored in the DT2 file (already scaled)
                [1 indexPos+1], ...              
                [1 size(KSRawdata,2)])
        end
        
        indexPos = indexPos + size(KSRawdata,2);
        clear neuralData
    end
    
    %% Set and apply filters. 
    % TODO: This section will be improved in a near future
    if in.ApplyHighPassFilter
        HighPassfilterSettings = designfilt('highpassfir','FilterOrder',3,'CutoffFrequency',350,'SampleRate',sampleRate);
    end
    
    if in.ApplyHighPassFilter == true
        % path of the final single channel files
        filePathProssFilt = fullfile(in.folderSingleChannels,'Filtered');
        mkdir(filePathProssFilt); % create folder for filtered channels
    
        Deuteron_highPassFilter(numChannels,in.folderSingleChannels,filePathProssFilt,HDF5chunkSize,HighPassfilterSettings);
    end

    %% Files for matrix compilation
    if in.ApplyHighPassFilter == true 
        out.myFiles       = ls(fullfile(filePathProssFilt,'*.h5')); % files that will be used to compile the final matrix
        out.dataDirectory = fullfile(filePathProssFilt);
    else
        out.myFiles       = ls(fullfile(in.folderSingleChannels,'*.h5'));
        out.dataDirectory = fullfile(in.folderSingleChannels);
    end
    
    %% If events are used to crop the matrix  
    if in.retrieveEvents == true
        out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
    else 
        out.events = false;
    end

end

end