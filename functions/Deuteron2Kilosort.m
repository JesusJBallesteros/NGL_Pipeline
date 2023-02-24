function out = Deuteron2Kilosort(opt, sessions, ss)
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
% Last Change:    24.02.2023 (Jesus)
 
%% Get already existing Parameters
Files           = sessions.info{ss}.files;
ext             = sessions.info{ss}.fileformat;
numberOfAdcBits = sessions.info{ss}.numADCBits;
numChannels     = sessions.info{ss}.numChannels;
numFiles        = length(sessions.info{ss}.files);
HDF5chunkSize   = sessions.info{ss}.HDF5chunkSize;

param.offset            = 2^(numberOfAdcBits-1);
param.voltageResolution = 1.95e-7;

%% Create folders in case they don't exist.
mkdir(opt.FolderSingleChannels);	    % create folder for single channel files
mkdir(opt.FolderProcDataMat);        % create folder for data matrix

%% Total session Data, divided per channels
% Create an empty file per channel:
% space is pre-allocated to save every sample of neural data (HDF5 files with infinite slots)
for i = 1:numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(opt.FolderSingleChannels,IndChnl);

    if ~isfile([fileName,'.h5'])
        % Create appropiate h5 file per channel
        if ~strcmp(ext, 'DF1'), out.datatype = 'int16';
        else,                   out.datatype = 'single';
        end

        h5create([fileName '.h5'], ...
                 ['/channel_' num2str(i)], ...
                 [1 Inf], 'ChunkSize', ...
                 [1 HDF5chunkSize], ...
                 'Datatype', out.datatype);
    end
end

% list with all pre-allocated files
allFileNames = ls([opt.FolderSingleChannels,'\Ch*']); 

%% Open each neural data file, resize data for detection with Kilosort,
% Allocate data to its respective single-channel file.
% Differentiate between old and new Deuteron Formats
if ~strcmp(ext, 'DF1')
    disp('Format is FLAT. Deprecating.')
    indexPos = 0;
    for i = 1:numFiles
        % Neural data points are 16 bit words
        fid = fopen(fullfile(opt.PathRaw, Files(i).name));
            data = fread(fid, 'uint16');
        fclose(fid);

        % Shape as channels x samples
        data = reshape(data', numChannels, []); 
        KSRawdata = int16(int32(data) - int32(intmax('uint16')/2)); % resize for detection with Kilosort
        clear data

        % distribute each row of data to its respective single-channel file
        for b = 1:numChannels
            h5write(fullfile(opt.FolderSingleChannels,allFileNames(b,:)), ...
                ['/channel_' num2str(b)], ...	 % dataset name
                KSRawdata(b,:), ...              % dataset: data of a channel stored in the DT2 file (already scaled)
                [1 indexPos+1], ...              
                [1 size(KSRawdata,2)])
        end
        indexPos = indexPos+size(KSRawdata,2);

    end
    
    %% Set and apply filters. 
    % TODO: This section will be improved in a near future
%     if in.ApplyHighPassFilter
%         HighPassfilterSettings = designfilt('highpassfir','FilterOrder',3,'CutoffFrequency',350,'SampleRate',sampleRate);
%     end
%     
%     if in.ApplyHighPassFilter == true
%         % path of the final single channel files
%         filePathProssFilt = fullfile(in.folderSingleChannels,'Filtered');
%         mkdir(filePathProssFilt); % create folder for filtered channels
%     
%         Deuteron_highPassFilter(numChannels,in.folderSingleChannels,filePathProssFilt,HDF5chunkSize,HighPassfilterSettings);
%     end
    
    %% Files for matrix compilation
    out.myFiles       = ls(fullfile(opt.FolderSingleChannels,'*.h5'));
    out.dataDirectory = fullfile(opt.FolderSingleChannels);
    
    %% If events are used to crop the matrix  
    if opt.RetrieveEvents == true
    %         out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
    else 
        out.events = false; 
    end

else % ext = DF1
    %% Allocate data to its respective single-channel file
    disp('Format is BLOCK. NEW')

    % Open each neural data file, resize data for detection with Kilosort,
    % Allocate data to its respective single-channel file.
    stream      = 1;
    tmpdata     = [];
    for i = 1:numFiles         
        if ~strcmp(Files(i).name(1:4),'NEUR')
            % Skips Event files (do not contain data)
            continue
        else
            fid = fopen(fullfile(opt.PathRaw, Files(i).name), 'r');
            data = Deuteron_extractData(stream, fid, param);
            fclose(fid);

            tmpdata = [tmpdata data.neuralData];
        end
    end
    
    %% Sort by channels.
    neuralDataMat = reshape(tmpdata, numChannels, []);
    clear tmpdata

    %% Set and apply filters. 
    % TODO: This section will be improved in a near future
    %     if in.ApplyHighPassFilter
    %         HighPassfilterSettings = designfilt('highpassfir','FilterOrder',3,'CutoffFrequency',350,'SampleRate',sampleRate);
    %     end
    %     
    %     if in.ApplyHighPassFilter == true
    %         % path of the final single channel files
    %         filePathProssFilt = fullfile(in.folderSingleChannels,'Filtered');
    %         mkdir(filePathProssFilt); % create folder for filtered channels
    %     
    %         Deuteron_highPassFilter(numChannels,in.folderSingleChannels,filePathProssFilt,HDF5chunkSize,HighPassfilterSettings);
    %     end

    %% Files for matrix compilation
    if opt.ApplyHighPassFilter == true 
    %    out.myFiles       = ls(fullfile(filePathProssFilt,'*.h5')); % files that will be used to compile the final matrix
        out.dataDirectory = fullfile(filePathProssFilt);
    else
        out.myFiles       = ls(fullfile(opt.FolderSingleChannels,'*.h5'));
        out.dataDirectory = fullfile(opt.FolderSingleChannels);
    end
        
    %% If events are used to crop the matrix  
    %     if in.retrieveEvents == true
    %         out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
    %     else 
    %         out.events = false;
    %     end

    % distribute each row of data to its respective single-channel file
    for b = 1:numChannels
        h5write(fullfile(opt.FolderSingleChannels,allFileNames(b,:)), ...
            ['/channel_' num2str(b)], ...	% dataset name
            neuralDataMat(b,:), ...             % dataset: data of a channel stored in the DT2 file (already scaled)
            [1 1], ...              
            [1 size(neuralDataMat,2)])
    end
    clear neuralDataMat
end

%% Pre-define data matrix and subtract averaged matrix.
disp('Single files created. Creating Full binary and h5 files.')

% Create HDF5 file to compile data matrix
if opt.h5
    h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
             '/allChnMat', [numChannels Inf], ...
             'ChunkSize', [1 HDF5chunkSize], ...
             'Datatype', out.datatype);
end

% binary files process faster if the data is appended to it in chunks  
if opt.bin
    fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 
end

% Find out total number of samples per channel by reading .h5 file metadata
out.h5info = h5info(fullfile(out.dataDirectory, out.myFiles(1,:)),'/channel_1');
out.maxSz  = out.h5info.Dataspace.Size(2);

% Adjustments to the writing of the matrix
if opt.RetrieveEvents
    % TODO: ??
    %     input.startEvent   = events{1,2}; % first trial
    %     input.endEvent     = events{end,2};
    %     output.ChunkStart  = input.startEvent:input.stpSz:output.maxSz-mod(output.maxSz,input.stpSz);
    %     output.ChunkStart(output.ChunkStart>input.endEvent) = [];
    disp('To retrieve events is still not functional here ...');
    opt.StartEvent = 1;
    out.ChunkStart = opt.StartEvent:opt.StpSz:out.maxSz-mod(out.maxSz,opt.StpSz);
else
    opt.StartEvent = 1;
    out.ChunkStart = opt.StartEvent:opt.StpSz:out.maxSz-mod(out.maxSz,opt.StpSz);
end
      
% Fill both matrixes
% Write each channel stepwise into a matrix (hdf5) and into a binary file
disp('Filling full matrices from single channel files...')
for j = opt.StartEvent:opt.StpSz:out.ChunkStart(end)
    sngChn = cell(numChannels,1);

    for i = 1:numChannels
        sngChn{i,1} = h5read(fullfile(opt.FolderSingleChannels, out.myFiles(i,:)), ...
                            ['/channel_' num2str(i)], [1 j], [1 opt.StpSz]);
        if opt.h5
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
                '/allChnMat', sngChn{i,1}, [i j-(out.ChunkStart(1)-1)], [1 opt.StpSz]);
        end

    end
    
    if opt.bin
        % Write bin file
        fwrite(fidDataMat, cell2mat(sngChn), 'int16');
    end
end
    
fclose(fidDataMat);

if ~opt.h5
    delete(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']));
end

if ~opt.bin
    delete(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.bin']));
end

D2K.in = opt;
D2K.out = out;
save('D2K.mat','D2K','-mat');

end