function [myFiles,dataDirectory,sampleRate,numChannels,events]= Deuteron2Kilosort(pathRaw,folderSingleChannels,...
     ApplyNotchFilter,ApplyHighPassFilter,ApplyBandPassFilter,retrieveEvents,dllFolder)
%
% This function is a dependency of the script pipelineNeurData2Kilosort,
% only necessary if the recording system in use is Deuteron
% It compiles the data save in DT2 files in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
% Filtering and the elimination of experimentally irrelevant recording 
% periods only seem to be necessary when recording with Deuteron.
%
% DEPENDENCIES:
% *GetMetaData :  detects logger file type and builds struct with essential 
%                information for data processing (this function is 
%                downloadable through Deuteron website)
%
% *importEventsDeuteronWithDLL: interaction between Matlab and Deuteron's .NET
%                 app for event reading. Retrieves the event codes in 
%                 order, with corresponding time points (as minutes from midnight),
%                 sample number, and pins with detected rising edge
%
% *highPassFilter500,  *bandPassfilter_C1_500_C2_7500 and *notchFilter0_500
%       These functions filter the data of each channel. The filtered 
%       channels are also saved, without deleting the original ones.
%       If filtered channels exist, the final matrix is compiled with them. 
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    May 2022
%
% 20.09.2022 :    Sara    : version to use in a common pipeline with Intan and Deuteron

%% ======================== Total session Data, divided per channels
% corresponding files: get all DT2 files in folder with raw data
rawDataFiles = dir([pathRaw '\*.DT2']);
sessionFiles = {rawDataFiles.name};

% Get meta data from Deuteron:
% Checks which type of logger was used and set respective parameter:
% This must be the original extension of the file.
ext          = sessionFiles{1,1}(end-2:end);
metaData     = GetMetaData(ext);
numChannels  = metaData.numChannels;
sampleRate   = metaData.fSample;
HDF5chunkSize = 300*sampleRate;
% Create an empty file per channel:
% space is pre-allocated to save every sample of neural data (HDF5 files with infinite slots)
for i = 1:numChannels
    IndChnl  = ['Channel_', sprintf( ['%0',num2str(numel(num2str(numChannels))),'d'], i )];
    fileName = fullfile(folderSingleChannels,IndChnl);
    h5create([fileName '.h5'],['/channel_' num2str(i)],[1 Inf],'ChunkSize',[1 HDF5chunkSize],'Datatype','int16')
end
allFileNames = ls(folderSingleChannels); % list with all pre-allocated files
allFileNames(1:2,:) = [];

%% ======================== Allocate data to its respective single-channel file
% open each neural data file
% resize data for detection with Kilosort
% allocate data to its respective single-channel file
indexPos = 0;
for i = 1:length(sessionFiles)
    fileName = sessionFiles{i};             % open each file from the recording session
    myFile = fullfile(pathRaw, fileName);   % open a file from the recording session
    
    fid = fopen(myFile);
    data = fread(fid, 'uint16'); % each data point of neural data is a 16 bit word
    fclose(fid);
    
    dataMatrix = reshape(data', metaData.numChannels, []); % data are now in form of channels x samples
    KSRawdata = int16(int32(dataMatrix) - int32(intmax('uint16')/2)); % resize for detection with Kilosort
    
    % distribute each row of data to its respective single-channel file
    
    for b = 1:numChannels
        h5write(fullfile(folderSingleChannels,allFileNames(b,:)),...
            ['/channel_' num2str(b)],...	% dataset name
            KSRawdata(b,:),...              % dataset: data of a channel stored in the DT2 file (already scaled)
            [1 indexPos+1],...              
            [1 size(KSRawdata,2)])
    end
    
    indexPos = indexPos+size(KSRawdata,2);
    clear KSRawdata dataMatrix
end
%% ======================== Set filters ---> This section will be improved in a near future
%set filters
if ApplyHighPassFilter
    HighPassfilterSettings = designfilt('highpassfir','FilterOrder',20,'CutoffFrequency',500,'SampleRate',32000);
end
if ApplyBandPassFilter
    % bandPassfilterSettings = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',500,...
    %     'CutoffFrequency2',7500,'SampleRate',32000);
    bandPassfilterSettings = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',500,...
        'CutoffFrequency2',3500,'SampleRate',32000); % the biopsy version
end
if  ApplyNotchFilter
    bndWLowConfirmFilter = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',501, ...
        'DesignMethod','butter','SampleRate',32000);
end
%% ======================== Apply filters
if ApplyHighPassFilter == true | ApplyNotchFilter == true | ApplyBandPassFilter == true
    % path of the final single channel files
    filePathProssFilt = fullfile(pathProc,sessionDate,recordingSystem,'oneFilePerChannelFiltered');
    if ~exist(filePathProssFilt, 'dir')
        mkdir(pathProc,fullfile(sessionDate,recordingSystem,'oneFilePerChannelFiltered')); % create folder for filtered channels
    end
    if ApplyHighPassFilter == true
        highPassFilter500(numChannels,folderSingleChannels,filePathProssFilt,HDF5chunkSize,HighPassfilterSettings);
    elseif ApplyBandPassFilter == true
        bandPassfilter_C1_500_C2_7500(numChannels,folderSingleChannels,filePathProssFilt,HDF5chunkSize,bandPassfilterSettings);
    else % ApplyNotchFilter == true
        notchFilter0_500(numChannels,folderSingleChannels,filePathProssFilt,HDF5chunkSize,bndWLowConfirmFilter);
    end
end
%% ======================== Files for matrix compilation

if ApplyHighPassFilter == true | ApplyNotchFilter == true | ApplyBandPassFilter == true
    myFiles       = ls(fullfile(filePathProssFilt,'*.h5')); % files that will be used to compile the final matrix
    dataDirectory = fullfile(filePathProssFilt);
else
    myFiles       = ls(fullfile(folderSingleChannels,'*.h5'));
    dataDirectory = fullfile(folderSingleChannels);
end
%% ======================== If events are used to crop the matrix  
if retrieveEvents == true
    events = importEventsDeuteronWithDLL(dllFolder,pathRaw,sampleRate);
else 
    events = false; 
end
end