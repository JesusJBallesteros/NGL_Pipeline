%% Jesus' Pipeline to analyze FieldTrip formatted files
%%% pipelineNeurData2Kilosort
%
%   Common pipeline for Intan and Deuteron recording systems. Prepare 
%   recorded data for spike sorting with Kilosort
%
% DEPENDENCIES
%  * read_Intan_RHD2000_file : function for reading data recording header
%  * Deuteron2Kilosort       : function to compile recorded data in a 
%                              single file per channel. Can also filter the
%                              data and retrieve event codes. 
%       - List of dependencies of this function (open function for
%       details):
%           *GetMetaData
%           *importEventsDeuteronWithDLL
%           *highPassFilter500
%           *bandPassfilter_C1_500_C2_7500
%           *notchFilter0_500
%
% INPUTS:
%        none, all variables are defined in the script
% 
% OUTPUT:
%        Binary file, channels(rows) per sample (columns), with channels
%        in increasing order ? as required for processing with Kilosort 
% 
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    23.09.2022
%
% Modified to add general scripting by Jesus 19.10.2022

%% ADD-on to fit regular scripting
% TES: {'20220609' '20220809' '20220922'}
% FRN: {'20181029' '20181105' '20181106'}
% FAT: {'20220909'}
% 427: {'20220929'}
input.animal   = '427';   % string. ie 'TES', 'FAT' , 'FRN', '427' ...
input.dates    = {'20220929'}; % cell array {'yyyymmdd' ...} or string 'all'
input.mainfolder = 'C:\Code\Scripts\ephys-data-pipeline'; % string. Main pipeline folder
input.datafolder = "D:\Experiments\";   % char string. Main data folder



%% START ORIGINAL SCRIPT 
%Notes:
% the most suitable folder structure
%   user/data/raw/bird/session/recording system
%   user/data/processed/bird/session/recording system/dataMatrix
%                                                    /dataMatrixAveraged
%                                                    /oneFilePerChannel
%                                                    /oneFilePerChannelFiltered

%% ======================== Processing parameters

% define subject and date of recording:
bird                    = 'Fatboy';  %these build the directories in later functions, so stay in folder naming conventions
sessionDate             = '20220318';
recordingSystem         = 'Deuteron'; %recording System
dataFormat              = 'oneFilePerChannel'; % format in which data is saved (Intan) or desired file format (Deuteron)
corePthRaw              = 'E:\Sara\data\raw'; %path constructor - raw data
corePthProc             = 'D:\Sara\data';%path constructor - processed data 
%Paths:
pathRaw                 = fullfile(corePthRaw,bird,sessionDate,recordingSystem);%path to folder containing the raw files
pathProc                = fullfile(corePthProc,'preprocessed',bird);%path constructor - processed data 
folderSingleChannels    = fullfile(pathProc,sessionDate,recordingSystem,'oneFilePerChannel');% path of the single channel files
folderProcDataMat       = fullfile(pathProc,sessionDate,recordingSystem,'dataMatrix');% path for processed matrix 
savFileName             = ['dataMatrix_',bird,sessionDate]; 

createAvrgDatMat        = true; % if true, another matrix (and respective binary file) are created 
                                %with the average of all channels
                                %subtracted from every channel 
stpSz                   = 1000000; % relative to HDF5file: chunks in which 
                                   %each channel is read and written into the final matrix
if createAvrgDatMat
    folderProcDataMatAveraged = fullfile(pathProc,sessionDate,recordingSystem,'dataMatrixAveraged'); % path for averaged matrix 
    savFileNameAvrg     = ['dataMatrixAveraged_',bird,sessionDate];
end

% make directories if they don't yet exist
if ~exist(folderProcDataMat, 'dir')
    mkdir(pathProc,fullfile(sessionDate,recordingSystem,'dataMatrix'));         % create folder for data matrix
end
if ~exist(folderProcDataMatAveraged, 'dir') && createAvrgDatMat
    mkdir(pathProc,fullfile(sessionDate,recordingSystem,'dataMatrixAveraged'));	% create folder for avereged data matrix
end

if strcmpi(recordingSystem,'Deuteron') % Deuteron specific parameters 
    if ~exist(folderSingleChannels, 'dir')
        mkdir(pathProc,fullfile(sessionDate,recordingSystem,'oneFilePerChannel'));	% create folder for single channel files
    end
    retrieveEvents          = false; %possibility to load event codes to bild a restriced matrix
    %(e.g., the matrix starts at the first 'itiON' and ends at 'end'
    %experiment, removing paradigm irrelevant periods) 
    dllFolder               = 'C:\Sara\Deuteron'; %forder where the the file Event_File_Reader_7_2.dll is 
    %(provided by Deuteron with the Event_File_Reader app)
    ApplyNotchFilter        = false;% Filtering
    ApplyHighPassFilter     = false;
    ApplyBandPassFilter     = false;
end

%% ======================== Prepare files for matrix compilation; set final parameters 

if strcmpi(recordingSystem,'Intan')
    fileName            = 'info.rhd';           % file name of INTAN header file (.rhd)
    read_Intan_RHD2000_file(pathRaw,fileName);  % read INTAN header file
    myFiles             = dir(fullfile(pathRaw,'*.dat')); %get files
    dataDirectory       = pathRaw; 
    numChannels         = length(amplifier_channels); 
    sampleRate          = frequency_parameters.amplifier_sample_rate;  
    numFilesPerChannel  = 1;                    % if we ever change the way we save data, this will change
    channelOrder        = 1:numChannels; % Kilosort rearranges the rows of the input matrix according to the a channel map (which is developed in another file).
% Therefore, the matrix should be compiled with the channels in an increasing order.
else
    [myFiles,dataDirectory,sampleRate,numChannels,events]= Deuteron2Kilosort(pathRaw,folderSingleChannels,...
        ApplyNotchFilter,ApplyHighPassFilter,ApplyBandPassFilter,retrieveEvents,dllFolder); 
    % converts Deuteron DT2 files into the 'oneFilePerChannel' format. Can
    % also filter the data and retrieve event codes. 
end

HDF5chunkSize           = 300*sampleRate; %chunk size of HDF5 file (e.g., 5 minutes is 300s *numChannels000Hz = 9600000)
% this chunk size works well. There is still a question of whether this
% value is optimal. Once it is, this variable no longer requires user
% input.
%% ======================== Pre-define data matrix and average subtracted matrix

h5create(fullfile(folderProcDataMat,[savFileName '.h5']),... % create HDF5 file to compile data matrix
    '/allChnMat',[numChannels Inf],'ChunkSize',[1 HDF5chunkSize],'Datatype','int16')

if createAvrgDatMat
    h5create(fullfile(folderProcDataMatAveraged,[savFileNameAvrg '.h5']),'/avgSubtracted',...
        [numChannels Inf],'ChunkSize',[1 stpSz],'Datatype','int16') % create HDF5 file to compile averaged data matrix
end
if strcmpi(recordingSystem,'Intan')
    maxSz = myFiles(channelOrder(1,1),1).bytes/(numFilesPerChannel*2); %total number of samples per channel 
else
    infoHlp = h5info(fullfile(dataDirectory,myFiles(1,:)));
    maxSz = infoHlp.Datasets.Dataspace.Size(2);%total number of samples per channel 
end

%% ======================== Adjustments to the writing of the matrix

if strcmpi(recordingSystem,'Deuteron') % For Deuteron only 
    if retrieveEvents == true
        startEvent = events{1,2}; %first trial
        endEvent = events{end,2};
        ChunkStart = startEvent:stpSz:maxSz-mod(maxSz,stpSz);
        ChunkStart(ChunkStart>endEvent)=[];
    else
        startEvent = 1;
        ChunkStart = startEvent:stpSz:maxSz-mod(maxSz,stpSz);
    end
end
%% ======================== Fill both matrixes
fidDataMat = fopen(fullfile(folderProcDataMat,[savFileName '.bin']), 'a'); %binary files process faster if the data is appended to it in chunks  
if createAvrgDatMat
    fidDataMatAvg = fopen(fullfile(folderProcDataMatAveraged,[savFileNameAvrg '.bin']), 'a');
end
if strcmpi(recordingSystem,'Intan')
    % Writing each channel into the HDF5 file in int16 (readable for kilosort)
    for i = 1 : numChannels
        fid = fopen(fullfile(pathRaw, myFiles(channelOrder(i),1).name), 'r'); % read .dat files as they are recorded by intan
        h5write(fullfile(folderProcDataMat,[savFileName '.h5']),'/allChnMat',...
            fread(fid,[1,maxSz], 'int16'),...
            [i 1],[1 maxSz]); % write each channel as a whole into the matrix (hdf5)
        fclose(fid);
    end
    % Writing HDF5 file step by step into binary
    for j=1:stpSz:maxSz-mod(maxSz,stpSz)
        if createAvrgDatMat
            meanForSub = int16(zeros(numChannels,stpSz)); % mean, to calculate averaged matrix is computed over one chunk of data at a time
        end
        Chunk = int16(zeros(numChannels,stpSz));
        for k=1:numChannels
            Chunk(k,:) = h5read(fullfile(folderProcDataMat,[savFileName '.h5']),'/allChnMat',[k j],[1 stpSz]);
            if createAvrgDatMat
                meanForSub(k,:) = Chunk(k,1)/numChannels;   
            end
        end
        fwrite(fidDataMat,Chunk,'int16');
        if createAvrgDatMat
            subtrAverage(meanForSub,numChannels,Chunk,folderProcDataMatAveraged,savFileNameAvrg,j,stpSz,fidDataMatAvg)
            % write into a new matrix and a new binary file the average subtracted channel data  
        end
    end
else %Deuteron
    %write each channel stepwise into a matrix (hdf5) and into a binary file
    for j=startEvent:stpSz:ChunkStart(end)
        if createAvrgDatMat
            meanForSub = int16(zeros(numChannels,stpSz)); % mean, to calculate averaged matrix is computed over one chunk of data at a time
        end
        sngChn = cell(numChannels,1);
        for i = 1:numChannels
            sngChn{i,1} = h5read(fullfile(dataDirectory,myFiles(i,:)),['/channel_' num2str(i)],[1 j],[1 stpSz]);
            h5write(fullfile(folderProcDataMat,[savFileName '.h5']),'/allChnMat',...
                sngChn{i,1},...
                [i j-(ChunkStart(1)-1)],[1 stpSz])                % compile channels -> unaltered matrix to load into kilosort
            if createAvrgDatMat
                meanForSub(i,:) = sngChn{i,1}/numChannels;   
            end
        end
        if createAvrgDatMat
            subtrAverage(meanForSub,numChannels,sngChn,folderProcDataMatAveraged,savFileNameAvrg,j,stpSz,fidDataMatAvg,ChunkStart)
            % write into a new matrix and a new binary file the average subtracted channel data  
        end
        fwrite(fidDataMat,cell2mat(sngChn),'int16');
    end
end
fclose(fidDataMat);
if createAvrgDatMat
    fclose(fidDataMatAvg);
end



function subtrAverage(meanForSub,numChannels,sngChn,folderProcDataMatAveraged,savFileNameAvrg,j,stpSz,fidDataMatAvg,ChunkStart)
% write into a new matrix and a new binary file the average subtracted channel data  

    meanForSub = int16(sum(meanForSub,1));%average of all channels 
    sngChnAvg = cell(numChannels,1);
    for k=1:numChannels
        if iscell(sngChn)
            sngChnAvg{k,1} = sngChn{k,1}-meanForSub; %average subtracted
        else
            sngChnAvg{k,1} = sngChn(k,1)-meanForSub;
        end 
        if nargin == 9
            h5write(fullfile(folderProcDataMatAveraged,[savFileNameAvrg '.h5']),'/avgSubtracted',...
                sngChnAvg{k,1},...
                [k j-(ChunkStart(1)-1)],[1 stpSz]) % write average subtracted data into a new matrix    
        else
            h5write(fullfile(folderProcDataMatAveraged,[savFileNameAvrg '.h5']),'/avgSubtracted',...
                sngChnAvg{k,1},[k j],[1 stpSz]) % write average subtracted data into a new matrix 
        end
    end
    fwrite(fidDataMatAvg,cell2mat(sngChnAvg),'int16'); % write average subtracted data into a new binary file
    
end
%% END Original script