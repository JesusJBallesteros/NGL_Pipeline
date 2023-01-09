function out = Intan2Kilosort_fileperch(in)
% This function is a dependency of the script Intan2Kilosort_wrapper,
% only necessary if the recording system in use is Intan
% It compiles the data save as filepertype format in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    Ene 2023 (Jesus)
 
%% Total session Data, divided per channels
% Create an empty file per channel:
% space is pre-allocated to save every sample of neural data (HDF5 files with infinite slots)
for i = 1:in.numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(in.folderSingleChannels, IndChnl);

    if ~isfile([fileName,'.h5'])
        h5create([fileName '.h5'],['/channel_' num2str(i)],[1 Inf],'ChunkSize',[1 in.HDF5chunkSize],'Datatype','int16');
    end
end

% list with all pre-allocated files
allFileNames = ls([in.folderSingleChannels,'\Ch*']); 

%% Allocate data to its respective single-channel file
% open each neural data file, resize data for detection with Kilosort,
% allocate data to its respective single-channel file.
for i = 1:length(in.myFiles)
    fid = fopen(fullfile(in.pathRaw, in.myFiles(i).name));
    data = fread(fid, [1 inf], 'int16'); % each data point of neural data is a 16 bit word
    fclose(fid);
    
    % data in form of channels x samples
    data = reshape(data', 1, []);

    % resize for Kilosort
    KSRawdata = int16(int32(data) - int32(intmax('int16')/2));
    
    % distribute each row of data to its respective single-channel file
    h5write(fullfile(in.folderSingleChannels, allFileNames(i,:)), ['/channel_' num2str(i)], KSRawdata, [1 1], [1 size(KSRawdata,2)]);
end

clear data KSRawdata

%% Files for matrix compilation
out.myFiles       = ls(fullfile(in.folderSingleChannels,'*.h5'));
out.dataDirectory = fullfile(in.folderSingleChannels);

%% If events are used to crop the matrix  
% if in.retrieveEvents == true
%     out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
% else 
%     out.events = false; 
% end

end