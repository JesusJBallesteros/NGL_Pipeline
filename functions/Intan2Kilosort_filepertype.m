function out = Intan2Kilosort_filepertype(opt)
% This function is a dependency of the script Intan2Kilosort_wrapper,
% only necessary if the recording system in use is Intan
% It compiles the data save as filepertype format in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    Feb 2023 (Jesus)
 
%% Total session Data, divided per channels
% Create an empty file per channel:
% space is pre-allocated to save every sample of neural data (HDF5 files with infinite slots)
for i = 1:opt.numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(opt.FolderSingleChannels,IndChnl);

    if ~isfile([fileName,'.h5'])
        h5create([fileName '.h5'],['/channel_' num2str(i)],[1 Inf],'ChunkSize',[1 opt.HDF5chunkSize],'Datatype','int16')
    end
end

% list with all pre-allocated files
allFileNames = ls([opt.FolderSingleChannels,'\Ch*']); 

%% Allocate data to its respective single-channel file
% Opens neural data file, resizes data for detection with Kilosort,
% allocates data to its respective single-channel file.
fid = fopen(fullfile(opt.PathRaw, opt.myFiles.name));
data = fread(fid, [opt.numChannels inf], 'int16'); % each data point of neural data is a 16 bit word
fclose(fid);

% Data already comes as channels x samples from INTAN
% resize for detection with Kilosort
KSRawdata = int16(int32(data) - int32(intmax('int16')/2)); 

% distribute each row of data to its respective single-channel file
for b = 1:opt.numChannels
    h5write(fullfile(opt.FolderSingleChannels,allFileNames(b,:)), ['/channel_' num2str(b)], KSRawdata(b,:), [1 1], [1 size(KSRawdata,2)]);
end
clear KSRawdata data

%% Files for matrix compilation
out.myFiles       = ls(fullfile(opt.FolderSingleChannels,'*.h5'));
out.dataDirectory = fullfile(opt.FolderSingleChannels);

%% If events are used to crop the matrix  
% if in.retrieveEvents == true
%     out.events = Deuteron_importEventsWithDLL(in.dllFolder,in.pathRaw, sampleRate);
% else 
%     out.events = false; 
% end

end