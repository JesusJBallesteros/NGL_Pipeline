function out = Intan2Kilosort_filepertype(opt)
% This function is a dependency of the script Intan2Kilosort_wrapper,
% only necessary if the recording system in use is Intan.
% Optimized to use only the necessary number of samples (instead of Inf).
% It compiles the data save as filepertype format in a HDF5file per channel
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    28th Feb 2023 (Jesus)
 
%% Total session Data, divided per channels

% Create an empty file per channel:
% Pre-allocated memory accounts for exact 'num_samples' of neural data (vs.
% HDF5 files with infinite slots).
for i = 1:opt.numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(opt.FolderSingleChannels, IndChnl);

    if ~isfile([fileName,'.h5'])
        h5create([fileName '.h5'], ...
                ['/channel_' num2str(i)],...
                [1 opt.num_samples], ...
                'ChunkSize', [1 opt.HDF5chunkSize], ...
                'Datatype','int16');
    end
end

% List with all pre-allocated files
allFileNames = ls([opt.FolderSingleChannels,'\Ch*']); 

%% Allocate data to its respective single-channel file
% Opens neural data file and coverts the ADC steps to microvolts, then
% writes one file per channel in .h5 format. 
fid = fopen(fullfile(opt.PathRaw, opt.myFiles.name));
    % each data point of neural data is a 16 bit word
    data = fread(fid, [opt.numChannels opt.num_samples], 'int16=>int16'); 
fclose(fid);

% Data already comes as channels x samples from INTAN. Convert to microvolts.
% By using int16 instead of double, we round to single digit microvolt
% values. No practical effect vs the double, where we would keep down to
% the hundredth of microvolt.
data = data * 0.195; 

% Distribute each row of data to its respective single-channel file
for b = 1:opt.numChannels
    h5write(fullfile(opt.FolderSingleChannels, allFileNames(b,:)),...
            ['/channel_' num2str(b)],...
            data(b,:), [1 1], ...
            [1 opt.num_samples]);
end
clear data

%% Files for matrix compilation
out.myFiles       = ls(fullfile(opt.FolderSingleChannels,'*.h5'));
out.dataDirectory = fullfile(opt.FolderSingleChannels);

end