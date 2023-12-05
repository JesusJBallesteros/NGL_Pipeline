function out = Intan2Kilosort_fileperch(opt)
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
% To obtain the number of samples per file, first read file info. In this
% case, it is a file per channel but all should be equal. Read the first
% channel.
fileinfo = dir(opt.myFiles(1).name);

% Divide the file size by bytes each int16 word takes
opt.num_samples = fileinfo.bytes/2; % int16 = 2 bytes

% Create an empty file per channel:
% Pre-allocated memory accounts for exact 'num_samples' of neural data (vs.
% HDF5 files with infinite slots).
for i = 1:opt.numChannels
    IndChnl  = ['Channel_',sprintf('%03d',i)];
    fileName = fullfile(opt.FolderSingleChannels, IndChnl);

    if ~isfile([fileName,'.h5'])
        h5create([fileName '.h5'], ...
                ['/channel_' num2str(i)],...
                [1 opt.num_samples],...
                'ChunkSize', [1 opt.HDF5chunkSize], ...
                'Datatype','int16');
    end
end

% list with all pre-allocated files
allFileNames = ls([opt.FolderSingleChannels,'\Ch*']); 

%% Allocate data to its respective single-channel file
% Opens neural data file and coverts the ADC steps to microvolts, then
% writes one file per channel in .h5 format. 
for i = 1:length(opt.myFiles)
    fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
        % each data point of neural data is a 16 bit word. Read as int16 and keep
        % it that way. 
        data = fread(fid, [1 opt.num_samples], 'int16=>int16'); 
    fclose(fid);
    
    % Data already comes as channels x samples from INTAN. Convert to microvolts.
    % By using int16 instead of double, we round to single digit microvolt
    % values. No practical effect vs the double, where we would keep down to
    % the hundredth of microvolt.
    data = data * 0.195;
    
    % distribute each row of data to its respective single-channel file.
    h5write(fullfile(opt.FolderSingleChannels, ...
        allFileNames(i,:)), ...
        ['/channel_' num2str(i)], ...
        data, [1 1], ...
        [1 opt.num_samples]);
end
clear data

%% Files for matrix compilation
out.myFiles         = ls(fullfile(opt.FolderSingleChannels,'*.h5'));
out.dataDirectory   = fullfile(opt.FolderSingleChannels);
out.num_samples     = opt.num_samples;

end