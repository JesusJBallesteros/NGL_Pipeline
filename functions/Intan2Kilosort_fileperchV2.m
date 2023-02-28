function Intan2Kilosort_fileperchV2(opt)
% This function is a dependency of the script Intan2Kilosort_wrapperV2,
% only necessary if the recording system in use is Intan.
% Optimized to use only the necessary number of samples (instead of Inf).
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    28th Feb 2023 (Jesus)
 
%% Pre-define .h5 and .bin files
% Create complete HDF5 file matching the size needs.
h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
        '/allChnMat', [opt.numChannels opt.num_samples], ...
        'ChunkSize', [1 opt.HDF5chunkSize], ...
        'Datatype', 'int16')

% Also create a bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 

%% Get data from INTAN data files, as single channels.
% Opens neural data file and coverts the ADC steps to microvolts, then
% writes into the full matrix in a per channel and chunked fashion. 
data_mat = cell(opt.numChannels,1);
    
% Write each channel stepwise into a matrix (hdf5) and into a binary file
% Total number of samples per channel is 'opt.num_samples'
ChunkStart = 1:opt.StpSz:opt.num_samples-mod(opt.num_samples,opt.StpSz);

disp('Writting the .h5 file.');
for i = 1:opt.numChannels

    % each data point of neural data is a 16 bit word. Read as int16 and keep
    % it that way. 
    fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
    data = fread(fid, [1 opt.num_samples], 'int16=>int16'); 
    fclose(fid);
    
    % Data already comes as channels x samples from INTAN. Convert to microvolts.
    % By using int16 instead of double, we round to single digit microvolt
    % values. No practical effect vs the double, where we would keep down to
    % the hundredth of microvolt.
    data_mat{i,1} = data * 0.195;
   
    for j = 1:opt.StpSz:ChunkStart(end)
        if opt.h5
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ... % filename
                   '/allChnMat', ...                        % dataset
                   data_mat{i,1}(1,j:(j+opt.StpSz)-1), ...  % data to be written
                   [i j-(ChunkStart(1)-1)], ...             % chunk start
                   [1 opt.StpSz]);                          % chunk size
        end

    end
end

disp('Done writting the .h5 file. Now writting the .bin file');
fwrite(fidDataMat, cell2mat(data_mat), 'int16');
fclose(fidDataMat);

end