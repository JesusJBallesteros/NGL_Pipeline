function Intan2Kilosort_filepertype(opt)
% This function is a dependency of the script Intan2Kilosort_wrapper,
% only necessary if the recording system in use is Intan.
% Optimized to use only the necessary number of samples (instead of Inf).
% It compiles the data save as filepertype format in a HDF5file per channel
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    28th Feb 2023 (Jesus)
%
% Version 01.03.2023 Jesus
 
%% Pre-define .h5 and .bin files
    % Create complete HDF5 file matching the size needs.
    h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
            '/allChnMat', [opt.numChannels opt.num_samples], ...
            'ChunkSize', [1 opt.HDF5chunkSize], ...
            'Datatype', 'int16')

% Also create a bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 

%% Open neural data file and coverts the ADC steps to microvolts, then
% writes one file per channel in .h5 format.
% Each data point of neural data is a 16 bit word
fid = fopen(fullfile(opt.PathRaw, opt.myFiles.name));
    data = fread(fid, [opt.numChannels opt.num_samples], 'int16=>int16'); 
fclose(fid);

% Data already comes as channels x samples from INTAN. Convert to microvolts.
% By using int16 instead of double, we round to single digit microvolt
% values. No practical effect vs the double, where we would keep down to
% the hundredth of microvolt.
data = data * 0.195; 

% Write each channel stepwise into a matrix (hdf5) and into a binary file
% Total number of samples per channel is 'opt.num_samples'.
ChunkStart = 1:opt.StpSz:opt.num_samples-mod(opt.num_samples,opt.StpSz);

disp('Writting the .h5 file.');
    for i = 1:opt.numChannels
        for j = 1:opt.StpSz:ChunkStart(end)
            % compile channels -> unaltered matrix to load into kilosort
            h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ... % filename
                   '/allChnMat', ...                        % dataset
                   data(i,j:(j+opt.StpSz)-1), ...           % data to be written
                   [i j-(ChunkStart(1)-1)], ...             % chunk start
                   [1 opt.StpSz]);                          % chunk size
        end
    end

disp('Done writting the .h5 file. Now writting the .bin file');
fwrite(fidDataMat, data, 'int16');
fclose(fidDataMat);

end