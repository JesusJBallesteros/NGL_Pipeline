function Intan2Kilosort_filepertype(opt)
% Dependency of Intan2Kilosort_wrapper. Only to use in INTAN.
% It compiles the data, saves as HDF5 and .bin files.
% Only for INTAn data saved as file per type. TO DEPRECATE.
%
% Last 07.05.2026 Jesus
 
%% Pre-define .h5 and .bin files
% Create complete HDF5 file matching the size needs.
h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ... % path
        '/allChnMat', [opt.numChannels opt.num_samples], ... % name, size
        'ChunkSize', [1 opt.HDF5chunkSize], ... % sample chunk
        'Datatype', 'int16') % data type

% Also create a bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 

%% Open neural data file
% By using int16 instead of double, we round to single digit microvolt
% values. No practical effect vs the double, where we would keep down to
% the hundredth of microvolt.
fid = fopen(fullfile(opt.PathRaw, opt.myFiles.name)); % open
    % Data already comes as channels x samples from INTAN. 
    data = fread(fid, [opt.numChannels opt.num_samples], 'int16=>int16'); % read data as 16bit integer
fclose(fid);

% Convert to microvolts.
data = data * 0.195; 

% Write channel-wise into hdf5 and binary files
% Set the sample of each data chunk. Only full chunks letting a non-full
% chunk for the end [mod(...)]
rest = mod(opt.num_samples,opt.StpSz);
ChunkStart = [1 opt.num_samples-rest];

disp('Writting the .h5 file.');
for i = 1:opt.numChannels % per channel
    for j = 1:opt.StpSz:ChunkStart(end) % per chunk
        % compile channels -> unaltered matrix to load into kilosort
        h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ... % filename
               '/allChnMat', ...                        % dataset
               data(i,j:(j+opt.StpSz)-1), ...           % data to be written
               [i j], [1 opt.StpSz]);                   % chunk start, chunk size
    end

    % last Chunk per channel
    h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ... % filename
           '/allChnMat', ...                                    % dataset
           data(i, ChunkStart(end):ChunkStart(end)+rest), ...   % data to be written
           [i ChunkStart(end)], [1 rest]);                      % chunk start, chunk size
end
disp('Done writting the .h5 file. Now writting the .bin file');

% Write the bin file, all at once.
fwrite(fidDataMat, data, 'int16');
fclose(fidDataMat);
disp('Done writting the .bin file.');
end