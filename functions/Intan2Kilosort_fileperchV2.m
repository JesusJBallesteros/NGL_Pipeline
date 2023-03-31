function Intan2Kilosort_fileperchV2(opt)
% This function is a dependency of the script Intan2Kilosort_wrapperV2,
% only necessary if the recording system in use is Intan.
% Optimized to use only the necessary number of samples (instead of Inf).
%
% Dependencies: bandfilter
%               variety of read/write commands
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
%
% Version 07.03.2023 Jesus
 
% TODO 
% Downsampling highpass data to perhaps 15kHz?
% Name dataset to an useful denomination?
opt.filename = fullfile(opt.FolderProcDataMat, [opt.SavFileName + '.h5']); 
opt.dataset = '/allChnMat'; % for now, as before.

%% Pre-define .h5 and .bin files
% Create complete HDF5 file matching the size needs.
if opt.h5
    h5create(opt.filename,                       ... % filename
             opt.dataset,                        ... % dataset name
             [opt.numChannels opt.num_samples],  ... % prepare data dimensions (nCh x samples).
             'ChunkSize', [1 opt.HDF5chunkSize], ... % prepare to write in chunks in time dimension
             'Datatype', 'int16');                   % Set data precision
end

% Also create a bin file. 
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName + '.bin']), 'a'); 

%% Get data from INTAN data files, as single channels.
% Opens neural data file and coverts the ADC steps to microvolts, then
% write channel-by-channel into .h5 file and the whole array into a .bin
% file. Filters the data if necessary/requested.

% Empty cell array with numChannels
data_mat = cell(opt.numChannels,1);
    
% Determine the starting sample number of each chunk to write in the .h5 file
ChunkStart = 1:opt.StpSz:opt.num_samples-mod(opt.num_samples,opt.StpSz);

disp('Reading data, filtering if necessary and writting the .h5 file if requested.');
for i = 1:opt.numChannels

    % Each sample of neural data is a 16 bit word. Read as int16,
    % and keep it that way. 
    fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
        data = fread(fid, [1 opt.num_samples], 'int16=>int16'); 
        % By using int16, we round to single digit microvolt values. No practical
        % effect vs the double, where we would keep down to the hundredth of microvolt.
    fclose(fid);
    
    % Data comes as channels x samples from INTAN. Convert to microvolts.
    data_mat{i,1} = data * 0.195;

    % Proceed with filter. Set variables in output if you want to have the
    % exact values applied during filtering.
    if opt.set_filter
        % 'bandfilter' has a whole description inside the fuction. Check
        % for optional arguments and mechanism of work. It like data in
        % 'double' precision, so we convert it within the line. Afterwards
        % needs to be reverted to 'int16'.
        [tmp, ~, ~] = bandFilter(double(data_mat{i,1}), [], opt.highpass, opt.sampleRate);

        % Back to 'int16'.
        data_mat{i,1} = int16(tmp);

        % No downsampling needed for highpass bands. 
        % (However, since we highpass only up to 7500 Hz, we could, in theory,
        % reduce the data to a half by downsampling to 15000 Hz.
    end
   
    if opt.h5
        % We have a full channel ready, write it to the .h5 file in chunks.
        for j = 1:opt.StpSz:ChunkStart(end)
            h5write(opt.filename,                        ... % filename
                    opt.dataset,                         ... % dataset.
                    data_mat{i,1}(1, j:(j+opt.StpSz)-1), ... % data to be written (ch x samples)
                    [i j-(ChunkStart(1)-1)],             ... % for channel i, starting sample
                    [1 opt.StpSz]);                          % amount of samples to write
        end
    end
end
clear tmp data fid

disp('Writting the .bin file.');
% We have a full array of channels ready, write it to the .bin file at once.
fwrite(fidDataMat, cell2mat(data_mat), 'int16');
fclose(fidDataMat);

end