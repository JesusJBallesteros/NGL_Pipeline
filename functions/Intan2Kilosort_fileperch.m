function Intan2Kilosort_fileperch(opt)
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
% Version 13.04.2023 Jesus
 
% TODO 
% Downsampling highpass data to perhaps 15kHz?

%% Pre-define .h5 and .bin files
%  % In case we want individual .h5 files back. But should not be necessary.
%  % Create individual .h5 files for memory use reduction
% for i = 1:opt.numChannels
%     IndChnl  = ['Channel_', sprintf( ['%0',num2str(numel(num2str(opt.numChannels))),'d'], i )];
%     fileName = fullfile(opt.FolderProcDataMat,IndChnl);
%     h5create([fileName '.h5'], ...
%             ['/channel_' num2str(i)], ...
%             [1 Inf], ...
%             'ChunkSize', [1 opt.HDF5chunkSize], ...
%             'Datatype','int16')
% end
% allh5files = ls([opt.FolderProcDataMat,'\Ch*']); % list with all Ch files

% Create complete HDF5 file matching the size needs.
% Name dataset to an useful denomination?
if opt.h5
    opt.filename = fullfile(opt.FolderProcDataMat, [opt.SavFileName + ".h5"]); 
    opt.dataset = '/allChnMat'; % for now, as before.
    
    h5create(opt.filename,                       ... % filename
             opt.dataset,                        ... % dataset name
             [opt.numChannels opt.num_samples],  ... % prepare data dimensions (nCh x samples).
             'ChunkSize', [1 opt.HDF5chunkSize], ... % prepare to write in chunks in time dimension
             'Datatype', 'int16');                   % Set data precision
end


%% Get data from INTAN data files, as single channels.
% Opens neural data file and coverts the ADC steps to microvolts, then
% write channel-by-channel into .h5 file and the whole array into a .bin
% file. Filters the data if necessary/requested.

% Empty cell array
data_mat = cell(1,1);
    
% Determine the starting sample number of each chunk to write in the .h5 file
% ChunkStart = 1:opt.StpSz:opt.num_samples-mod(opt.num_samples,opt.StpSz);

disp('Reading data, filtering if necessary and writting the .h5 file if requested.');
% indexPos = 0;
for i = 1:opt.numChannels

% Create or open a bin file. Append data at end.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]), 'a+'); 

    % Each sample of neural data is a 16 bit word. Read as int16,
    % and keep it that way. 
    fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
        data = fread(fid, [1 opt.num_samples], 'int16=>int16'); 
        % By using int16, we round to single digit microvolt values. No practical
        % effect vs the double, where we would keep down to the hundredth of microvolt.
    fclose(fid);
    
    % Data comes as channels x samples from INTAN. Convert to microvolts.
    data_mat{1,1} = data * 0.195;
    clear data

    % Proceed with filter. Set variables in output if you want to have the
    % exact values applied during filtering.
    if opt.set_filter
        % 'bandfilter' has a whole description inside the fuction. Check
        % for optional arguments and mechanism of work. It like data in
        % 'double' precision, so we convert it within the line. Afterwards
        % needs to be reverted to 'int16'.
        [tmp, ~, ~] = bandFilter(double(data_mat{1,1}), [], opt.highpass, opt.sampleRate);

        % Back to 'int16'.
        data_mat{1,1} = int16(tmp);

        % No downsampling needed for highpass bands. 
        % (However, since we highpass only up to 7500 Hz, we could, in theory,
        % reduce the data to a half by downsampling to 15000 Hz.
    end
   
%      % In case we want individual .h5 files back. But should not be necessary.
%      % Distribute each row of data to its respective single-channel file
%     h5write(fullfile(opt.FolderProcDataMat, allh5files(i,:)),...
%             ['/channel_' num2str(i)], ...	% dataset name
%             data_mat{i,1}, ...              % dataset: data of a channel stored in the DT2 file (already scaled)
%             [1 indexPos+1], ...              
%             [1 size(data_mat{i,1}, 2)]);
% 
%     indexPos = indexPos+size(data_mat{i,1},2);

    fprintf('Writting the .bin file: Channel %d.\n', i);
    % Write the current channel to the .bin file.
    fwrite(fidDataMat, data_mat{1,1}, 'int16');
    fclose(fidDataMat);
end

% Close the .bin file

end