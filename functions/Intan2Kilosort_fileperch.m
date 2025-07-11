function Intan2Kilosort_fileperch(opt)
% This function is a dependency of the script Intan2Kilosort_wrapperV2,
% only necessary if the recording system in use is Intan.
% Optimized to use only the necessary number of samples (instead of Inf).
%
% Dependencies: bandfilter
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
%
% Version 24.08.2024 Jesus
 
%% Pre-define .h5 and .bin files
% Create a complete HDF5 file matching the size needs.
% opt.filename = fullfile(opt.FolderProcDataMat, [opt.SavFileName + ".h5"]); 
% if isfile(opt.filename)
%     delete(opt.filename);
% end
% 
% opt.dataset = '/allChnMat'; % for now, as before.
% 
% if opt.HDF5chunkSize > opt.num_samples
%     opt.HDF5chunkSize = opt.num_samples;
% end
% 
% h5create(opt.filename,                       ... % filename
%          opt.dataset,                        ... % dataset name
%          [opt.numChannels opt.num_samples],  ... % prepare data dimensions (nCh x samples).
%          'ChunkSize', [1 opt.HDF5chunkSize], ... % prepare to write in chunks in time dimension
%          'Datatype', 'int16');                   % Set data precision

%% Get data from INTAN data files, as single channels.
% Opens neural data file and coverts the ADC steps to microvolts, then
% write channel-by-channel into .h5 file and the whole array into a .bin
% file. Filters the data if necessary/requested.  
data = int16(zeros(opt.numChannels,opt.num_samples));

% Create or open a bin file. Append data at end.
opt.binfilename = fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]);
if isfile(opt.binfilename)
    delete(opt.binfilename);
end

fidDataMat = fopen(opt.binfilename, 'a');

disp('Reading data, filtering if necessary.');
for i = 1:opt.numChannels
    % Each sample of neural data is a 16 bit word. Read as int16,
    % and keep it that way. 
    fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
        tempdata = fread(fid, [1 opt.num_samples], 'int16=>int16');
    fclose(fid);
  
    % Data comes as channels x samples from INTAN. Convert to microvolts.
    % By using int16, we round to single digit microvolt values. No practical
    % effect vs the double, where we would keep down to the hundredth of microvolt.
    tempdata = tempdata * 0.195;
   
    % Back to 'int16'.
    data(i,:) = int16(tempdata);
    clear tmp 
       
%     % write each channel as a whole into the matrix (hdf5)
%     h5write(opt.filename, ... % filename
%             opt.dataset,  ... % dataset name
%             data,         ... % data to write as ch x smp
%             [i 1],        ... % channels to write
%             [1 opt.num_samples]); % samples to write
end

% Proceed with filters. Set variables in output if you want to have the
% exact values applied during filtering.

if opt.CAR
    % In principle, data from a single HS on a single region.
    disp('Re-referencing by Common Average Referencing (CARing).')
    data = ft_preproc_rereference(data, 'all', 'median');
end

% 'bandfilter' has a whole description inside the fuction. Check
% for optional arguments and mechanism of work. It likes data in
% 'double' precision, so we convert it within the line. Afterwards
% needs to be reverted to 'int16'.
%    [tmp, ~, ~] = bandFilter(double(tempdata), [], opt.highpass, opt.sampleRate);

txt = sprintf('Highpass filter set at %d Hz. It may take a moment.\n', opt.highpass);
fprintf(txt);
% Keep memory usage low doing one channel at a time.
for i = 1:opt.numChannels
    fprintf('- Filtering channel %d of %d.\n', i, opt.numChannels);
    
    % Detrend channel (remove DC)
    disp('Detrending...')
    data(i,:) = ft_preproc_detrend(data(i,:));

    % Highpass channel (Butterwort, 6th order, back&forth)
    disp('Filtering...')
    [data(i,:), ~, ~] = ft_preproc_highpassfilter(data(i,:), opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
end

%% Write bin file.
% disp('Writting .bin file. Depending on the data size, this may take a while...')
% iter = 0;
% lastchunk = mod(opt.num_samples, opt.StpSz)+1;
% chunks = 1:opt.StpSz:opt.num_samples-lastchunk;
% 
% for j = chunks
%     iter = iter +1;
%     fprintf('Chunk %d/%d. \n', iter, length(chunks))
%     Chunk = int16(zeros(opt.numChannels, opt.StpSz));
% 
%     for k = 1:opt.numChannels
%         if j < opt.num_samples-mod(opt.num_samples, opt.StpSz)+1
%             Chunk(k,:) = h5read(opt.filename, ... % filename
%                                 opt.dataset,  ... % dataset name
%                                 [k j],        ... % chunk and channel to write       
%                                 [1 opt.StpSz]);   % samples to write
%         else
%             % Last chunk, normally smaller than stepSize
%             Chunk(k,:) = h5read(opt.filename, ... % filename
%                                 opt.dataset,  ... % dataset name
%                                 [k j],        ... % chunk and channel to write    
%                                 [1 lastchunk]);   % samples to write
%         end
%     end
%     
%     fwrite(fidDataMat, Chunk, 'int16');
% end
% 
% % Close the .bin file
% fclose(fidDataMat);
% 
% delete(opt.filename);

fwrite(fidDataMat, data, 'int16');
fclose(fidDataMat);

end