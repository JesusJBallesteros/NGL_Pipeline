function Intan2Kilosort_fileperch(opt)
% This function is a dependency of the script Intan2Kilosort_wrapperV2,
% only necessary if the recording system in use is Intan. Writes the processed 
% data into a .bin file in time-based chunks.
%
% If the full data matrix fits within available RAM, it is processed
% in-memory. Otherwise, a temporary matfile on disk is used as a buffer,
% and channels are processed one at a time to avoid OOM errors. In disk
% mode, CAR cannot be applied (requires all channels in RAM simultaneously)
% and will be skipped with a warning.
%
% Version 27.03.2026 Jesus

%% Check if the full data matrix fits in available RAM
bytesRequired  = opt.numChannels * opt.num_samples * 2; % int16 = 2 bytes

% Windows: MATLAB built-in
[~, sys] = memory();
bytesAvailable = sys.PhysicalMemory.Available;
useRAM         = (bytesRequired * 3) <= bytesAvailable;
fprintf('Overestimated RAM needed: %.2f GB\n', bytesRequired*3/1e9);

if useRAM
    disp('Sufficient RAM available — processing in memory.');
else
    disp('Insufficient RAM — falling back to disk-based (matfile) processing.');
    warning('If CAR was requested it will not be applied in disk mode.');
end

%% Prepare .bin output file
opt.binfilename = fullfile(opt.FolderProcDataMat, opt.SavFileName + ".bin");
if isfile(opt.binfilename)
    delete(opt.binfilename);
end

%%  IN-MEMORY
%  Allocate full matrix, read all channels, apply CAR, filter, write chunks
if useRAM
    data = int16(zeros(opt.numChannels, opt.num_samples));

    disp('Reading data, filtering if necessary.');
    for i = 1:opt.numChannels
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
            tempdata = fread(fid, [1 opt.num_samples], 'int16=>int16');
        fclose(fid);
        data(i,:) = int16(tempdata * 0.195);
    end

    % Common Average Referencing. needs all channels, done before filtering
    if opt.CAR
        disp('Re-referencing by Common Average Referencing (CARing).')
        data = ft_preproc_rereference(data, 'all', 'median');
    end

    % Per-channel filtering (full signal required for twopass and detrend)
    for i = 1:opt.numChannels
        chandata = double(data(i,:));

        fprintf('Detrending Ch %d\n', i)
        chandata = ft_preproc_detrend(chandata);

        if opt.highpass > 0
            fprintf('Highpassing Ch %d at %d Hz\n', i, opt.highpass)
            [chandata, ~, ~] = ft_preproc_highpassfilter(chandata, opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
        end
        if opt.lowpass < 9500
            fprintf('Lowpassing Ch %d at %d Hz\n', i, opt.lowpass)
            [chandata, ~, ~] = ft_preproc_lowpassfilter(chandata, opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
        end

        data(i,:) = int16(chandata);
    end

    % Write to .bin in chunks
    fidDataMat    = fopen(opt.binfilename, 'a');
    chunkSize     = opt.StpSz;
    numFullChunks = floor(opt.num_samples / chunkSize);
    lastChunkSize = mod(opt.num_samples, chunkSize);
    totalChunks   = numFullChunks + (lastChunkSize > 0);

    disp('Writing .bin file in chunks...')
    % fwrite writes column-major: each column (= one time sample across all
    % channels) is written contiguously. Chunks must therefore be shaped as
    % [nChannels x chunkSize] so that concatenated chunks reproduce the 
    % nChannels x nSamples layout on read-back
    for j = 1:numFullChunks
        sampleStart = (j - 1) * chunkSize + 1;
        sampleEnd   =  j      * chunkSize;
        fwrite(fidDataMat, data(:, sampleStart:sampleEnd), 'int16');
        fprintf('  Chunk %d/%d written\n', j, totalChunks);
    end
    if lastChunkSize > 0
        sampleStart = numFullChunks * chunkSize + 1;
        fwrite(fidDataMat, data(:, sampleStart:end), 'int16');
        fprintf('  Chunk %d/%d written (%d samples)\n', totalChunks, totalChunks, lastChunkSize);
    end
    fclose(fidDataMat);

%% IN DISK
%  Use a temporary matfile. Goes one channel at a time so the full matrix
%  never fills the RAM. CAR not possible
else
    tmpMatPath = fullfile(opt.FolderProcDataMat, opt.SavFileName + "_tmp.mat");
    if isfile(tmpMatPath)
        delete(tmpMatPath);
    end
    % Pre-allocate the matfile by writing the last element
    % This extends the file to full size without loading anything into RAM
    mf = matfile(tmpMatPath, 'Writable', true);
    mf.data(opt.numChannels, opt.num_samples) = int16(0);

    disp('Reading and filtering data channel by channel (disk mode).');
    for i = 1:opt.numChannels
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
            tempdata = fread(fid, [1 opt.num_samples], 'int16=>int16');
        fclose(fid);

        chandata = double(int16(tempdata * 0.195));

        fprintf('Detrending Ch %d\n', i)
        chandata = ft_preproc_detrend(chandata);

        if opt.highpass > 0
            fprintf('Highpassing Ch %d at %d Hz\n', i, opt.highpass)
            [chandata, ~, ~] = ft_preproc_highpassfilter(chandata, opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
        end
        if opt.lowpass < 9500
            fprintf('Lowpassing Ch %d at %d Hz\n', i, opt.lowpass)
            [chandata, ~, ~] = ft_preproc_lowpassfilter(chandata, opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
        end

        % Write processed channel directly to disk
        mf.data(i, 1:opt.num_samples) = int16(chandata);
    end

    % Read back from matfile in chunks and write to .bin
    % Each chunk loads only opt.StpSz samples across all channels
    fidDataMat    = fopen(opt.binfilename, 'a');
    chunkSize     = opt.StpSz;
    numFullChunks = floor(opt.num_samples / chunkSize);
    lastChunkSize = mod(opt.num_samples, chunkSize);
    totalChunks   = numFullChunks + (lastChunkSize > 0);

    disp('Writing .bin file in chunks from matfile...')
    % fwrite writes column-major: each column (= one time sample across all
    % channels) is written contiguously. Chunks must therefore be shaped as
    % [nChannels x chunkSize] so that concatenated chunks reproduce the 
    % nChannels x nSamples layout on read-back
    for j = 1:numFullChunks
        sampleStart = (j - 1) * chunkSize + 1;
        sampleEnd   =  j      * chunkSize;
        fwrite(fidDataMat, mf.data(:, sampleStart:sampleEnd), 'int16');
        fprintf('  Chunk %d/%d written\n', j, totalChunks);
    end
    if lastChunkSize > 0
        sampleStart = numFullChunks * chunkSize + 1;
        fwrite(fidDataMat, mf.data(:, sampleStart:opt.num_samples), 'int16');
        fprintf('  Chunk %d/%d written (%d samples)\n', totalChunks, totalChunks);
    end
    fclose(fidDataMat);

    % Clean up temporary matfile
    clear mf;
    delete(tmpMatPath);
end
disp('Done writing .bin file.')

end % main function