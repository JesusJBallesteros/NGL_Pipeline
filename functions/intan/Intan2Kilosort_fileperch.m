function Intan2Kilosort_fileperch(opt)
% Intan2Kilosort_fileperch  Convert INTAN fileperch amp*.dat to Kilosort .bin file.
%
% PURPOSE:
%   Reads one-file-per-channel INTAN recordings (amp-A-000.dat, amp-A-001.dat, …),
%   scales to µV (×0.195), optionally applies detrend + high-pass + low-pass
%   filtering, and writes a flat int16 interleaved binary file for Kilosort. 
%   Uses an adaptive RAM-vs-disk strategy to handle large datasets.
%
% USAGE:
%   Intan2Kilosort_fileperch(opt)
%   Called from Intan2Kilosort_wrapper; do not call directly.
%
% INPUT:
%   opt  - options struct set by Intan2Kilosort_wrapper; must contain:
%            .numChannels    number of channels to process
%            .num_samples    samples per channel (from file size)
%            .PathRaw        raw data folder
%            .myFiles        dir-struct list of amp*.dat files
%            .FolderProcDataMat  output folder
%            .SavFileName    session name (output file = <name>.bin)
%            .set_filter     1 = apply filtering pipeline; 0 = raw pass-through
%            .sampleRate     raw sample rate (Hz)
%            .lowpass        spike-band low-pass cutoff (Hz), [] = off
%            .highpass       high-pass cutoff (Hz), 0 = off
%            .CAR            common-average referencing flag (0 = off)
%            .StpSz          chunk size in samples for .bin write loop
%
% OUTPUT:
%   <opt.SavFileName>.bin written to opt.FolderProcDataMat
%   Format: int16, channels interleaved [nChannels × nSamples], no header.
%
% MEMORY STRATEGY:
%   If 3× the data matrix fits in available physical RAM, processes in memory
%   (faster, supports CAR). Otherwise falls back to disk (matfile) mode,
%   processing one channel at a time (CAR not possible in disk mode).
%
% Version 07.05.2026 (Jesus)

%% Check if the full data matrix fits in available RAM
bytesRequired  = opt.numChannels * opt.num_samples * 2; % int16 = 2 bytes

% check memory in Windows. MATLAB built-in
[~, sys] = memory();
bytesAvailable = sys.PhysicalMemory.Available;
fprintf('Overestimated RAM needed: %.2f GB\n', bytesRequired*3/1e9);
useRAM = (bytesRequired * 3) <= bytesAvailable;

%% Prepare .bin output file
opt.binfilename = fullfile(opt.FolderProcDataMat, opt.SavFileName + ".bin");
if isfile(opt.binfilename)
    delete(opt.binfilename);
end

%% IN-MEMORY
% Allocate full matrix, read all channels, apply CAR, filter, write chunks
if useRAM
    disp('Sufficient RAM available — processing in memory.');
    data = int16(zeros(opt.numChannels, opt.num_samples));

    disp('Reading data, filtering if necessary.\n');
    for i = 1:opt.numChannels
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
            tempdata = fread(fid, [1 opt.num_samples], 'int16=>int16');
        fclose(fid);
        data(i,:) = int16(tempdata * 0.195);
    end

    % Common Average Referencing. needs all channels, done before filtering
    if opt.CAR
        disp('Re-referencing by Common Average Referencing (CAR).\n')
        data = ft_preproc_rereference(data, 'all', 'median');
    end

    % Per-channel filtering (full signal required for twopass and detrend)
    for i = 1:opt.numChannels
        chandata = double(data(i,:));

        fprintf('Detrending Ch %d\n', i)
        chandata = ft_preproc_detrend(chandata);

        if opt.set_filter == 1
            if opt.highpass > 0
                fprintf('Highpassing Ch %d at %d Hz\n', i, opt.highpass)
                [chandata, ~, ~] = ft_preproc_highpassfilter(chandata, opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
            end
            if opt.lowpass < 9500 && opt.lowpass > 0
                fprintf('Lowpassing Ch %d at %d Hz\n', i, opt.lowpass)
                [chandata, ~, ~] = ft_preproc_lowpassfilter(chandata, opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
            end
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
        fprintf(' Chunk %d/%d written\n', j, totalChunks);
    end

    if lastChunkSize > 0
        sampleStart = numFullChunks * chunkSize + 1;
        fwrite(fidDataMat, data(:, sampleStart:end), 'int16');
        fprintf(' Chunk %d/%d written (%d samples)\n', totalChunks, totalChunks, lastChunkSize);
    end
    fclose(fidDataMat);

%% IN DISK
%  Use a temporary matfile. Goes one channel at a time so the full matrix
%  never fills the RAM. CAR not possible
else
    disp('Insufficient RAM, using in-disk (matfile) processing.');
    warning('If CAR was requested, it will not be applied.');

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

        %chandata = double(int16(tempdata * 0.195)); % reason: (tempdata*0.195) as double is immediately truncated to int16 and widened to double again.truncation discards sub-integer precision that would otherwise survive to the filter
        chandata = double(tempdata) * 0.195; % modified 05.05.2026

        fprintf('Detrending Ch %d\n', i)
        chandata = ft_preproc_detrend(chandata);

        if opt.set_filter == 1
            if opt.highpass > 0
                fprintf('Highpassing Ch %d at %d Hz\n', i, opt.highpass)
                [chandata, ~, ~] = ft_preproc_highpassfilter(chandata, opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
            end
            if opt.lowpass < 9500 &&  opt.lowpass > 0
                fprintf('Lowpassing Ch %d at %d Hz\n', i, opt.lowpass)
                [chandata, ~, ~] = ft_preproc_lowpassfilter(chandata, opt.sampleRate, opt.lowpass, 6, 'but', 'twopass');
            end
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