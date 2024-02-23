function Deuteron2Kilosort(opt)
% This function is a dependency of the script Deuteron2Kilosort_wrapper.
% It compiles the data saved DF/DT2 in mat file channel by channel.
% It performs high-pass filtering when the data is not already like that.
%
% DEPENDENCIES:
%  Deuteron_extractData: For new .DF1 format, data is extracted by this
%                 function.
%  bandFilter: to filter the highpass signal, since signal is intended to
%                 Kilosort.

% Jesus 22.02.2024

% if ~isfield(opt,'h5'),             opt.h5                 = false;         end

%% Pre-define .h5 (Deprecated)
% filename = fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']); 
% dataset = '/allChnMat'; % Name dataset to an useful denomination?
% if opt.h5
%     % Create complete HDF5 file matching size needs.
%     h5create(filename,                          ... % filename.
%              dataset,                           ... % dataset name.
%              [opt.numChannels Inf],              ... % prepare data dimensions (nCh x samples).
%              'ChunkSize', [1 opt.HDF5chunkSize], ... % prepare to write chunks in time dimension.
%              'Datatype', 'int16');                   % data precision.
% end

%% Check bin file existence and completion.
% If an error happens during processing, the bin file persists created but with zero size
if isfile(fullfile(opt.FolderProcDataMat, [opt.SavFileName + ".bin"]))
    bininfo = dir(fullfile(opt.FolderProcDataMat, [opt.SavFileName + ".bin"]));
    if bininfo.bytes > 0
        disp('An actual .bin file exists in this directory. Conversion Skipped.')
        return
    end
    disp('A failed .bin file was found in this directory. Will be overwritten.')
end

%% Create an new .bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]), 'a'); 

%% DT2 format conversion.
if strcmp(opt.ext, 'DT2')
    disp('Format is DT2. Deprecating.')

    % Prepare a cell array with nChannels
    data_mat = cell(opt.numChannels,1);

    % Initialize variable to chunk the writting. 
    indexPos  = 0;

    % Go over individual files
    for i = 1:length(opt.myFiles)
        
        % Neural data points are unsigned 16 bit words. Read.
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
            tempdata = fread(fid, 'uint16');
        fclose(fid);

        % Remove trailing zeroes in last file
        if i == length(opt.myFiles)
            tempdata(tempdata==0) = [];
        end

        % Because all channels come concatenated, we need to reshape as
        % channels x samples, knowing the nChannels.
        tempdata = reshape(tempdata', opt.numChannels, []); 

        % Convert ADC steps into microvolts, so conversion to int16 is
        % possible without loss.
        tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);
        
        % Get nSamples coming from this file. Should stay constant, until
        % last file which normally will be smaller.
        nSamples = size(tempdata,2);

        % Distribute each channel to its respective slot in new array
        for b = 1:opt.numChannels
            % if opt.h5
            %     h5write(filename,       ... % filename.
            %             dataset,        ... % dataset name.
            %             tempdata(b,:),  ... % data of a channel stored in the DT2 file (already scaled)
            %             [b indexPos+1], ... % Write channel b, from starting sample
            %             [1 nSamples]);      %  and this amount of samples.
            % end
            
            % Write channel into general cell array.
            data_mat{b,1} = [data_mat{b,1} tempdata(b,:)];
        end

        % Index for next chunk's starting sample
        indexPos = indexPos + nSamples;

    end
    clear tempdata fid

    % Convert cell to array, and to single
    data_mat = cell2mat(data_mat);

end 
%% DF1 format conversion.
if strcmp(opt.ext, 'DF1')
    data_mat        = [];   % Create empty variable to store all data (do not pre-allocate the whole matrix)
    opt.stream      = 1;    % Pass variable to read continuous neural signals.
    
    % Open each channel file and read it, resize data to fit Kilosort expectations,
    % and concatenate it to the matrix consecutively.
    disp('Format is DF1. Proceeding. It may take a moment.')
    for i = 1:length(opt.myFiles)
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r'); % Open file.
            tempdata = Deuteron_extractData(fid, opt); % read data with Deuteron's script.
        fclose(fid); % Close file

        % Convert ADC bit steps into microvolts to convert to int16 without loss.
        tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);

        try
            data_mat = [data_mat tempdata]; % Concatenate.
        catch
            disp('found non-matching file') % If size does not match (tipically last chunck).
        end
    end
    clear tempdata fid

    % Reshape to channels x samples.
    data_mat = reshape(data_mat, opt.numChannels, []);

end

%% Common methods of preprocessing. DC substraction, Referencing and filter.
if opt.CAR
    % Subtract the mean from each channel
    data_mat = single(data_mat);
    data_mat = data_mat - mean(data_mat, 1);

    % CAR, common average referencing by median.
    data_mat = data_mat - median(data_mat, 2); % subtract median across channels
end

% Highpass filter.
filt_data_mat = int16([]);
txt = sprintf('Filtering between %d and %d Hz. It may take a moment.\n', opt.highpass(1), opt.highpass(2));
fprintf(txt);

% Keep memory usage low doing one channel at a time.
for b = 1:opt.numChannels
    [filt_data_mat(b,:), ~, ~] = bandFilter(double(data_mat(b,:)), [], opt.highpass, opt.sampleRate);
    filt_data_mat(b,:) = int16(filt_data_mat(b,:));
end

% %% HDF5 method. Distribute each row of data to its respective single-channel file
% if opt.h5
%     for b = 1:opt.numChannels
%         h5write(filename,         ... % filename.
%                 dataset,          ... % dataset name.
%                 data_mat(b,:),        ... % channel b, complete
%                 [b b],                ... %   into slot b
%                 [1 size(data_mat,2)]);    %   as long as it is.
%     end
% end

%% Write bin file.
fwrite(fidDataMat, filt_data_mat, 'int16');
fclose(fidDataMat);

end  