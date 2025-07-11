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

% Jesus 12.06.2024

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

%% Common methods of preprocessing. Re-Referencing, DC substraction and filter.
if opt.CAR
    % In principle, data from a single HS on a single region.
    disp('Re-referencing by Common Average Referencing (CARing).')
    data_mat = ft_preproc_rereference(data_mat, 'all', 'median');
end

% Enforce int16
filt_data_mat = int16([]);

txt = sprintf('Highpass filter set at %d Hz. It may take a moment.\n', opt.highpass);
fprintf(txt);

% Keep memory usage low doing one channel at a time.
for b = 1:opt.numChannels
    fprintf('- Filtering channel %d of %d.\n', b, opt.numChannels);
    
    % Detrend channel (remove DC)
    disp('Detrending...')
    filt_data_mat(b,:) = ft_preproc_detrend(data_mat(b,:));

    % Highpass channel (Butterwort, 6th order, back&forth)
    disp('Filtering...')
    [filt_data_mat(b,:), ~, ~] = ft_preproc_highpassfilter(data_mat(b,:), opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
                
end

%% Write bin file.
fwrite(fidDataMat, filt_data_mat, 'int16');
fclose(fidDataMat);

end  