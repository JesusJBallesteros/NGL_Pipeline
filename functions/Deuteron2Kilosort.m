function Deuteron2Kilosort(opt)
% This function is a dependency of the script Deuteron2Kilosort_wrapper,
% only necessary if the recording system in use is Deuteron
% It compiles the data save in DT2 opt.myFiles in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
% Filtering and the elimination of experimentally irrelevant recording 
% periods only seem to be necessary when recording with Deuteron.
%
% DEPENDENCIES:
%  Deuteron_extractData: For new .DF1 format, data is extracted by this
%                 function.
%  bandFilter: to filter the highpass signal, since signal is intended to
%                 Kilosort.
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    21.12.2023 (Jesus)

% % Pre-define .h5
% % Name dataset to an useful denomination?
% filename = fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']); 
% dataset = '/allChnMat'; % for now, as before.
% opt.h5 = true;
% if opt.h5
%     % Create complete HDF5 file matching size needs.
%     h5create(filename,                          ... % filename.
%              dataset,                           ... % dataset name.
%              [opt.numChannels Inf],              ... % prepare data dimensions (nCh x samples).
%              'ChunkSize', [1 opt.HDF5chunkSize], ... % prepare to write chunks in time dimension.
%              'Datatype', 'int16');                   % data precision.
% end

%% Check bin file existence and if it was a failure or not.
if isfile(fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]))
    bininfo = dir(fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]));
    if bininfo.bytes > 0
        disp('An actual .bin file exists in this directory. Skipping.')
        return
    end
    disp('A failed .bin file was found in this directory. Will be overwritten.')
end

% Create an empty .bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]), 'a'); 

%% DT2 FORMAT
if strcmp(opt.ext, 'DT2')
    disp('Format is FLAT. Deprecating.')

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

    % Because this data goes to Kilosort, apply highpass filter to data.
    if opt.set_filter
        % Proceed in a channel by channels basis (lower memory use)
        filt_data_mat = int16([]);
        disp('Filtering.')
        for b = 1:opt.numChannels
            % Proceed with filter
            [filt_data_mat{b,1}, ~, ~] = bandFilter(double(data_mat{b,1}), [], opt.highpass, opt.sampleRate);
            filt_data_mat{b,1} = int16(filt_data_mat{b,1});
        end
    else
        filt_data_mat = data_mat;
    end
    clear data_mat

    % Write bin file
    fwrite(fidDataMat, cell2mat(filt_data_mat), 'int16');
    fclose(fidDataMat);

%% DF1 FORMAT
elseif strcmp(opt.ext, 'DF1')
    disp('Format is BLOCK.')

    % Create variable to store all data
    data_mat    = [];

    % Set stream to continuous neural signals
    stream      = 1;
    
    % Allocate data to its respective single-channel file
    % Open each neural data file, resize data for detection with Kilosort,
    % Allocate data to its respective single-channel file.
    for i = 1:length(opt.myFiles)
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
            tempdata = Deuteron_extractData(stream, fid, opt);
        fclose(fid);

        % Convert ADC steps into microvolts, so conversion to int16 is
        % possible without loss.
        tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);

        try
            data_mat = [data_mat tempdata];
        catch
            disp('found non-matching file')
        end
    end
    clear tempdata fid

    % Reshape to channels x samples.
    data_mat = reshape(data_mat, opt.numChannels, []);

    if opt.set_filter
        % To keep memory usage low, we proceed in a channel by channels basis
        filt_data_mat = int16([]);
        disp('Filtering.')
        for b = 1:opt.numChannels
            % Proceed with filter
            [filt_data_mat(b,:), filt.high.a, filt.high.b] = bandFilter(double(data_mat(b,:)), [], opt.highpass, opt.sampleRate);
            filt_data_mat(b,:) = int16(filt_data_mat(b,:));
        end
    else
        filt_data_mat = data_mat;
    end

    % % distribute each row of data to its respective single-channel file
    % if opt.h5
    %     for b = 1:opt.numChannels
    %         h5write(filename,         ... % filename.
    %                 dataset,          ... % dataset name.
    %                 data_mat(b,:),        ... % channel b, complete
    %                 [b b],                ... %   into slot b
    %                 [1 size(data_mat,2)]);    %   as long as it is.
    %     end
    % end

    % Write bin file
    fwrite(fidDataMat, filt_data_mat, 'int16');
    fclose(fidDataMat);
end  

end