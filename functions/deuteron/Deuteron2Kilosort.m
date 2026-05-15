function Deuteron2Kilosort(opt)
% Deuteron2Kilosort  Convert Deuteron DF1 files to a Kilosort-ready .bin file.
%
% PURPOSE:
%   Reads NEUR*.DF1 files block by block via Deuteron_extractData, converts
%   ADC counts to µV (int16), applies optional CAR and DC removal, applies
%   a high-pass filter if opt.highpass is set, and writes the result as a
%   single int16 binary file (<SavFileName>.bin) in channels × samples layout
%   as required by Kilosort.
%
% USAGE:
%   Deuteron2Kilosort(opt)
%   Called from Deuteron_PipelineWrapper when opt.bin = true.
%
% INPUTS:
%   opt  - options struct; relevant fields:
%            .ext               file format ('DF1')
%            .myFiles           dir-struct of NEUR*.DF1 files
%            .numChannels       electrode count (set by EventProcess via opt.channelOrder)
%            .sampleRate        raw sample rate (Hz, typically 32000)
%            .highpass          high-pass cutoff (Hz); [] = skip filtering
%            .CAR               common-average re-referencing flag (0 = off)
%            .voltageResolution, .offset  ADC-to-µV conversion parameters
%            .PathRaw           raw data folder path
%            .FolderProcDataMat output folder for the .bin file
%            .SavFileName       session name used for output filename
%
% OUTPUT:
%   <SavFileName>.bin written to opt.FolderProcDataMat; int16, channels × samples,
%   little-endian. Skipped if a non-empty .bin already exists (re-run safety).
%
% NOTES:
%   - For Kilosort 4 a high-pass filter (300–600 Hz) is strongly recommended;
%     set opt.highpass accordingly (e.g. 300).
%   - ADC conversion: int16((voltageResolution × (data − offset)) × 1e6)
%
% CALLS:
%   Deuteron_extractData, ft_preproc_rereference, ft_preproc_detrend,
%   ft_preproc_highpassfilter
%
% Last modified 13.05.2026 (Jesus)

%% Check bin file existence and completion.
% If an error happens during processing, the bin file persists created but with zero size
if isfile(fullfile(opt.FolderProcDataMat, opt.SavFileName + ".bin"))
    bininfo = dir(fullfile(opt.FolderProcDataMat, opt.SavFileName + ".bin"));
    if bininfo.bytes > 0
        disp('An actual .bin file exists in this directory. Conversion Skipped.')
        return
    end
    disp('A failed .bin file was found in this directory. Will be overwritten.')
end

%% Create an new .bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,opt.SavFileName + ".bin"), 'a'); 

% %% DT2 format conversion.
% if strcmp(opt.ext, 'DT2')
%     disp('Format is DT2. Deprecating.')
% 
%     % Prepare a cell array with nChannels
%     data_mat = cell(opt.numChannels,1);
% 
%     % Initialize variable to chunk the writting. 
%     indexPos  = 0;
% 
%     % Go over individual files
%     for i = 1:length(opt.myFiles)
% 
%         % Neural data points are unsigned 16 bit words. Read.
%         fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
%             tempdata = fread(fid, 'uint16');
%         fclose(fid);
% 
%         % Remove trailing zeroes in last file
%         if i == length(opt.myFiles)
%             tempdata(tempdata==0) = [];
%         end
% 
%         % Because all channels come concatenated, we need to reshape as
%         % channels x samples, knowing the nChannels.
%         tempdata = reshape(tempdata', opt.numChannels, []); 
% 
%         % Convert ADC steps into microvolts, so conversion to int16 is
%         % possible without loss.
%         tempdata = int16((opt.voltageResolution * (tempdata - opt.offset)) * 1000000);
% 
%         % Get nSamples coming from this file. Should stay constant, until
%         % last file which normally will be smaller.
%         nSamples = size(tempdata,2);
% 
%         % Distribute each channel to its respective slot in new array
%         for b = 1:opt.numChannels            
%             % Write channel into general cell array.
%             data_mat{b,1} = [data_mat{b,1} tempdata(b,:)];
%         end
% 
%         % Index for next chunk's starting sample
%         indexPos = indexPos + nSamples;
% 
%     end
%     clear tempdata fid
% 
%     % Convert cell to array, and to single
%     data_mat = cell2mat(data_mat);
% 
% end 

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
            if i==length(opt.myFiles)
                disp('The last data chunk could not be concatenated')
            else
                disp('A non-matching file was found during data creation')
            end
        end
    end
    clear tempdata fid

    % Reshape to channels x samples.
    data_mat = reshape(data_mat, opt.numChannels, []);
end

%% Save raw data matrix for NWB conversion (before any filtering).
% This file is consumed by Deuteron2NWB and deleted after a successful NWB
% write. If opt.doNWB is false the file is never created.
if opt.doNWB
    rawMatPath = fullfile(opt.FolderProcDataMat, [char(opt.SavFileName), '_raw.mat']);
    if ~isfile(rawMatPath)
        disp('Saving raw data matrix for later NWB conversion...')
        sampleRate_raw = opt.sampleRate; %  saved alongside data
        save(rawMatPath, 'data_mat', 'sampleRate_raw', '-v7.3');
        disp('Raw data matrix saved.')
    end
end

%% Common methods of preprocessing. Re-Referencing, DC substraction and filter.
if opt.CAR
    % In principle, data from a single HS on a single region.
    disp('Re-referencing by Common Average Referencing (CAR).')
    data_mat = ft_preproc_rereference(data_mat, 'all', 'median');
end

% Enforce int16
filt_data_mat = int16([]);

% Keep memory usage low doing one channel at a time.
for b = 1:opt.numChannels
    fprintf('- Channel %d of %d.\n', b, opt.numChannels);
    
    % Detrend channel (remove DC)
    disp('Detrending...')
    filt_data_mat(b,:) = ft_preproc_detrend(data_mat(b,:));

    % Highpass channel (Butterwort, 4th order, back&forth)
    if ~isempty(opt.highpass)
        txt = sprintf('Highpass at %d Hz.\n', opt.highpass);
        fprintf(txt);
        [filt_data_mat(b,:), ~, ~] = ft_preproc_highpassfilter(data_mat(b,:), opt.sampleRate, opt.highpass, 4, 'but', 'twopass');
    end       
end

%% Write bin file.
fwrite(fidDataMat, filt_data_mat, 'int16');
fclose(fidDataMat);

end  