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
% Last modified 08.05.2026 (Jesus)

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