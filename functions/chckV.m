function [info] = chckV(varargin)
% chckV  Detect recording format and extract basic session metadata.
%
% PURPOSE:
%   Called from prepforsession while in the session's raw-data directory.
%   Identifies the recording format by characteristic file presence, then
%   reads enough metadata to populate the 'info' struct used by all
%   downstream pipeline functions.
%
% USAGE:
%   info = chckV()
%   Must be called from within the session's raw-data folder (cd there first).
%
% FORMAT DETECTION (in priority order):
%   EVENTLOG.NLE present  → Deuteron flat format   (DT2)
%   EVENT000.DF1 present  → Deuteron block format  (DF1)
%   info.rhd present      → INTAN RHX format       (fileperch | filepertype)
%     fileperch  : multiple amp*.dat files (one per channel)
%     filepertype: single amp*.dat file (all channels interleaved)
%   *FTcont.mat present   → Pre-processed FieldTrip (FieldTrip)
%   settings.xml present  → INTAN traditional/legacy format (tradFormat)
%   None of the above     → error; info.fileformat = 'NAN'
%
% OUTPUT:
%   info - struct with fields (populated vary by format):
%     .fileformat           (char)  one of: 'DT2','DF1','fileperch',
%                                   'filepertype','tradFormat','FieldTrip','NAN'
%     .files                (dir-struct) relevant data files
%     .nfiles               (scalar) number of data files
%     .nChannels            (scalar) active recording channel count
%     .amplifier_sample_rate (scalar, Hz)
%     .numADCBits           (scalar) ADC resolution bits
%     .voltageRes           (scalar) voltage resolution
%     .HDF5chunkSize        (scalar) chunk size for HDF5 output (300 s × Fs)
%     .bandpass             (char, INTAN only) 'low' or 'amp'
%     .INTAN_hdr            (struct, added by findSetting) full RHD header
%
% CALLS:
%   Deuteron_GetMetaData (for Deuteron formats)
%   findSetting (called separately from prepforsession after chckV returns)
%
% Last modified 07.11.2023 (Jesus)

err = 0;
% Evaluate if file exist with full name and assign case.
if isfile('EVENTLOG.NLE') 
    formatis = 1;
elseif isfile('EVENT000.DF1') 
    formatis = 2; 
elseif isfile('info.rhd') 
    formatis = 3; 
else
    info.files = dir('*FTcont.mat'); % list files matching Allego's
    if ~isempty(info.files)
        formatis = 4; 
    else
        if isfile('settings.xml') 
            % Possibly, old data in traditional format. Here for specific
            % case, and to be deprecated
           formatis = 5; 
        else 
           formatis = 0;
           err = 1;
        end
    end
end

switch formatis
    case 1
        % For this format, we list the files with neural data and generate 
        % metadata with the function 'Deuteron_GetMetaData'.
        info.fileformat = 'DT?'; % First try most common 'DT?'
        info.files = dir(['*.' info.fileformat]); % list files
        
        if isempty(info.files) 
            info.fileformat = 'DAT'; % Try with 'DAT'
            info.files = dir(['*.' info.fileformat]); % list files again
        end
    
        if ~isempty(info.files)
            [~,~,ext] = fileparts(info.files(1).name); % extract actual extension
            info.fileformat = ext(2:end); % remove dot
    
            % Corresponding files: get all DT2 files in folder.
            % Get meta data from Deuteron:
            % Checks which type of logger was used and sets some parameters:
            metaData                = Deuteron_GetMetaData(info);
            info.nChannels          = metaData.numChannels; % Coded as default here. If necessary, overrided later on.
            info.numADCBits         = metaData.numADCBits;
            info.voltageRes         = metaData.voltageRes;
            info.amplifier_sample_rate         = metaData.fSample;
            info.HDF5chunkSize      = 300*info.amplifier_sample_rate;
    
        else % Still empty for some reason
            err = 1;
        end
    case 2
        % For this format, we list the files with neural data and extract some
        % metadata with the function 'Deuteron_GetMetaData'. No further check needed.
        info.fileformat = 'DF1';
        info.files = dir(['N*.' info.fileformat]);
        if ~isempty(info.files)
            % Extract metadata
            metaData           = Deuteron_GetMetaData(info);
            info.nChannels     = metaData.numChannels;
            info.numADCBits    = metaData.numADCBits;
            info.voltageRes    = metaData.voltageRes;
            info.amplifier_sample_rate    = metaData.fSample;
            info.HDF5chunkSize = 300*info.amplifier_sample_rate;
        else
            err = 1;
        end

    case 3
        % Here we look for files of each of the bandpass to use as source. 
        % 'amp' should always exist. Would be used as primary source of data
        metaData = readstruct('settings.xml');
        for b = 1:size(metaData.SignalGroup,2)
            if strlength(metaData.SignalGroup(b).PrefixAttribute)==1
                % Read only info from analog inputs, labeled "A", "B", etc
                if size(metaData.SignalGroup(b).Channel,2)>68
                    % 2x32 channel HeadStage in same port
                    nchan(b) = size(metaData.SignalGroup(b).Channel,2)-6-2; % 6 AUX + 2 VDD
                else % either 32ch or 64ch HS. In any case,
                    nchan(b) = size(metaData.SignalGroup(b).Channel,2)-3-1; % 3 AUX + 1 VDD
                end
            end
        end

        % Add all channels detected across ports
        info.nChannels     = sum(nchan);
                  
       for i = 1:2
           % look for low band data first (legacy need from past for LFP)
           info.files = dir('low*.dat');
           if ~isempty(info.files)
              % if they exists, will work with them
              info.bandpass = 'low';
              break
           else
              % If not, expected standard now, go for wideband
              info.files = dir('amp*.dat');
              info.bandpass = 'amp';
           end
       end
       % Check how many of those files exist
       info.nfiles = length(info.files);
        
       % If there are many files, is 'fileperch', if there is only one, is 'filepertype'
       % If we have several types with file per channel format, 'fileperch' applies anyways.
       if info.nfiles > 1, info.fileformat = 'fileperch';
       else,               info.fileformat = 'filepertype';  end

    case 4
        info.fileformat = 'FieldTrip';
        % TODO: Figure out how to work with this files.
        % (most likely, after Allego's self preprocessing tool?)

    case 5
       metaData = readstruct('settings.xml');
       for b = 1:size(metaData.SignalGroup,2)
          if strlength(metaData.SignalGroup(b).PrefixAttribute)==1
            % Read only info from analog inputs, labeled "A", "B", etc
            if size(metaData.SignalGroup(b).Channel,2)>68
                % 2x32 channel HeadStage in same port
                nchan(b) = size(metaData.SignalGroup(b).Channel,2)-6-2; % 6 AUX + 2 VDD
            else % either 32ch or 64ch HS. In any case,
                nchan(b) = size(metaData.SignalGroup(b).Channel,2)-3-1; % 3 AUX + 1 VDD
            end
          end
       end
       info.nChannels     = sum(nchan);                  
       info.files = dir('*.rhd');
       info.nfiles = length(info.files);
       info.fileformat = 'tradFormat';
       
    case 0
        % Invalid formats or error
        warning('The format of the sessions could not be determined.')
end

if err
    warning('Something went wrong with this session.')
    info.fileformat    = 'NAN'; % Flag for error with the file format
    info.nChannels     = [];
    info.numOfADCBits  = [];
    info.voltageRes    = [];
    info.amplifier_sample_rate    = [];
    info.HDF5chunkSize = [];
    return
end

end