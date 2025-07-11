function [info] = chckV(varargin)
% Once in the session folder, checks for existence of any of the following
% characteristic files, determinant of the recording format:
%   'EVENTLOG.NLE': characteristic of Deuteron flat format
%   'EVENT000.DF1': characteristic of Deuteron block format
%   'info.rhd': characteristic of INTAN file-per-channel or file-per-type
%                formats. They differenciate based on number of files.
%   '*xdat.json': Base format from Allego. Prob can be converted to .H5 and
%                take it from there as well.
%
% Any other would just throw a warning and skip the session.
%
% Calls to: 'Deuteron_GetMetaData'
%
% 07.11.2023. Jesus

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
        formatis = 0;
        err = 1; 
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
       %  'amp' should always exist. Would be used as ultimate source of
       %  data if 'low' does not. If 'low' exist, the loop breaks and takes
       %  the indexed file list with such extension.
       
       metaData = readstruct('settings.xml');
       for b = 1:size(metaData.SignalGroup,2)
        if strlength(metaData.SignalGroup(b).PrefixAttribute)==1
            nchan(b) = size(metaData.SignalGroup(b).Channel,2)-3-1; % always 3 AUX + 1 VDD
        end
       end
       info.nChannels     = sum(nchan);
           
       % This info is extracted later on, 'findSetting.m'
           % info.amplifier_sample_rate = metaData.SampleRateHertzAttribute;
       
       for i = 1:2
           info.files = dir('low*.dat');
           if ~isempty(info.files)
              info.bandpass = 'low';
              break
           else
              info.files = dir('amp*.dat');
              info.bandpass = 'amp';
           end
       end
       info.nfiles = length(info.files);
        
       % How many files exist for this sessions. If we have several types 
       % with file per channel format, 'fileperch' applies anyways.
       % If there are many files, is 'fileperch', if there is only one, is 'filepertype'
       if info.nfiles > 1, info.fileformat = 'fileperch';
       else,               info.fileformat = 'filepertype';  end
    case 4
        info.fileformat = 'FieldTrip';
        % TODO: Figure out how to work with this files.
        % (most likely, after Allego's self preprocessing tool?)
    
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