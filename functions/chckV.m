function [info] = chckV()
% Onc in the session folder, checks for existence of any of the following
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
% 27.01.2023. Jesus

    if isfile('EVENTLOG.NLE') 
        % For this format, we list the files with neural data and extract some
        % metadata with the function 'Deuteron_GetMetaData'. Apparently, it can
        % be DT2/4/8 or DAT. We only have .DT2
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
            info.numChannels        = 32; % Coded as default here. If necessary, overrided later on.
            info.numADCBits         = metaData.numADCBits;
            info.voltageRes         = metaData.voltageRes;
            info.sampleRate         = metaData.fSample;
            info.HDF5chunkSize      = 300*info.sampleRate;
    
        else % Still empty for some reason
            warning('Something went wrong with this Deuteron flat format session.')
            info.fileformat    = 'NAN'; % Flag for error with the file format
            info.numChannels   = [];
            info.numOfADCBits  = [];
            info.voltageRes    = [];
            info.sampleRate    = [];
            info.HDF5chunkSize = [];
            return
        end
    
    elseif isfile('EVENT000.DF1') 
        % For this format, we list the files with neural data and extract some
        % metadata with the function 'Deuteron_GetMetaData'. No further check needed.
        info.fileformat = 'DF1';
        info.files = dir(['*.' info.fileformat]);
        if ~isempty(info.files)
            % Extract metadata
            metaData           = Deuteron_GetMetaData(info);
            info.numChannels   = metaData.numChannels;
            info.numADCBits    = metaData.numADCBits;
            info.voltageRes    = metaData.voltageRes;
            info.sampleRate    = metaData.fSample;
            info.HDF5chunkSize = 300*info.sampleRate;
    
        else
            warning('Something went wrong with this Deuteron block format session.')
            info.fileformat    = 'NAN'; % Flag for error with the file format
            info.numChannels   = [];
            info.numADCBits    = [];
            info.voltageRes    = [];
            info.sampleRate    = [];
            info.HDF5chunkSize = [];
            return
        end
    
    elseif isfile('info.rhd') 
       % Here we look for files of each of the bandpass to use as source. 
       %  'amp' should always exist. Would be used as ultimate source of
       %  data if 'low' does not. If 'low' exist, the loop breaks and takes
       %  the indexed file list with such extension.
       for i=1:2
           info.files = dir('low*.dat');
           info.nfiles = length(info.files);
           if ~isempty(info.files)
              break
           else
              info.files = dir('amp*.dat');
           end
       end
        
       % How many files exist for this sessions. If we have several types 
       % with file per channel format, 'fileperch' applies anyways.
       % If there are many files, is 'fileperch', if there is only one, is 'filepertype'
       if info.nfiles > 1, info.fileformat = 'fileperch';
       else,               info.fileformat = 'filepertype';  end
    
    else
        info.files = dir('*xdat.json'); % list files matching Allego's
        if ~isempty(info.files)
            info.fileformat = 'Allego';
            % TODO: Figure out how to work with this files.
            % (most likely, after Allego's self preprocessing tool?)

        else
            warning('The format of the sessions could not be determined.')
            info.fileformat    = 'NAN'; % Flag for error with the file format
            info.numChannels   = [];
            info.numADCBits    = [];
            info.voltageRes    = [];
            info.sampleRate    = [];
            info.HDF5chunkSize = [];
        end
    end
end