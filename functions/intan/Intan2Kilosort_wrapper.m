function Intan2Kilosort_wrapper(sessions, varargin)
% Intan2Kilosort_wrapper  Create flat binary .bin file for Kilosort spike sorting.
%
% PURPOSE:
%   Wrapper for .bin file creation. Reads INTAN amp*.dat (or high*.dat)
%   files, scales ADC values to µV (×0.195), optionally filters, and writes
%   a flat int16 interleaved binary file (channels × samples) for Kilosort.
%   Dispatches to format-specific implementations based on sessions.info.fileformat.
%
% USAGE:
%   Intan2Kilosort_wrapper(sessions, opt)
%   Intan2Kilosort_wrapper(sessions)    % opt defaulted
%
% INPUTS:
%   sessions  - struct from input.sessions(x), must contain:
%                 .info.fileformat        ('fileperch' | 'filepertype' | 'tradFormat')
%                 .info.nChannels         electrode count from header
%                 .info.nfiles            number of .dat files found
%                 .info.amplifier_sample_rate (Hz)
%   opt       - options struct; relevant fields:
%                 .myFiles      list of .dat files to read
%                 .set_filter   1 = apply detrend+highpass+lowpass; 0 = skip
%                 .numChannels  overrides sessions.info.nChannels if mismatch
%                 .sampleRate   sample rate (Hz)
%                 .StpSz        chunk size in samples (default = 300 × sampleRate)
%                 .num_samples  total samples per channel (set per format)
%
% OUTPUT:
%   <session>.bin file written to opt.FolderProcDataMat (set by prepforsession)
%   File format: int16, [nChannels × nSamples]
%
% CALLS:
%   Intan2Kilosort_fileperch, Intan2Kilosort_filepertype, Intan2Kilosort_tradFormat
%
% Last modified 07.05.2026 (Jesus)

if nargin < 2,  opt = struct();
else, opt = varargin{1};
end

if contains(sessions.info.fileformat, 'fileper')
    % Collect parameters that not need to necessarily defaulted to a given value. 
    % To proceed, list files depending on filetype.
    opt.myFiles = dir('high*.dat');
    opt.set_filter = 0;
    
    % If no highpass files are found, it will use the raw 'amp' data and filtering will be applied.
    if isempty(opt.myFiles)
        opt.myFiles = dir('amp*.dat');
        opt.set_filter = 1;
    end
    
    % How many channels, from Intan_hdr.
    opt.numChannels = sessions.info.nChannels;
    
    % Check that numChannels and the number of files coincide. Fix if necessary
    if sessions.info.nfiles < opt.numChannels
        % Report
        disp('Found less channel files than expected by header. Using number of files as truth.')
        opt.numChannels = sessions.info.nfiles;
        sessions.info.nChannels = opt.numChannels;
    end
    
    % Sample rate, from Intan_hdr.
    opt.sampleRate  = sessions.info.amplifier_sample_rate;
    
    % ChunkSize (e.g. 5 minutes = 300 s @30000 Hz = 9600000 samples).
    opt.StpSz = 300*opt.sampleRate;
        
    % To obtain the number of samples per file, first read file info. Can be a 
    % file per channel or only one for all, it does not matter. Read the first.
    fileinfo = dir(opt.myFiles(1).name);

    if strcmp(sessions.info.fileformat,'filepertype')
        % To get the number of samples, divide the file size by number of 
        % channels, times bytes that each int16 word takes (int16 = 2 bytes).
        opt.num_samples = fileinfo.bytes/(opt.numChannels * 2); 
        Intan2Kilosort_filepertype(opt);
    
    elseif strcmp(sessions.info.fileformat,'fileperch')
        % To get the number of samples, divide the file size by the bytes 
        % each int16 word takes (int16 = 2 bytes).
        opt.num_samples = fileinfo.bytes/2;
        Intan2Kilosort_fileperch(opt);
    end
    
elseif strcmp(sessions.info.fileformat,'tradFormat')
    opt.myFiles = sessions.info.files;
    opt.set_filter = 0;
        
    % How many channels, from Intan_hdr.
    opt.numChannels = sessions.info.nChannels;

    Intan2Kilosort_tradFormat(opt);
end

end