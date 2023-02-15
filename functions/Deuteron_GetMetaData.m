function metaData = Deuteron_GetMetaData(info)
% GetMetaData Checks for what type of file this is.
%   Builds a struct with essential information for data processing
%   This information is found in the instruction manual of each type of logger.
    switch (info.fileformat)
        case 'DT2'
            numberOfChannels = [];
            numberOfADCBits = 16;
            voltageResolution = 0.2e-6;
            fSample = 32e3;
        case 'DT4'
            numberOfChannels = [];
            numberOfADCBits = 16;
            voltageResolution = 0.2e-6;
            fSample = 32e3;
        case 'DT8'
            numberOfChannels = [];
            numberOfADCBits = 15;
            voltageResolution = 0.42e-6;
            fSample = 4e3;
        case 'DAT'
            numberOfChannels = [];
            numberOfADCBits = 12;
            voltageResolution = 3.3e-6;
            fSample = 31.25e3;
        case 'DF1'
            numberOfChannels = [];
            numberOfADCBits = 16;
            voltageResolution = 1.95e-7;
            fSample = 32000;
        otherwise
            error('Invalid file extension. Please choose a file with an extension of ''DAT'', ''DT2'', ''DT4'', or ''DT8''.')
        
    end
    metaData = struct('numChannels', numberOfChannels, 'numADCBits', numberOfADCBits,... 
        'voltageRes', voltageResolution, 'fSample', fSample, 'ext', info.fileformat);

end