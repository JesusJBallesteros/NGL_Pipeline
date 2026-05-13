function metaData = Deuteron_GetMetaData(info)
% Deuteron_GetMetaData  Return hardcoded acquisition parameters for a Deuteron format.
%
% PURPOSE:
%   Provides the hardware constants needed to interpret raw ADC values:
%   channel count, ADC bit depth, voltage resolution (V/bit), and sample rate.
%   Values are taken directly from Deuteron hardware documentation.
%   Channel count for DF1 is returned as [] because the actual value must be
%   inferred later from the event log (via Deuteron_ExtractEvents / opt.channelOrder).
%
% USAGE:
%   metaData = Deuteron_GetMetaData(info)
%   Called from findSessions during session inventory.
%
% INPUT:
%   info  - struct with field:
%             .fileformat   'DF1' (current) or 'DT2' (legacy)
%
% OUTPUT:
%   metaData  - struct with fields:
%                 .numChannels    channel count ([] for DF1; 32 for DT2)
%                 .numADCBits     ADC bit depth (16 for all formats)
%                 .voltageRes     V/bit scaling (1.95e-7 for DF1; 0.2e-6 for DT2)
%                 .fSample        sample rate in Hz (32000 for all formats)
%                 .ext            file format string (passed through)
%
% Last modified 08.05.2026 (Jesus)

switch (info.fileformat)
    case 'DF1' % NEW FORMAT, predominant
        numberOfChannels = []; % Could be 32 or 64 
        numberOfADCBits = 16;
        voltageResolution = 1.95e-7;
        fSample = 32000;
    case 'DT2' % OLD FORMAT, deprecating
        numberOfChannels = 32;
        numberOfADCBits = 16;
        voltageResolution = 0.2e-6;
        fSample = 32e3;

%  DEPR    % Following are never seen formats
%     case 'DT4'
%         numberOfChannels = [];
%         numberOfADCBits = 16;
%         voltageResolution = 0.2e-6;
%         fSample = 32e3;
%     case 'DT8'
%         numberOfChannels = [];
%         numberOfADCBits = 15;
%         voltageResolution = 0.42e-6;
%         fSample = 4e3;
%     case 'DAT'
%         numberOfChannels = [];
%         numberOfADCBits = 12;
%         voltageResolution = 3.3e-6;
%         fSample = 31.25e3;
% DEPR

    otherwise
        error('Invalid file extension. Please choose a file with an extension of ''DAT'', ''DT2'', ''DT4'', or ''DT8''.')   
end

metaData = struct('numChannels', numberOfChannels, 'numADCBits', numberOfADCBits,... 
    'voltageRes', voltageResolution, 'fSample', fSample, 'ext', info.fileformat);

end