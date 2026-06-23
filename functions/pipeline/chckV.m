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
%     .hardwareFilters      (struct, INTAN only) settings.xml-derived
%                                   record of the hardware filter chain
%                                   .notchFreq (Hz, 0=off, NaN=unset)
%                                   .hpfFreq / .lpfFreq / .dspFreq
%                                   (NaN where the attribute wasn't set).
%                                   Attribute names verified against RHX 3.4
%                                   (NotchFilterFreqHertz, Desired{Lower,Upper}
%                                   BandwidthHertz, DesiredDSPCutoffFreqHertz).
%     .recordingTime        (char, INTAN only) reserved; currently always ''.
%                                   settings.xml does not carry the wall-
%                                   clock recording time and info.rhd has
%                                   no absolute timestamp either; we leave
%                                   the field for a future source (file
%                                   mtime / folder name parser).
%     .recordingNotes       (cellstr, INTAN only) non-empty Note1/Note2/Note3
%                                   from settings.xml GeneralConfig.
%     .rhxVersion           (char, INTAN only) IntanRHX Version root attribute,
%                                   e.g. '3.4.0'.
%     .controllerType       (char, INTAN only) IntanRHX Type root attribute,
%                                   e.g. 'ControllerRecordUSB3'.
%     .INTAN_hdr            (struct, added by findSetting) full RHD header
%
% CALLS:
%   Deuteron_GetMetaData (for Deuteron formats)
%   findSetting (called separately from prepforsession after chckV returns)
%
% Last modified 18.06.2026 (Jesus) - case 3: read fileformat from
%                                     settings.xml GeneralConfig.FileFormat.
%                                     Verified attribute names against a real
%                                     RHX 3.4 settings.xml + info.rhd:
%                                     SampleRateHertz is at the root; filter
%                                     attrs use a 'Hertz' suffix; Notch can
%                                     be the string "None". Added
%                                     info.controllerType from root Type.
%                                     info.recordingTime is now a reserved
%                                     empty char (no source in settings.xml
%                                     or info.rhd).

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
        % INTAN RHX with info.rhd present. Use settings.xml as the
        % primary source of metadata; fall back to file-count heuristics
        % only when an attribute is missing (older RHX versions).
        metaData = readstruct('settings.xml');

        % --- File format from GeneralConfig.FileFormat -------------
        %   'OneFilePerChannel'    -> 'fileperch'    (amp_*.dat per channel)
        %   'OneFilePerSignalType' -> 'filepertype'  (single amp.dat)
        %   'Traditional'          -> 'tradFormat'   (handled by case 5;
        %                                              fall through if we
        %                                              get here defensively).
        %   missing / unknown      -> warn and fall back to amp*.dat count.
        info.fileformat = '';
        if isfield(metaData, 'GeneralConfig') ...
                && isfield(metaData.GeneralConfig, 'FileFormatAttribute')
            fmtAttr = char(string(metaData.GeneralConfig.FileFormatAttribute));
            switch lower(fmtAttr)
                case 'onefileperchannel'
                    info.fileformat = 'fileperch';
                case 'onefilepersignaltype'
                    info.fileformat = 'filepertype';
                case 'traditional'
                    info.fileformat = 'tradFormat';
                otherwise
                    warning('NGL:chckV:unknownFmtAttr', ...
                        ['settings.xml GeneralConfig.FileFormat=''%s'' is ', ...
                         'unrecognised; falling back to amp*.dat count heuristic.'], ...
                        fmtAttr);
            end
        else
            warning('NGL:chckV:missingFmtAttr', ...
                ['settings.xml has no GeneralConfig.FileFormat attribute; ', ...
                 'falling back to amp*.dat count heuristic. ', ...
                 'This suggests an older RHX version — consider re-exporting.']);
        end

        % --- Amplifier sample rate (settings.xml root attribute) ------
        % Verified against RHX 3.4 (June 2026): SampleRateHertz lives on
        % the document root <IntanRHX>, NOT inside GeneralConfig.
        % findSetting (called later from prepforsession) also fills this
        % from info.rhd; we set it here so any caller that runs chckV in
        % isolation gets it directly.
        if isfield(metaData, 'SampleRateHertzAttribute')
            info.amplifier_sample_rate = double(metaData.SampleRateHertzAttribute);
        elseif isfield(metaData, 'GeneralConfig') ...
                && isfield(metaData.GeneralConfig, 'SampleRateHertzAttribute')
            info.amplifier_sample_rate = double(metaData.GeneralConfig.SampleRateHertzAttribute);
        end
        if isfield(info, 'amplifier_sample_rate') && ~isempty(info.amplifier_sample_rate)
            info.HDF5chunkSize = 300 * info.amplifier_sample_rate;
        end

        % --- Provenance metadata (user notes, recorder version, type).
        %     Pass-through fields; not consumed by analysis but stamped
        %     onto preprocInfo for forensic continuity. The wall-clock
        %     recording time has NO source in settings.xml or info.rhd
        %     (verified against RHX 3.4); leave it empty until a better
        %     source is wired up (file mtime / folder name).
        info.recordingTime   = '';   % no source in settings.xml; reserved
        info.recordingNotes  = {};   % cellstr of non-empty Note1/Note2/Note3 (GeneralConfig)
        info.rhxVersion      = '';   % e.g. '3.4.0' (root attribute)
        info.controllerType  = '';   % e.g. 'ControllerRecordUSB3' (root attribute)
        meta_gc = struct();
        if isfield(metaData, 'GeneralConfig'), meta_gc = metaData.GeneralConfig; end

        for n = 1:3
            attr = sprintf('Note%dAttribute', n);
            v    = localPickAttr({meta_gc, metaData}, attr);
            if ~isempty(v), info.recordingNotes{end+1, 1} = v; end
        end

        info.rhxVersion     = localPickAttr({metaData, meta_gc}, 'VersionAttribute');
        info.controllerType = localPickAttr({metaData, meta_gc}, 'TypeAttribute');

        % --- Hardware filter settings (informational; downstream sanity
        %     checks against opt.lowpass / opt.highpass / opt.linefilter
        %     should consult these). Attribute names verified against
        %     RHX 3.4 — all carry a '*Hertz' suffix and live in
        %     GeneralConfig. NotchFilterFreqHertz can be the string
        %     "None"/"Off" rather than 0, which we coerce explicitly.
        info.hardwareFilters = struct( ...
            'notchFreq', NaN, ...  % 50 / 60 Hz notch (0 = off, NaN = unset)
            'hpfFreq',   NaN, ...  % desired hardware HPF (lower bandwidth)
            'lpfFreq',   NaN, ...  % desired hardware LPF (upper bandwidth)
            'dspFreq',   NaN);     % desired DSP cutoff (settling HPF)
        if isfield(metaData, 'GeneralConfig')
            gc = metaData.GeneralConfig;
            if isfield(gc, 'NotchFilterFreqHertzAttribute')
                info.hardwareFilters.notchFreq = localNotchToHz(gc.NotchFilterFreqHertzAttribute);
            end
            if isfield(gc, 'DesiredLowerBandwidthHertzAttribute')
                info.hardwareFilters.hpfFreq = double(gc.DesiredLowerBandwidthHertzAttribute);
            end
            if isfield(gc, 'DesiredUpperBandwidthHertzAttribute')
                info.hardwareFilters.lpfFreq = double(gc.DesiredUpperBandwidthHertzAttribute);
            end
            if isfield(gc, 'DesiredDSPCutoffFreqHertzAttribute')
                info.hardwareFilters.dspFreq = double(gc.DesiredDSPCutoffFreqHertzAttribute);
            end
        end

        %  Per-port channel counts 
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

       % Fall-back: if GeneralConfig.FileFormat couldn't be read above,
       % use the legacy amp*.dat count heuristic. fileperch when >1 file,
       % filepertype when ==1 file.
       if isempty(info.fileformat)
           if info.nfiles > 1, info.fileformat = 'fileperch';
           else,               info.fileformat = 'filepertype';  end
       end

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

% ----------------------------------------------------------------------
function f = localNotchToHz(raw)
% RHX writes NotchFilterFreqHertz as the string "None" / "Off" when the
% notch is disabled, or a numeric "50" / "60" when on. Coerce to Hz with
% 0 for the off-states (so downstream "notch active?" checks are simple
% f > 0) and NaN only when we genuinely can't tell.
    s = strtrim(char(string(raw)));
    if isempty(s)
        f = NaN;
        return
    end
    switch lower(s)
        case {'none','off','0'}
            f = 0;
        otherwise
            f = str2double(s);
    end
end

% ----------------------------------------------------------------------
function v = localPickAttr(structs, attrName)
% Return the first non-empty value of `attrName` across the candidate
% structs (in order). Coerces to char and trims. Empty char when none.
    v = '';
    for k = 1:numel(structs)
        s = structs{k};
        if isstruct(s) && isfield(s, attrName)
            raw = s.(attrName);
            if ~isempty(raw)
                v = strtrim(char(string(raw)));
                if ~isempty(v), return, end
            end
        end
    end
end