function input = findSetting(input)
% findSetting  Read INTAN RHD header and settings.xml for the current session.
%
% PURPOSE:
%   Called from prepforsession (while in the session's raw-data folder) for
%   INTAN 'fileperch' and 'filepertype' formats. Reads 'info.rhd' via
%   mod_read_Intan_RHD2000_file and, if present, parses 'settings.xml' to
%   extract amplifier sample rate, lowpass downsample factor, firmware
%   version, and the full parameter tree. Results are merged into
%   input.sessions(x).info so all downstream functions have header access.
%
% USAGE:
%   input = findSetting(input)
%   Requires input.run to be set and CWD to be the session's raw-data folder.
%
% INPUT:
%   input  - struct with input.run and input.sessions(input.run(1)).info
%            already populated by chckV()
%
% OUTPUT:
%   input  - input.sessions(input.run(1)).info extended with:
%              .INTAN_hdr             full RHD header struct
%              .nChannels             active channel count (from header)
%              .amplifier_sample_rate (Hz) from settings.xml or header
%              .lowpass_downsample    downsample factor (or [] if unavailable)
%              .lowpass_sample_rate   amplifier_sample_rate / downsample factor
%              .Version               firmware/software version string
%              .Name                  system name
%              .Children              raw settings.xml parameter tree
%
% NOTES:
%   - settings.xml is present in recordings made with INTAN RHX software.
%   - Older datasets (RHD GUI) may only have info.rhd; lowpass fields will
%     be [] in that case.
%   - The helper function getField (local) traverses nested structs via a
%     dot-path string, e.g. 'Children(2).Attributes(85).Value'.
%
% CALLS:
%   mod_read_Intan_RHD2000_file, parseXML
%
% Last modified 05.05.2026 (Jesus)
sessions = input.sessions(input.run(1));

    if isfile('info.rhd')
        %  Uses a modified Intan function, to make the information output more
        %  straightforward into an 'INTAN_hdr' variable.
        [sessions.info.INTAN_hdr] = mod_read_Intan_RHD2000_file('info.rhd');
        
        % Number of channels.
        sessions.info.nChannels = length(sessions.info.INTAN_hdr.amplifier_channels);
        
        % If available, parse info from 'settings.xml' (for easy access).
        %  Contains metadata that may be worth to keep. Some data is relocated
        %  to have easier access.
        if isfile('settings.xml')
            settingStruct = parseXML('settings.xml');

            % sessions.info.amplifier_sample_rate   = str2double(settingStruct.Attributes(1).Value); 
            % sessions.info.lowpass_downsample      = str2double(settingStruct.Children(2).Attributes(85).Value); 
            % sessions.info.Version                 = settingStruct.Attributes(3).Value;
            % sessions.info.Name                    = settingStruct.Name;
            % sessions.info.Children                = settingStruct.Children;
            
            % Set parameters to look up in the settings.xml file
            % as     'desired_field_name',    'path.inthe.settings',                 @function_to_convert_value
            specs = {'amplifier_sample_rate', 'Attributes(1).Value',                 @str2double % Gets amp sampl_rate
                     'lowpass_downsample',    'Children(2).Attributes(85).Value',    @str2double % Gets downsmpl factor 
                     'Version',              'Attributes(3).Value',                  [] % system version
                     'Name',                 'Name',                                 [] % system name
                     'Children',             'Children',                             [] % Set of parameters from recording
                    };
            
            for i = 1:size(specs,1)
                sessions.info.(specs{i,1}) = getField(settingStruct, specs{i,2}, specs{i,3});
            end

            % We calculate lowpass sampling rate
            sessions.info.lowpass_sample_rate     = sessions.info.amplifier_sample_rate / sessions.info.lowpass_downsample;
        else
            % Old datasets do not necessarily have an associatted .xml file
            % We have to look for basic info in the header. Some cannot be
            % populated yet.
            sessions.info.amplifier_sample_rate   = sessions.info.INTAN_hdr.frequency_parameters.amplifier_sample_rate  ;
            sessions.info.lowpass_downsample      = []; 
            sessions.info.Version                 = str2double(sessions.info.INTAN_hdr.version);
            sessions.info.Name                    = [];
            sessions.info.Children                = [];
            
            % Lowpass sampling rate not available
            sessions.info.lowpass_sample_rate     = [];
        end
    end

input.sessions(input.run(1)) = sessions;
end

function val = getField(s, pathStr, converter)
% getField  Extract a value from a nested struct using a path string.
    if nargin < 3, converter = []; end

    val = s;
    parts = strsplit(pathStr, '.');
    for k = 1:numel(parts)
        name   = regexp(parts{k}, '^\w+',              'match', 'once');
        idxStr = regexp(parts{k}, '(?<=\()\d+(?=\))',  'match', 'once');
        val = val.(name);
        if ~isempty(idxStr)
            val = val(str2double(idxStr));
        end
    end

    if ~isempty(converter)
        val = converter(val);
    end
end