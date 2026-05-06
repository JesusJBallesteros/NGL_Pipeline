function sessions = findSetting(sessions)
% Find out INTAN settings and header file. Extracts the info.
% Using 'specs' and 'getField' as a more comprehensive parameter field, in
% case things change in the future. Still rigid in a way.
%
% Version 05.05.2026 Jesus

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
            
            sessions.info = struct();
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