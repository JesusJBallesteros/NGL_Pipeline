function sessions = findSetting(sessions)
% Find out INTAN settings and header file. Extracts the info.
%
% Version 01.03.2023 Jesus

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
            sessions.info.amplifier_sample_rate   = str2double(settingStruct.Attributes(1).Value); 
            sessions.info.lowpass_downsample      = str2double(settingStruct.Children(2).Attributes(85).Value); 
            sessions.info.Version                 = settingStruct.Attributes(3).Value;
            sessions.info.Name                    = settingStruct.Name;
            sessions.info.Children                = settingStruct.Children;
            
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