function sessions = findSetting(sessions, ss)
% Find out INTAN settings and header file. Extracts the info.
%
% Version 01.03.2023 Jesus

if isfile('info.rhd')
    %  Uses a modified Intan function, to make the information output more
    %  straightforward into an 'INTAN_hdr' variable.
    [sessions.info(ss).INTAN_hdr] = mod_read_Intan_RHD2000_file('info.rhd');
    
    % Number of channels.
    sessions.info(ss).nchannels = length(sessions.info(ss).INTAN_hdr.amplifier_channels);
    
    % If available, parse info from 'settings.xml' (for easy access).
    %  Contains metadata that may be worth to keep. Some data is relocated
    %  to have easier access.
        if isfile('settings.xml')
            settingStruct = parseXML('settings.xml');
            sessions.info(ss).amplifier_sample_rate   = str2double(settingStruct.Attributes(1).Value); 
            sessions.info(ss).lowpass_downsample      = str2double(settingStruct.Children(2).Attributes(85).Value); 
            sessions.info(ss).Version                 = settingStruct.Attributes(3).Value;
            sessions.info(ss).Name                    = settingStruct.Name;
            sessions.info(ss).Children                = settingStruct.Children;
            
            % We calculate lowpass sampling rate
            sessions.info(ss).lowpass_sample_rate     = sessions.info(ss).amplifier_sample_rate / sessions.info(ss).lowpass_downsample;
        else
            % Old datasets do not necessarily have an associatted .xml file
            % We have to look for basic info in the header. Some cannot be
            % populated yet.
            sessions.info(ss).amplifier_sample_rate   = sessions.info(ss).INTAN_hdr.frequency_parameters.amplifier_sample_rate  ;
            sessions.info(ss).lowpass_downsample      = []; 
            sessions.info(ss).Version                 = str2double(sessions.info(ss).INTAN_hdr.version);
            sessions.info(ss).Name                    = [];
            sessions.info(ss).Children                = [];
            
            % Lowpass sampling rate not available
            sessions.info(ss).lowpass_sample_rate     = [];
        end

end