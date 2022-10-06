function convert2nwb()
% WHAT IT IS: Substitute input layer to avoid using Jupyter Notebook GUI.
% It bypasses 'IntanToNWB.ipynb' and 'ConverterUI.py' scripts, so they are
% not longer necessary.
%
% WHAT IT DOES: It navigates to the Python folder, runs the routine, copies the .nwb file
% back to the original data folder and navigates back to it.
%
% TODO: Add inputs from Matlab and pass them to the script for easyness
%
% WHAT IT CONTAINS: 'my_IntanToNWB.py' script with default settings as:
%   intan_filename="info.rhd"       % Default header filename when multiple files are created with INTAN
%   nwb_filename="info.nwb"         % Default new file name. Suggested to be changed afterwards.
%   session_description=None        % Do not add descriptions, since the INTAN python wrapper seems to not really like it.
%   blocks_per_chunk=2000           % This can be modified to change how much data is proccesed per bulk.
%   use_compression=False           % No compression
%   compression_level=4             % Left default, not used if above False
%   lowpass_description=None        % Not defined
%   highpass_description=None       % Not defined
%   merge_files=None                % No merging header files
%   subject=None                    % No metadata info at subject level
%   manual_start_time=datetime.now()% Use current time
%
% The Python Script is reproduced here. Copy it in a new file and add it to
% the folder where the IntanToNWB scripts are as '.py' file. Modify the
% name in pyrunfile("new_name.py"). 
% 'my_IntanToNWB.py' file starts below (remove % % symbols when pasting):
% % from ConvertIntanToNWB import *
% % from datetime import datetime
% % """ Callback function for when the user begins conversion.
% % 
% % Parameters
% % ----------
% % None
% % 
% % Returns
% % -------
% % None
% % """
% %         
% % # Get NWB settings and pass these to NWB conversion
% % intan_filename="info.rhd"
% % nwb_filename="info.nwb"
% % session_description=None
% % blocks_per_chunk=2000
% % use_compression=False
% % compression_level=4
% % lowpass_description=None
% % highpass_description=None
% % merge_files=None
% % subject=None
% % manual_start_time=datetime.now()
% % 
% % convert_to_nwb(intan_filename=intan_filename,
% %                nwb_filename=nwb_filename,
% %                session_description=session_description,
% %                blocks_per_chunk=blocks_per_chunk,
% %                use_compression=use_compression,
% %                compression_level=compression_level,
% %                lowpass_description=lowpass_description,
% %                highpass_description=highpass_description,
% %                merge_files=merge_files,
% %                subject=subject,
% %                manual_start_time=manual_start_time)
% File ends above.
%
% By Jesus J. Ballesteros 09.2022

%% 000. Verify Python Configuration
pe = pyenv

%% Run the routine
disp('- Conversion in progress...');
pyrunfile("my_IntanToNWB.py");

end