function Deuteron2KilosortV2(opt)
% This function is a dependency of the script Deuteron2Kilosort_wrapper,
% only necessary if the recording system in use is Deuteron
% It compiles the data save in DT2 opt.myFiles in a HDF5file per channel
% It can also perform other specificities, such as filtering and the
% consideration of event codes in the final matrix compilation.
% Filtering and the elimination of experimentally irrelevant recording 
% periods only seem to be necessary when recording with Deuteron.
%
% DEPENDENCIES:
%  importEventsDeuteronWithDLL: interaction between Matlab and Deuteron's .NET
%                 app for event reading. Retrieves the event codes in 
%                 order, with corresponding time points (as minutes from midnight),
%                 sample number, and pins with detected rising edge
%
% VERSION HISTORY:
% Author:         Aylin, Lukas & Sara
% Version:        1
% Last Change:    01.03.2023 (Jesus)
 
%% Pre-define .h5 and .bin opt.myFiles
% The number of samplesof the full file depends on how many single files we have
% created.
if strcmp(opt.ext,'DT2')
    nSamplesTotal = 262144*length(opt.myFiles);
elseif strcmp(opt.ext,'DF1')
    % I haven't figure out how much space the neural data takes on each file of this format, yet.
    nSamplesTotal = Inf; 
end

if opt.h5
    % Create complete HDF5 file matching the size needs.
    h5create(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
            '/allChnMat', [opt.numChannels nSamplesTotal], ...
            'ChunkSize', [1 opt.HDF5chunkSize], ...
            'Datatype', 'int16')
end

% We always want the .bin file.
fidDataMat = fopen(fullfile(opt.FolderProcDataMat,[opt.SavFileName '.bin']), 'a'); 

%% Open each neural data file, resize data for detection with Kilosort,
% Allocate data to its respective single-channel file.
% Differentiate between old and new Deuteron Formats
if ~strcmp(opt.ext, 'DF1')
    disp('Format is FLAT. Deprecating.')
    data_mat = cell(opt.numChannels,1);

    indexPos  = 0;
    for i = 1:length(opt.myFiles)
        % Neural data points are 16 bit words
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name));
            data = fread(fid, 'uint16');
        fclose(fid);

        % Because all channels come concatenated, we need to reshape as
        % channels x samples.
        data = reshape(data', opt.numChannels, []); 

        % Convert ADC steps into microvolts, so conversion to int16 is
        % possible without loss.
        data = int16((opt.voltageResolution * (data - opt.offset)) * 1000000);
        
        % Get nSamples coming from this file. Shouls stay constant. (only
        % for i==1)?
        nSamples = size(data,2);

        for b = 1:opt.numChannels
            % compile channels -> unaltered matrix to load into kilosort
            % distribute each row of data to its respective single-channel file

            % Write channel into cell for further writing into .bin
            data_mat{b,1} = [data_mat{b,1} data(b,:)]; 

            if opt.h5
                h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
                        '/allChnMat', ...	 % dataset name
                        data(b,:), ...       % dataset: data of a channel stored in the DT2 file (already scaled)
                        [b indexPos+1], ...  
                        [1 nSamples]);
            end
        end
        indexPos = indexPos+nSamples;
    end
    clear data

    % Write bin file
    fwrite(fidDataMat, cell2mat(data_mat), 'int16');
    fclose(fidDataMat);

else % opt.ext = DF1
    %% Allocate data to its respective single-channel file
    % Open each neural data file, resize data for detection with Kilosort,
    % Allocate data to its respective single-channel file.
    disp('Format is BLOCK. NEW')

    stream      = 1;
    tmpdata     = [];

    for i = 2:length(opt.myFiles)-1      
        fid = fopen(fullfile(opt.PathRaw, opt.myFiles(i).name), 'r');
            data = Deuteron_extractData(stream, fid, opt);
        fclose(fid);

        % Convert ADC steps into microvolts, so conversion to int16 is
        % possible without loss.
        data = int16((opt.voltageResolution * (data - opt.offset)) * 1000000);

        tmpdata = [tmpdata data];
    end
    clear data fid

    % Reshape to sort as channels x samples.
    tmdata = reshape(tmpdata, opt.numChannels, []);

    % distribute each row of data to its respective single-channel file
    if opt.h5
        for b = 1:opt.numChannels
            h5write(fullfile(opt.FolderProcDataMat, [opt.SavFileName '.h5']), ...
                    '/allChnMat', ...	     % dataset name
                    tmpdata(b,:), ...  % dataset: data of a channel, scaled
                    [b b], ...              
                    [1 size(tmpdata,2)])
        end
    end

    % Write bin file
    fwrite(fidDataMat, tmpdata, 'int16');
    fclose(fidDataMat);
end  

end