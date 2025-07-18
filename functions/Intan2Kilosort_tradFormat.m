function Intan2Kilosort_tradFormat(opt)
% Function created in case that data is collected in traditional format
% from Intan, by mistake. NOT RECOMMENDED.
% Version 14.07.2025 Jesus

[INTAN_hdr] = mod_read_Intan_RHD2000_file(fullfile(opt.myFiles(1).folder, opt.myFiles(1).name),1);

%% Get data from INTAN data files, as single channels.
% Opens neural data file and coverts the ADC steps to microvolts, then
% write channel-by-channel into .h5 file and the whole array into a .bin
% file. Filters the data if necessary/requested.  
data = INTAN_hdr.amplifier_data;

% Get Sample rate, from Intan_hdr.
opt.sampleRate  = INTAN_hdr.frequency_parameters.amplifier_sample_rate;
clear INTAN_hdr

% Make noisy channels NaNs
if ~isempty(opt.noise)
    data(opt.noise,:) = nan(size(data(opt.noise,:)));
end

% In principle, only for data from a single active zone...
if opt.CAR
    disp('Re-referencing by Common Average Referencing (CARing).')
    data = ft_preproc_rereference(data, 'all', 'median', true);
end


% Make noisy channels zeros
if ~isempty(opt.noise)
    data(opt.noise,:) = zeros(size(data(opt.noise,:)));
end

txt = sprintf('Highpass filter set at %d Hz. It may take a moment.\n', opt.highpass);
fprintf(txt);

% Keep memory usage low doing one channel at a time.
for i = 1:opt.numChannels
    fprintf('- Filtering channel %d of %d.\n', i, opt.numChannels);
    
    % Detrend channel (remove DC)
    disp('Detrending...')
    data(i,:) = ft_preproc_detrend(data(i,:));

    % Highpass channel (Butterwort, 6th order, back&forth)
    disp('Filtering...')
    [data(i,:), ~, ~] = ft_preproc_highpassfilter(data(i,:), opt.sampleRate, opt.highpass, 6, 'but', 'twopass');
end

% Create or open a bin file. Append data at end.
opt.binfilename = fullfile(opt.FolderProcDataMat,[opt.SavFileName + ".bin"]);
if isfile(opt.binfilename)
    delete(opt.binfilename);
end

fidDataMat = fopen(opt.binfilename, 'a');

fwrite(fidDataMat, data, 'int16');
fclose(fidDataMat);

end