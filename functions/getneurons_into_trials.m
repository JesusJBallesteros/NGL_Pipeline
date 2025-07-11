function getneurons(results, trials)
% Reads the results extracted from KS+Phy and obtains the spike times
% corresponding to valid clusted IDs. It sorts the spikes per cluster into
% its corresponding trials, if existing.

st   = results.spike.st;
clu  = results.spike.clu;
cids = results.spike.cids;

% Find out how many 'neurons'. For KS, clusters with valid ID after
% curation.
nneurons = length(cids);

% Create variable
neurons = cell(nneurons,1);

% Fill variable
for i = 1:nneurons
    neurons{i} = cell(trials,1); % trials = number of trials per session 
    
    for n = 1:trials
        % column vector with all spiketimes per trial
        neurons{i,1}{n,1} = st(clu==cids(i)); % TODO. Needs to index for time as well
    end
end

save(fullfile(opt.FolderProcDataMat, 'neurons.mat'), "neurons", '-mat');

end