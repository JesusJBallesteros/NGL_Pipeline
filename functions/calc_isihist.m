function [isihist] = calc_isihist(spike)
% TODO description

%% Get relevant info
nclust          = numel(spike.label); % number of clusters

% Set bins for histograms
bins  = [0:0.5:200]; % in miliseconds

%% Run per cluster
for cl = 1:nclust
    % ISI is diff ts(n)-ts(n-1)
    isi = [NaN diff((spike.timestamp{cl}(:)*1000)')]'; % timestamps from seconds to miliseconds
    isihist{cl}   = histcounts(isi,bins);    
end

end