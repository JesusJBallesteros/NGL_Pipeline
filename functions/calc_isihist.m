function [isihist] = calc_isihist(spike)
% TODO description

%% Get relevant info
nclust          = numel(spike.label); % number of clusters

% Set bins for histograms
bins  = [0:0.5:200];

%% Run per cluster
for cl = 1:nclust
    % ISI is diff ts(n)-ts(n-1)
    isi = [NaN diff(spike.timestamp{cl}(:)')]';
    isihist{cl}   = histcounts(isi,bins);    
end

%% TODO Give histogram as proportion
% isihist = isihist./repmat(sum(isihist,2,'omitnan'),1,size(isihist,2));

end