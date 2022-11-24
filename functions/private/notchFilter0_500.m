function notchFilter0_500(numChannels,filePathSingleChannels,filePathSingleChannelsFilt,HDF5chunkSize,bndWLowConfirmFilter)
fileNames = ls(filePathSingleChannels);
fileNames(1:2,:) = [];
for i = 1:numChannels 
    
    %estimate of total duration 0.33 s (per channel, using filtfilt) for dataset 220304
    filteredData = filter(bndWLowConfirmFilter,...
        double(... 
        h5read(...
        fullfile(filePathSingleChannels,fileNames(i,:)),['/channel_' num2str(i)]))); %data needs 
    % to be converted from int16 to double in order to apply this filter 
    
    h5create(fullfile(filePathSingleChannelsFilt,fileNames(i,:)),...
        ['/channel_' num2str(i)],[1 Inf],'ChunkSize',[1 HDF5chunkSize],'Datatype','int16')%create a new file per filtered channel
    h5write(fullfile(filePathSingleChannelsFilt,fileNames(i,:)),...
        ['/channel_' num2str(i)],int16(filteredData),[1 1],[1 size(filteredData,2)]) %save filtered data 
end
end