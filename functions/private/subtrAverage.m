function subtrAverage(meanForSub,numChannels,sngChn,folderProcDataMatAveraged,savFileNameAvrg,j,stpSz,fidDataMatAvg,ChunkStart)
% write into a new matrix and a new binary file the average subtracted channel data  

%average of all channels 
meanForSub = int16(sum(meanForSub,1));

sngChnAvg = cell(numChannels,1);

for k=1:numChannels

    %average subtracted
    if iscell(sngChn)
        sngChnAvg{k,1} = sngChn{k,1}-meanForSub; 
    else
        sngChnAvg{k,1} = sngChn(k,1)-meanForSub;
    end 

    % write average subtracted data into a new matrix 
    if nargin == 9
        h5write(fullfile(folderProcDataMatAveraged, [savFileNameAvrg '.h5']), '/avgSubtracted', ...
            sngChnAvg{k,1}, [k j-(ChunkStart(1)-1)], [1 stpSz]);
    else
         h5write(fullfile(folderProcDataMatAveraged, [savFileNameAvrg '.h5']), '/avgSubtracted', ...
            sngChnAvg{k,1}, [k j], [1 stpSz]);
    end

end

% write average subtracted data into a new binary file
fwrite(fidDataMatAvg,cell2mat(sngChnAvg),'int16'); 

end