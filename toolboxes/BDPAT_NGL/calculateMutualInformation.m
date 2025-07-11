function MI = calculateMutualInformation(pXY)
%% 
%function MI = calculateMutualInformation(pXY)
%
%This function calculates the mutual information (MI) of values in a
%matrix. Theinput matrix (pXY) needs to have an applicable data-format, 
% e.g. a classifier confusion matrix, which contains only values of joint
% probability, from which marginal probabilites are automatically
% calculated in this function. The function uses a base 2 logarithm to
% output MI in bits (see also Shannon-Hartley theorem). Range of MI is
% between 0 (no mutual information, i.e. no relation between the variables
% X and Y) and an upper limit based on the number of different groups g
% (e.g., classes of the classifier), with maxBit = log2((1/g)/((1/g)^2))
% e.g. with g = 3, the result falls into the interval [0 1.585]
% based on the matrix
% pXY = [0.1667 0.0417 0.125; 0.125 0.25 0.0833; 0.0417 0.0417 0.125];
% the result is 0.1727 bits
%
%INPUTS
%  * 'pXY'          : matrix containing all joint probabilities
%
%OUTPUTS
%   * 'MI'          : Mutual information in bits

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.0
% Last Change:    13.12.2023
%
% 13.12.2023, Lukas: v1.0.0 release version
%%
bitPerEl = nan(numel(pXY),1);
pXY(size(pXY,1)+1,:) = sum(pXY,1);
pXY(:,size(pXY,2)+1) = sum(pXY,2);

i = 0;
for r=1:size(pXY,1)-1
    for c=1:size(pXY,2)-1
        i = i+1;
        if ~pXY(r,c)==0
            bitPerEl(i,1) = pXY(r,c)*log2(pXY(r,c)/(pXY(r,4)*pXY(4,c)));
        else %in case joint probability is zero write 0
            %required due to 0*log2(0/x) = NaN
            bitPerEl(i,1) = 0;
        end
    end
end
MI = sum(bitPerEl);
end