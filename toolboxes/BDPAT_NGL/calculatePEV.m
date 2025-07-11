function effectSize = calculatePEV(tbl,value)
%%
%effectSize = calculatePEV(tbl,value)
%
% This function calculates the percent explained variance of a factor in a
% one-way or two-way ANOVA (omega2 and partialomega2 respectively).
% Calculation of eta2 is also possible based on an ANOVA result.
%
%INPUTS
%  * 'tbl'           : ANOVA table (-> anovan)
%  * 'value'         : string of effect size value (omega2, partialomega2, 
%                      eta2)
%
%OUTPUTS
%  * 'effectSize'           : effect size

% VERSION HISTORY:
% Author:         Lukas Hahn
% Version:        1.0.1
% Last Change:    12.12.2023
%
% 15.07.2019, Lukas: v1.0.0 release version
% 12.12.2023, Lukas: v1.0.1 updated documentation
%% calculate PEV
% tbl - ANOVA table (works for one way and omega2 and for two way and
%       partialomega2
%       feeds
%           k - degrees of freedom
%           SSQ - sum of squares of tested condition
%           MSerr - mean squares error
%           SSQtotal - sum of squares total
% value - string of effect size measure ('omega2' or 'partialomega2')
% for more information visit 
% https://gitlab.ruhr-uni-bochum.de/ngl/
% basicdataprocessingandanalysistools/-/wikis/calculatePEV
%%
if strcmp(value,'omega2')
    if ~all(size(tbl)==[4 7])
        error(['Check dimensions of ANOVA table!' newline ...
            'This function exclusively accepts ANOVA tables from the '...
            'function ''anovan''!'])
    else
    end
    SSQeff = tbl{2,2}; %sum of squares of factor 1 ('Sum Sq.')
    SSQtotal = tbl{4,2}; %total sum of squares ('Sum Sq.')
    keff = tbl{2,3}; %degrees of freedom of factor 1 ('d.f.')
    MSerr = tbl{3,5}; %mean error squares ('Mean Sq.')
    omega2 = (SSQeff-(keff*MSerr))/(SSQtotal+MSerr);
    effectSize = omega2;
elseif strcmp(value,'partialomega2')
    if ~all(size(tbl)==[6 7])
        error(['Check dimensions of ANOVA table!' newline ...
            'This function exclusively accepts ANOVA tables from the '...
            'function ''anovan''!'])
    else
    end
    partOmega2 = nan(3,1);
    for i=1:3
        keff = tbl{i+1,3}; %keff = level-1;
        Feff = tbl{i+1,6}; %F = MSeff/MSerr;
        N = tbl{6,3}+1;
        partOmega2(i,1) = (keff*(Feff-1))/(keff*(Feff-1)+N);
    end
    effectSize = partOmega2;
elseif strcmp(value,'eta2')
    if ~all(size(tbl)==[4 7])
        error(['Check dimensions of ANOVA table!' newline ...
            'This function exclusively accepts ANOVA tables from the '...
            'function ''anovan''!'])
    else
    end
    SSQeff = tbl{2,2}; %sum of squares of factor 1 ('Sum Sq.')
    SSQtotal = tbl{4,2}; %total sum of squares ('Sum Sq.')
    eta2 = SSQeff/SSQtotal;
    effectSize = eta2;
end
end