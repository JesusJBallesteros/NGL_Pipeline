function [As_across_bins,As_across_bins_index,assembly,Assemblies_all_orders]=CADopti(spM,MaxLags,BinSizes,ref_lag,alph,No_th,O_th,bytelimit)
% this function returns cell assemblies detected in spM spike matrix binned 
% at a temporal resolution specified in 'BinSizes' vector and testing for all 
% lags between '-MaxLags(i)' and 'MaxLags(i)'
%
% USAGE: [assembly]=Main_assemblies_detection(spM, MaxLags, BinSizes, ref_lag, alph, Dc, No_th, O_th, bytelimit)
%
% ARGUMENTS:
% spM     := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
% BinSizes:= vector of bin sizes to be tested;
% MaxLags:= vector of maximal lags to be tested. For a binning dimension of BinSizes(i) the program will test all pairs configurations with a time shift between -MaxLags(i) and MaxLags(i);
% (optional) ref_lag      := reference lag. Default value 2
% (optional) alph      := alpha level. Default value 0.05
% (optional) No_th     := minimal number of occurrences required for an assembly (all assemblies, even if significant, with fewer occurrences than No_th are discarded). Default value 0.
% (optional) O_th      := maximal assembly order (the algorithm will return assemblies of composed by maximum O_th elements).
% (optional) bytelimit := maximal size (in bytes) allocated for all assembly structures detected with a bin dimension. When the size limit is reached the algorithm stops adding new units. 
%
% RETURNS:
% assembly - structure containing assembly information:
%     assembly.parameters       - parameters used to run Main_assemblies_detection
%     assembly.bin{i} contains information about assemblies detected with 
%                     'BinSizes(i)' bin size tested for all lags between 
%                     '-MaxLags(i)' and 'MaxLags(i)'
% 
%        assembly.bin{i}.bin_edges - bin edges (common to all assemblies in assembly.bin{i})
%        assembly.bin{i}.n{j} information about the j-th assembly detected with BinSizes(i) bin size 
%                   elements: vector of units taking part to the assembly (unit order correspond to the agglomeration order)
%                        lag: vector of time lags. '.lag(z)' is the activation delay between .elements(1) and .elements(z+1)
%                         pr: vector of pvalues. '.pr(z)' is the pvalue of the statistical test between performed adding .elements(z+1) to the structure .elements(1:z)
%                       Time: assembly activation time. It reports how many times the complete assembly activates in that bin. .Time always refers to the activation of the first listed assembly element (.elements(1)), that doesn't necessarily corresponds to the first unit firing.
%               Noccurrences: number of assembly occurrence. '.Noccurrences(z)' is the occurrence number of the structure composed by the units .elements(1:z+1) 
%
% As_across_bins         - structure containing assembly information (exactly same information contained in "assembly" but collected across different temporal resolutions)
% As_across_bins_index   - information to link assemblies in "As_across_bins"
%                          back to the structure "assembly":
%                          assembly As_across_bins{i} is contained in assembly.bin{As_across_bins_index{i}(1)}.n{As_across_bins_index{i}(2)}.
%
%  © 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de


if nargin<4 || isempty(ref_lag), ref_lag=2; end  
if nargin<5 || isempty(alph), alph=0.05; end  
if nargin<6 || isempty(No_th), No_th=0; end      % no limitation on the number of assembly occurrences
if nargin<7 || isempty(O_th), O_th=Inf; end     % no limitation on the assembly order (=number of elements in the assembly)
if nargin<8 || isempty(bytelimit), bytelimit=Inf; end     % no limitation on assembly dimension

%%

nneu=size(spM,1); % number of units
testit=ones(1,length(BinSizes));
binM=cell(1,length(BinSizes));  
number_tests=0;
% matrix binning at all bins
for gg=1:length(BinSizes)
    int=BinSizes(gg);
    tb=min(spM(:)):int:max(spM(:));  
    
    binM{gg}=zeros(nneu,length(tb)-1,'uint8');
    number_tests=number_tests+nneu*(nneu-1)*(2*MaxLags(gg)+1)/2;
     for n=1:nneu
        [ binM{gg}(n,:),~] = histcounts(spM(n,:),tb);
     end 
     assembly.bin{gg}.n=[];
     assembly.bin{gg}.bin_edges=tb;
     if size(binM{gg},2)-MaxLags(gg)<100
        fprintf('Warning: testing bin size=%f. The time series is too short, consider taking a longer portion of spike train or diminish the bin size to be tested \n', int);
        testit(gg)=0;
     end
end

fprintf('order 1\n')
clear Assemblies_all_orders
O=1;
Dc=100; %length (in # bins) of the segments in which the spike train is divided to compute #abba variance (parameter k).


assembly_selected_xy=[];
% assembly_selected_xy=nan(nneu,nneu);
p_values=[];
% first order assembly
parfor w1=1:nneu   
    for w2=w1+1:nneu
        assemblybin=cell(1,length(BinSizes));  
        p_by_bin=[];
        for gg=1:length(BinSizes)              
            [assemblybin{gg}]=FindAssemblies_recursive_prepruned([binM{gg}(w1,:);binM{gg}(w2,:)],w1,w2,MaxLags(gg),Dc,ref_lag);
            p_values=[p_values; assemblybin{gg}.pr(end)];
            assemblybin{gg}.bin=BinSizes(gg);
            p_by_bin(gg)=assemblybin{gg}.pr;   
        end
        [~, b]=min(p_by_bin);        
        assembly_selected_xy=[assembly_selected_xy,assemblybin(b)];
    end
end
assembly_selected=assembly_selected_xy;

if ~isempty(assembly_selected)

%% Holm-Bonferroni 
x=1:length(p_values);
p_values=sort(p_values);
p_values_alpha=alph./(number_tests+1-x);
aus=find((p_values'-p_values_alpha)<0);
if isempty(aus)
    HBcorrected_p=0;
else
    HBcorrected_p=p_values(aus(end));   
end
ANfo=zeros(nneu,nneu);

for oo=length(assembly_selected):-1:1
    if assembly_selected{oo}.pr(end)>HBcorrected_p
         assembly_selected(oo)=[];
    else
         ANfo(assembly_selected{oo}.elements(1),assembly_selected{oo}.elements(2))=1; 
    end
end
Assemblies_all_orders{O}=assembly_selected;
      
%%
% higher orders
Oincrement=1;
while Oincrement && O<(O_th-1)
    O=O+1;
    fprintf('order %d\n',O)
    Oincrement=0;
    assembly_selected_aus=[];
    xx=1;
    for w1=1:size(assembly_selected,2)

        % bin at which to test w1
        ggg=find(BinSizes==assembly_selected{w1}.bin);

        % element to test with w1
        w1_elements=assembly_selected{w1}.elements;
        [~, w2_to_test]=find(ANfo(w1_elements,:)==1);   % I try to add only neurons that have significant first order cooccurrences with members of the assembly
        w2_to_test(ismember(w2_to_test,w1_elements))=[];  % I erase the one that are already in the assembly
        w2_to_test=unique(w2_to_test);

        for ww2=1:length(w2_to_test)
                w2=w2_to_test(ww2);
                spikeTrain2=binM{ggg}(w2,:)';
                [assemblybin_aus]=TestPair_ref(assembly_selected{w1},spikeTrain2,w2,MaxLags(ggg),Dc,ref_lag);
                p_values=[p_values; assemblybin_aus.pr(end)];
                number_tests=number_tests+2*MaxLags(ggg)+1;
                if assemblybin_aus.pr(end)<HBcorrected_p
                     assembly_selected_aus{xx}=assemblybin_aus;
                     assembly_selected_aus{xx}.bin=BinSizes(ggg);
                     xx=xx+1;
                     Oincrement=1;
                end
        end
    end

    if Oincrement
        %%% pruning within the same size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % between two assemblies with the same unit set arranged into different configurations I choose the most significant one

        na=length(assembly_selected_aus);                      % number assemblies
        nelement=size(assembly_selected_aus{1}.elements,2);    % number elements for assembly
        selection=nan(na,nelement+1+1);
        assembly_final=cell(1,na);     %max possible dimension
        nns=1;
        for i=1:na
            elem=sort(assembly_selected_aus{i}.elements);
            [ism,indx]=ismember(elem,selection(:,1:nelement),'rows');
            if ~ism
                assembly_final{nns}=assembly_selected_aus{i};      
                selection(nns,1:nelement)=elem;
                selection(nns,nelement+1)=assembly_selected_aus{i}.pr(end);
                selection(nns,nelement+2)=i;
                nns=nns+1;
            else
                if selection(indx,nelement+1)>assembly_selected_aus{i}.pr(end)
                    assembly_final{indx}=assembly_selected_aus{i};      
                    selection(indx,nelement+1)=assembly_selected_aus{i}.pr(end);
                    selection(indx,nelement+2)=i;
                end
            end
        end
        assembly_final(nns:end)=[];   
        assembly_selected=assembly_final;
        Assemblies_all_orders{O}=assembly_final;
        clear assembly_final
    end
    
    %% Holm-Bonferroni 
    x=1:length(p_values);
    p_values=sort(p_values);
    p_values_alpha=alph./(number_tests+1-x);
    aus=find((p_values'-p_values_alpha)<0);
%     HBcorrected_p=p_values(aus(end));
    if isempty(aus)
        HBcorrected_p=0;
    else
        HBcorrected_p=p_values(aus(end));   
    end


end




%% Holmâ€“Bonferroni 
x=1:length(p_values);
p_values=sort(p_values)';
p_values_alpha=alph./(number_tests+1-x);
aus=find((p_values-p_values_alpha)<0);
% HBcorrected_p=p_values(aus(end));
if isempty(aus)
    HBcorrected_p=0;
else
    HBcorrected_p=p_values(aus(end));   
end

for o=1:length(Assemblies_all_orders)
    for oo=length(Assemblies_all_orders{o}):-1:1
        if Assemblies_all_orders{o}{oo}.pr(end)>HBcorrected_p
             Assemblies_all_orders{o}(oo)=[];
        end
    end
end

%% pruning between differen assembly size
o=length(Assemblies_all_orders);
x=1;
for oo=length(Assemblies_all_orders{o}):-1:1
    Element_template{x}=Assemblies_all_orders{o}{oo}.elements;
    x=x+1;
end


for o=length(Assemblies_all_orders)-1:-1:1
    for oo=length(Assemblies_all_orders{o}):-1:1
        found=0;
        ooo=1;
        while ~found && ooo<x
            if ismember(Assemblies_all_orders{o}{oo}.elements,Element_template{ooo})
                Assemblies_all_orders{o}(oo)=[];
                found=1;
            else
               ooo=ooo+1;               
            end
        end
        if found==0
            Element_template{x}=Assemblies_all_orders{o}{oo}.elements;
            x=x+1;
        end
    end
end

%% reformat dividing by bins
for o=length(Assemblies_all_orders):-1:1
    for oo=length(Assemblies_all_orders{o}):-1:1
        bx=find(BinSizes==Assemblies_all_orders{o}{oo}.bin);
        assembly.bin{bx}.n=[assembly.bin{bx}.n,Assemblies_all_orders{o}(oo)];
    end
end
for gg=length(BinSizes):-1:1
    if isempty(assembly.bin{gg}.n)
        assembly.bin{gg}=[];
    end
end

fprintf('\n');
else
    for gg=length(BinSizes):-1:1

        assembly.bin{gg}=[];
Assemblies_all_orders=[];
end

end

assembly.parameters.alph=alph;
assembly.parameters.Dc=Dc;
assembly.parameters.No_th=No_th;
assembly.parameters.O_th=O_th;
assembly.parameters.bytelimit=bytelimit;
assembly.parameters.ref_lag=ref_lag;


[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);


end    


