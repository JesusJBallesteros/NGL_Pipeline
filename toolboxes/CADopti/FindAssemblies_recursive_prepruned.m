function [assemD]=FindAssemblies_recursive_prepruned(binMaus,w1,w2,maxlag,Dc,ref_lag)
%  © 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de

   
%% zero order  
assembly_in.n{1}.elements=w1;
assembly_in.n{1}.lag=[];
assembly_in.n{1}.pr=[];
assembly_in.n{1}.Time=binMaus(1,:)';
assembly_in.n{1}.Noccurrences=sum(binMaus(1,:));

assembly_in.n{2}.elements=w2;
assembly_in.n{2}.lag=[];
assembly_in.n{2}.pr=[];
assembly_in.n{2}.Time=binMaus(2,:)';
assembly_in.n{2}.Noccurrences=sum(binMaus(2,:));
    
%% first order: test over pairs

[assemD]=TestPair_ref(assembly_in.n{1},binMaus(2,:)',w2,maxlag,Dc,ref_lag);



end


      







