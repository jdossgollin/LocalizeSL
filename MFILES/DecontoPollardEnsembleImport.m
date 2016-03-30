function [RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap]=DecontoPollardEnsembleImport(DecontoPollardpath,targyears)

% [RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap]=DecontoPollardEnsembleImport(DecontoPollardpath,targyears)
%
% Imports Deconto and Pollard (2016) Antarctic ice sheet ensemble, to be passed to 
% DecontoPollardEnsembleGSLCompose.
%
% INPUTS
% ------
% DecontoPollardpath: path to directory with Deconto & Pollard input files
% targyears: target years, used by LocalizeSL
%
% See help for DecontoPollardEnsembleGSLCompose for example.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Mar 30 17:57:23 EDT 2016

datMembers=importdata(fullfile(DecontoPollardpath,'EnsembleDefinition.csv'));
ensembleids=datMembers.data(:,1);
ensembleset=logical(datMembers.data(:,2:end));
ensembleLab=datMembers.textdata(2:end);

RDscens={'RCP85','RCP45','RCP45','RCP26'};
RDscenmap=[1 2 3 4];
subsets={'','PIT'};
for uuu=1:length(subsets)
    subset=subsets{uuu};
    for rrr=1:length(RDscenmap)
        wpath=fullfile(DecontoPollardpath,[RDscens{rrr} subset]);
        files=dir(fullfile(wpath,'*.22'));
        ii=1;
        dat=importdata(fullfile(wpath,files(ii).name));
        RDt=1950+dat(:,1);

        clear RDWAIS RDEAIS;
        for ii=1:length(ensembleids)
            dat=importdata(fullfile(wpath,[num2str(ensembleids(ii)) '.22']));
            RDWAIS(:,ii)=dat(:,end-1)*1000;
            RDEAIS(:,ii)=dat(:,end-2)*1000;
            
        end
        RDrefyr=find(RDt==2000);
        RDWAIS=bsxfun(@minus,RDWAIS,RDWAIS(RDrefyr,:));
        RDEAIS=bsxfun(@minus,RDEAIS,RDEAIS(RDrefyr,:));
        RDWAIS2{uuu,rrr}=interp1(RDt,RDWAIS,targyears);
        RDEAIS2{uuu,rrr}=interp1(RDt,RDEAIS,targyears);
    end
end

