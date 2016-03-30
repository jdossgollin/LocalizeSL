function [sampsDP]=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap);

% [sampsDP]=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap);
%
% Generates substitute samples (sampsDP) to be used with LocalizeSL
% by substituting Deconto & Pollard (2016) AIS ensemble members
% (with equal probability).
%
% EXAMPLE
% -------
%
% [sampsGSLrise,sampsGSLcomponents,siteidsGSL,sitenamesGSL,targyearsGSL,scensGSL,colsGSL] ...
%        = LocalizeStoredProjections(0,corefile);
% [RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap] ...
%        = DecontoPollardEnsembleImport(DecontoPollardpath,targyears);
% sampsDP=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2, ...
%         ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap);
% clear RDWAIS2 RDEAIS2;
%
% for jjj=1:size(sampsDP,1);
%     for kkk=1:size(sampsDP,2)
%         clear substituteDP; substituteDP.samps = sampsDP{jjj,kkk};
%         [sampsGSLrise2s{jjj,kkk},sampsGSLcomponents2s{jjj,kkk}]=LocalizeStoredProjections(0,corefile, ...
%                                                                 RDscenmap,substituteDP);
%     end
% end
%
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Mar 30 17:57:29 EDT 2016

samps=sampsGSLcomponents{1};
for jjj=1:size(sampsGSLcomponents,2)
    samps(:,:,:,jjj)=sampsGSLcomponents{jjj};
end

for sss=1:length(ensembleLab);
    for qqq=1:length(subsets)
        subset=subsets{qqq};
        sampDP{qqq,sss}=samps;
        
        for rrr=1:length(RDscenmap)
            
            uGSL2=samps(:,:,:,RDscenmap(rrr));

            wRDWAIS2=RDWAIS2{qqq,rrr}(:,ensembleset(:,sss));
            wRDEAIS2=RDEAIS2{qqq,rrr}(:,ensembleset(:,sss));
            
            RDWAISsamps=repmat(wRDWAIS2,1,ceil(size(uGSL2,1)/size(wRDWAIS2,2)))';
            RDWAISsamps=RDWAISsamps(1:size(uGSL2,1),:);
            uGSL2(:,colsGSL.colAIS(1),:)=RDWAISsamps;
            
            RDEAISsamps=repmat(wRDEAIS2,1,ceil(size(uGSL2,1)/size(wRDEAIS2,2)))';
            RDEAISsamps=RDEAISsamps(1:size(uGSL2,1),:);
            uGSL2(:,colsGSL.colAIS(2),:)=RDEAISsamps;
            
            sampsDP{qqq,sss}(:,:,:,RDscenmap(rrr))=uGSL2;
        end
    end
end
