% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-11-12 11:00:08 -0500

targyears=[2000 targyearsGSL];

preferredEnsemble=1; % Pliocene 5-15 m
preferredBC=2; % include bias correction

% define subsets
icesublabels={'Full','NoCliff','SlowCliff','MedCliff','FastCliff','NoHydro','StrongHydro','LowHydro','MedHydro','Hydro'};
clear icesub;
icesub{1}=find(RD.RDparams(:,1)>=0);
icesub{end+1}=find(RD.RDparams(:,3)==0);
icesub{end+1}=find(RD.RDparams(:,3)==1);
icesub{end+1}=find(RD.RDparams(:,3)==3);
icesub{end+1}=find(RD.RDparams(:,3)==5);
icesub{end+1}=find(RD.RDparams(:,2)==0);
icesub{end+1}=find(RD.RDparams(:,2)==150);
icesub{end+1}=find(RD.RDparams(:,2)==50);
icesub{end+1}=find(RD.RDparams(:,2)==100);
icesub{end+1}=find(RD.RDparams(:,2)>0);

for sss=1:length(icesub)
for jjj=1:size(RD.RDWAIS2,1)
    for kkk=1:size(RD.RDWAIS2,2)
        RDWAISsub{sss}{jjj,kkk}=RDWAIS2{jjj,kkk}(:,icesub{sss});
        RDEAISsub{sss}{jjj,kkk}=RDEAIS2{jjj,kkk}(:,icesub{sss});
    end
end
sampsDPsub{sss}=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAISsub{sss},RDEAISsub{sss},ensembleLab(1),bcsets,ensembleids(icesub{sss}),ensembleset(icesub{sss},1),RDscens,RDscenmap);
end

for sss=1:length(icesub)
    clear substituteDP; substituteDP.samps = sampsDPsub{sss}{preferredBC,preferredEnsemble};
    [sampsGSLrise2sub{sss},sampsGSLcomponents2sub{sss}]=LocalizeStoredProjections(0,corefile,RDscenmap,substituteDP);
end

for sss=1:length(icesub)
    clear substituteDP; substituteDP.samps = sampsDPsub{sss}{preferredBC,preferredEnsemble};
    [sampsGSLrise2subNOBC{sss},sampsGSLcomponents2subNOBC{sss}]=LocalizeStoredProjections(0,corefile,RDscenmap,substituteDP);
end

% table

qlevs=[.5 .17 .83 .05 .95 .01 .99 .999];
scenlabs={'RCP 8.5','RCP 4.5','RCP 2.6'};
doscens=[1 3 4];
dot=[2050 2100 2200 2300];

fid=fopen('GMSLprojections-icesubs.tsv','w');
fprintf(fid,'\t%0.1f',qlevs*100);
fprintf(fid,'\n');

for uuu=1:length(icesub)
    fprintf(fid,[icesublabels{uuu} '\n']);
    for ss=1:length(doscens)
        sss=doscens(ss);
        wset=sampsGSLrise2sub{uuu}{1,sss};
        
        fprintf(fid,[scenlabs{ss} '\n']);
        for ttt=1:length(dot)
            sub=find(targyearsGSL==dot(ttt));
            fprintf(fid,'%0.0f',dot(ttt));
            fprintf(fid,'\t%0.0f',quantile(wset(:,sub),qlevs)/10);
            fprintf(fid,'\n');
        end
    end
end


fclose(fid);
%%%

