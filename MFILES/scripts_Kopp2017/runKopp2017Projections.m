% runKopp2017Projections
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-12-29 15:15:20 -0500

% set up path

rootdir='~/Dropbox/Code/LocalizeSL';
IFILES=fullfile(rootdir,'IFILES');
corefile=load(fullfile(IFILES,'SLRProjections170113GRIDDEDcore.mat'));
DecontoPollardpath=fullfile(IFILES,'DecontoPollard-AIS-160411');
addpath(fullfile(rootdir,'MFILES'));

selectedSites=corefile.targregions;

% pull GSL samples, Kopp et al. 2014 version
[sampsGSLrise,sampsGSLcomponents,siteidsGSL,sitenamesGSL,targyearsGSL,scensGSL,colsGSL] = LocalizeStoredProjections(0,corefile);
targyears=targyearsGSL;

% output K14 GSL projections

preferredEnsembleSetLabel='K14';
quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
WriteTableDecomposition(sampsGSLcomponents,quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL,['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
WriteTableSLRProjection(sampsGSLrise,quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,['GSLproj_' preferredEnsembleSetLabel '_']);

% output K14 local projections

fn='LSLproj_K14';
if exist([fn '.tsv'],'file');
    delete([fn '.tsv']);
end

WriteTableSLRProjectionAppend(sampsGSLrise,quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,fn);
for ppp=1:length(selectedSites)
    selectedSite=selectedSites(ppp);
    [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);

    % output quantiles of projections
    quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
    WriteTableSLRProjectionAppend(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,[fn]);
end

% now create DP16 samples

% if important DP16 for the first time, do this

%[RDWAIS2,RDEAIS2,ensembleLab,bcsets,ensembleids,ensembleset,RDscens,RDscenmap,RDtargyears]=DecontoPollardEnsembleImport(DecontoPollardpath,targyears);
%save(DecontoPollardpath,'RDWAIS2','RDEAIS2','ensembleLab','bcsets','ensembleids','ensembleset','RDscens','RDscenmap','RDtargyears');

% but we've already done that, so save time by using the pre-imported
% version

RD=load(DecontoPollardpath);
sub=find(ismember(RD.RDtargyears,targyears));
for jjj=1:size(RD.RDWAIS2,1)
    for kkk=1:size(RD.RDWAIS2,2)
        RDWAIS2{jjj,kkk}=RD.RDWAIS2{jjj,kkk}(sub,:);
        RDEAIS2{jjj,kkk}=RD.RDWAIS2{jjj,kkk}(sub,:);
    end
end
ensembleLab=RD.ensembleLab; bcsets=RD.bcsets; ensembleids=RD.ensembleids; ensembleset=RD.ensembleset; RDscens=RD.RDscens; RDscenmap=RD.RDscenmap;

sampsDP=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2,ensembleLab,bcsets,ensembleids,ensembleset,RDscens,RDscenmap);

runDP16ProjectionsOutput;
runDP16SubsetsGenerate;
runDP16SubsetsPlots;
runDP16GMSLPlots;
