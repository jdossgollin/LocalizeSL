% runDecontoPollard2016FullUpdate
%
% Loop through all sites and output quantiles of projections
% with Deconto & Pollard AIS components.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jan 26 17:37:49 EST 2017

workdir='workdir-170113griddedglobal';
if ~exist(workdir,'dir')
    mkdir(workdir)
end
cd(workdir);

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

preferredEnsembleSetLabel='K14';
quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
WriteTableDecomposition(sampsGSLcomponents,quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL,['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
WriteTableSLRProjection(sampsGSLrise,quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,['GSLproj_' preferredEnsembleSetLabel '_']);

% now produce standard local samples

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

% now create Deconto-Pollard samples

%[RDWAIS2,RDEAIS2,ensembleLab,bcsets,ensembleids,ensembleset,RDscens,RDscenmap,RDtargyears]=DecontoPollardEnsembleImport(DecontoPollardpath,targyears);
%save(DecontoPollardpath,'RDWAIS2','RDEAIS2','ensembleLab','bcsets','ensembleids','ensembleset','RDscens','RDscenmap','RDtargyears');

RD=load(DecontoPollardpath);
sub=find(ismember(RD.RDtargyears,targyears));
for jjj=1:size(RD.RDWAIS2,1)
    for kkk=1:size(RD.RDWAIS2,2)
        RDWAIS2{jjj,kkk}=RD.RDWAIS2{jjj,kkk}(sub,:);
        RDEAIS2{jjj,kkk}=RD.RDWAIS2{jjj,kkk}(sub,:);
    end
end
ensembleLab=RD.ensembleLab; bcsets=RD.bcsets; ensembleids=RD.ensembleids; ensembleset=RD.ensembleset; RDscens=RD.RDscens; RDscenmap=RD.RDscenmap;
clear RD;

sampsDP=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2,ensembleLab,bcsets,ensembleids,ensembleset,RDscens,RDscenmap);

for jjj=1:size(sampsDP,1);
    for kkk=1:size(sampsDP,2)
        clear substituteDP; substituteDP.samps = sampsDP{jjj,kkk};
        [sampsGSLrise2s{jjj,kkk},sampsGSLcomponents2s{jjj,kkk}]=LocalizeStoredProjections(0,corefile,RDscenmap,substituteDP);
    end
end

for preferredEnsemble=1:length(ensembleLab)
    for preferredBC=1:length(bcsets)
        preferredEnsembleSetLabel=[ensembleLab{preferredEnsemble} bcsets{preferredBC}];

        % output revised GSL

        WriteTableDecomposition(sampsGSLcomponents2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL,['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
        WriteTableSLRProjection(sampsGSLrise2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,['GSLproj_' preferredEnsembleSetLabel '_']);

        %%% now generate local projections

        fn=['LSLproj_' preferredEnsembleSetLabel ];
        if exist([fn '.tsv'],'file');
            delete([fn '.tsv']);
        end

        WriteTableSLRProjectionAppend(sampsGSLrise2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,fn);
        
        for ppp=1:length(selectedSites)

            selectedSite=selectedSites(ppp);
            
            % generate local samples
            substituteDP.samps = sampsDP{preferredBC,preferredEnsemble};
            [sampslocrise2,sampsloccomponents2,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,[],substituteDP);

            % output quantiles of projections
            quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
            WriteTableSLRProjectionAppend(sampslocrise2,quantlevs,siteids,sitenames,targyears,scens,fn);

        end
    end
end
