% runDecontoPollard2016FullUpdate
%
% Loop through all sites and output quantiles of projections
% with Deconto & Pollard AIS components.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 14:26:36 EST 2017

workdir='workdir-170113griddedglobal';
if ~exist(workdir,'dir')
    mkdir(workdir)
end
cd(workdir);

% set up path

rootdir='~/Dropbox/Code/LocalizeSL';
IFILES=fullfile(rootdir,'IFILES');
corefile=load('~/tmp/SLRProjections170113GRIDDEDcore.mat');
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

[RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap]=DecontoPollardEnsembleImport(DecontoPollardpath,targyears);
sampsDP=DecontoPollardEnsembleGSLCompose(sampsGSLcomponents,colsGSL,RDWAIS2,RDEAIS2,ensembleLab,subsets,ensembleids,ensembleset,RDscens,RDscenmap);
clear RDWAIS2 RDEAIS2;

for jjj=1:size(sampsDP,1);
    for kkk=1:size(sampsDP,2)
        clear substituteDP; substituteDP.samps = sampsDP{jjj,kkk};
        [sampsGSLrise2s{jjj,kkk},sampsGSLcomponents2s{jjj,kkk}]=LocalizeStoredProjections(0,corefile,RDscenmap,substituteDP);
    end
end

for preferredEnsemble=1:length(ensembleLab)
    for preferredSubset=1:length(subsets)
        preferredEnsembleSetLabel=[ensembleLab{preferredEnsemble} subsets{preferredSubset}];

        % output revised GSL

        WriteTableDecomposition(sampsGSLcomponents2s{preferredSubset,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL,['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
        WriteTableSLRProjection(sampsGSLrise2s{preferredSubset,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,['GSLproj_' preferredEnsembleSetLabel '_']);

        %%% now generate local projections

        fn=['LSLproj_' preferredEnsembleSetLabel ];
        if exist([fn '.tsv'],'file');
            delete([fn '.tsv']);
        end

        WriteTableSLRProjectionAppend(sampsGSLrise2s{preferredSubset,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,fn);
        
        for ppp=1:length(selectedSites)

            selectedSite=selectedSites(ppp);
            
            % generate local samples
            substituteDP.samps = sampsDP{preferredSubset,preferredEnsemble};
            [sampslocrise2,sampsloccomponents2,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,[],substituteDP);

            % output quantiles of projections
            quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
            WriteTableSLRProjectionAppend(sampslocrise2,quantlevs,siteids,sitenames,targyears,scens,fn);

        end
    end
end
