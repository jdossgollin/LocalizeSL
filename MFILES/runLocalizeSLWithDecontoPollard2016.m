% runLocalizeSLWithDecontoPollard2016
%
% Generate GSL and LSL projections using Deconto and Pollard (2016)'s AIS projections
% with other projections from Kopp et al. (2014). Assumes each of the ensemble
% members that meets Deconto-Pollard criteria is of equal probability.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Mar 30 17:58:32 EDT 2016

selectedSites = [235];

% set up path

addpath(pwd);
rootdir='~/Dropbox/Code/LocalizeSL';
IFILES=fullfile(rootdir,'IFILES');
corefile=fullfile(IFILES,'SLRProjections140523core.mat');
DecontoPollardpath=fullfile(IFILES,'DecontoPollard-AIS-160125');
addpath(fullfile(rootdir,'MFILES'));

% pull GSL samples, Kopp et al. 2014 version
[sampsGSLrise,sampsGSLcomponents,siteidsGSL,sitenamesGSL,targyearsGSL,scensGSL,colsGSL] = LocalizeStoredProjections(0,corefile);
targyears=targyearsGSL;

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

        quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
        WriteTableDecomposition(sampsGSLcomponents2s{preferredEnsemble,preferredSubset},quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL(RDscenmap),['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
        WriteTableSLRProjection(sampsGSLrise2s{preferredEnsemble,preferredSubset},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL(RDscenmap),['GSLproj_' preferredEnsembleSetLabel '_']);

%%% now generate local projections

        for ppp=1:length(selectedSites)

            selectedSite=selectedSites(ppp);
            % generate local samples
            clear substituteDP; substituteDP.samps = sampsDP{preferredEnsemble,preferredSubset};
            [sampslocrise2,sampsloccomponents2,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,[],substituteDP);

            % plot curves

            figure;
            hp1=PlotSLRProjection(sampslocrise2,targyears,[],scens);
            xlim([2000 2100]); ylim([0 200]);
            title(sitenames{1});

            % plot variance decomposition

            figure;
            hp2=PlotSLRProjectionVariance(sampsloccomponents2,targyears,cols,[],1);
            subplot(2,2,1); title([ sitenames{1} ' - RCP 8.5']);

            figure;
            hp3=PlotSLRProjectionVariance(sampsloccomponents2,targyears,cols,[],1,4);
            subplot(2,2,1); title([sitenames{1} ' - RCP 2.6']);

            % output quantiles of projections

            quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
            WriteTableSLRProjection(sampslocrise2,quantlevs,siteids,sitenames,targyears,scens,['LSLproj_' preferredEnsembleSetLabel '_']);

            % output decomposition
            WriteTableDecomposition(sampsloccomponents2,quantlevs,siteids,sitenames,targyears,cols,scens,['LSLproj_decomp_' preferredEnsembleSetLabel '_']);

        end
    end
end
