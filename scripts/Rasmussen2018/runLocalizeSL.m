% Example script for use of LocalizeSL
% with Rasmussen et al. (2018) projections.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-04-01 19:57:04 -0400

selectedSite = 12; % use PSMSL ID here to select site

% set up path

rootdir='~/Dropbox/Code/LocalizeSL';
corefile=load(fullfile(rootdir,'IFILES/SLRProjections180124GRIDDEDcore_Tscens.mat'));
addpath(fullfile(rootdir,'MFILES'));

% specify scenario labels and scenarios to use
% important since we differ from defaults here

scenlabs={'tmp15','tmp20','tmp25'};
selscens=[1 2 3];

% generate local samples

[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,selscens);

% plot curves

figure;
hp1=PlotSLRProjection(sampslocrise,targyears,[],scenlabs,selscens);
title(sitenames{1});

% plot variance decomposition

figure;
hp2=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1,1);
subplot(2,2,1); title([ sitenames{1} ' - 1.5 C']);

figure;
hp3=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1,2);
subplot(2,2,1); title([sitenames{1} ' - 2.0 C']);

% output quantiles of projections

quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens);

% output timing of height exceedances
WriteTableSLRHeightExceedanceTiming(sampslocrise,[],siteids,sitenames,targyears,scens,1);

% output Monte Carlo samples
WriteTableMC(sampslocrise,[],siteids,sitenames,targyears,scens);

% output Monte Carlo samples without background trend,
% to allow incorporation of alternative estimates of background trend

WriteTableMC(sampsloccomponents,1:23,siteids,sitenames,targyears,scens,'LSLProj_nobkgd_');

% output decomposition
WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens);

