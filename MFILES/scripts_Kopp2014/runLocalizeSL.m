% Example script for use of LocalizeSL
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-06-29 08:50:54 -0400

selectedSite = 12; % use PSMSL ID here to select site

% set up path

rootdir='~/Dropbox/Code/LocalizeSL';
corefile=load(fullfile(rootdir,'IFILES/SLRProjections140523core.mat'));
addpath(fullfile(rootdir,'MFILES'));

% generate local samples

[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);

% plot curves

figure;
hp1=PlotSLRProjection(sampslocrise,targyears,[],scens);
xlim([2000 2100]); ylim([0 200]);
title(sitenames{1});

% plot variance decomposition

figure;
hp2=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1);
subplot(2,2,1); title([ sitenames{1} ' - RCP 8.5']);

figure;
hp3=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1,4);
subplot(2,2,1); title([sitenames{1} ' - RCP 2.6']);

% output quantiles of projections

quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens);

% output Monte Carlo samples
WriteTableMC(sampslocrise,[],siteids,sitenames,targyears,scens);

% output Monte Carlo samples without background trend,
% to allow incorporation of alternative estimates of background trend

WriteTableMC(sampsloccomponents,1:23,siteids,sitenames,targyears,scens,'LSLProj_nobkgd_');

% output decomposition
WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens);

% pull GSL samples
[sampsGSLrise,sampsGSLcomponents,GSLsiteids,GSLsitenames,GSLtargyears,GSLscens,GSLcols] = LocalizeStoredProjections(0,corefile);
WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,GSLtargyears,GSLcols,GSLscens);
