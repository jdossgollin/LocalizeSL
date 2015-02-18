function [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(focussites,storefile)

% [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] =
%          LocalizedStoredProjections(focussites,storefile)
%
% Load stored GSL MC samples and generate local MC samples
%
% INPUTS
% ------
%
% focussites: vector of PSMSL IDs of sites of interest
% storefile: path of file with GSL samples and other output from
%            Kopp et al. 2014
%
% OUTPUTS
% -------
% 
% sampslocrise: m x n cell array, with rows corresponding to sites and
%               columns to scenarios; each cell is a p x q array, with
%               rows corresponding to samples and columns to time
%
% sampsloccomponents: m x n cell array, with rows corresponding to sites
%                     and columns to scenarios; each cell is a p x q x r
%                     array, with dimensions corresponding to samples,
%                     contributing processes, and time
%
% siteids: PSMSL ids of sites
% sitenames: names of sites
% targyears: years
% scens: scenarios
% cols: structure with column identifiers for contributing processes
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Feb 18 08:33:24 EST 2015

defval('focussites',12);
defval('storefile','SLRProjections140523core');
load(storefile,'scens','targregions','targregionnames','targyears','samps','seeds','OceanDynRegions','OceanDynYears','OceanDynMean','OceanDynStd','OceanDynN','ThermExpYears','ThermExpMean','ThermExpStd','OceanDynTECorr','rateprojs','rateprojssd','mergeZOSZOSTOGA','fpsite','quantlevs','colGIC','colGIS','colAIS','colLS','colTE');

cols.colGIC=colGIC; cols.colGIS=colGIS; cols.colAIS=colAIS; cols.colLS=colLS; cols.colTE=colTE;

focussites=focussites(:)';

[targregions2,ia,ib]=intersect(targregions,focussites);
[ib,ic]=sort(ib); ia=ia(ic); targregions2=targregions2(ic);
targregions2names=targregionnames(ia);

[quantlocrise, quantloccomponents,quantlocscalefactors,sampslocrise,sampsloccomponents,cols.colGIA,cols.colOD] =  ProjectLSL(scens,targregions(ia),targregionnames(ia),targyears,samps,seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean,OceanDynStd,OceanDynN, ThermExpYears, ThermExpMean, ThermExpStd, OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);

siteids=focussites;
for qq=1:length(siteids)
    sub=find(targregions==siteids(qq));
    sitenames{qq}=targregionnames{sub};
end
