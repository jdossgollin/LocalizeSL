function [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(focussites,storefile,selectscens,substitutep)

% [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] =
%          LocalizeStoredProjections(focussites,storefile,[selectscens],[substitutetp])
%
% Load stored GSL MC samples and generate local MC samples
%
% INPUTS
% ------
%
% focussites: vector of PSMSL IDs of sites of interest
%             (specify 0 if you want GSL results returned in same format)
% storefile: path of file with GSL samples and other output from
%            Kopp et al. 2014, or else structure with field assignments
%            from such file
% selectscens: indices of desired scens (1 = RCP 8.5, 2 = 6.0, 3 = 4.5, 4=2.6)
%              default = [1 2 3 4]
% substitutep: structure with alternative values to substitute
%              for variables imported from storefile
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
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-05-20 13:55:57 -0400

defval('focussites',12);
defval('storefile','SLRProjections140523core');
if isstruct(storefile)
    parseFields(storefile);
    else
load(storefile,'scens','targregions','targregionnames','targyears','samps','seeds','OceanDynRegions','OceanDynYears','OceanDynMean','OceanDynStd','OceanDynN','ThermExpYears','ThermExpMean', ...
     'ThermExpStd','OceanDynTECorr','rateprojs','rateprojssd','mergeZOSZOSTOGA','fpsite','quantlevs','colGIC','colGIS','colAIS','colLS','colTE');
end

defval('selectscens',1:length(scens));

if exist('substitutep')
    parseFields(substitutep);
end

cols.colGIC=colGIC; cols.colGIS=colGIS; cols.colAIS=colAIS; cols.colLS=colLS; cols.colTE=colTE;

focussites=focussites(:)';

if focussites==0
    % if request GSL
    for qqq=1:length(selectscens)
        sampslocrise{1,qqq} = squeeze(sum(samps(:,:,:,selectscens(qqq)),2)); 
        sampsloccomponents{1,qqq}=samps(:,:,:,selectscens(qqq));
    end
    siteids=0;
    sitenames{1}='GSL';
    cols.colGIA=[];
    cols.colOD=[];
    
else
    
    [targregions2,ia,ib]=intersect(targregions,focussites);
    if length(targregions2)==0
        error(['No matching sites found!']);
    end
    
    [ib,ic]=sort(ib); ia=ia(ic); targregions2=targregions2(ic);
    targregions2names=targregionnames(ia);

    [~, ~,~,sampslocrise,sampsloccomponents,cols.colGIA,cols.colOD] =  ProjectLSL(scens(selectscens),targregions(ia),targregionnames(ia),targyears,samps(:,:,:,selectscens),seeds,targregions,OceanDynRegions,OceanDynYears,OceanDynMean(:,:,selectscens),OceanDynStd(:,:,selectscens),OceanDynN(:,:,selectscens), ThermExpYears, ThermExpMean(:,selectscens), ThermExpStd(:,selectscens), OceanDynTECorr,  rateprojs,rateprojssd,mergeZOSZOSTOGA,fpsite,quantlevs,colGIC,colGIS,colAIS,colLS,colTE,focussites);

    siteids=focussites;
    for qq=1:length(siteids)
        sub=find(targregions==siteids(qq));
        sitenames{qq}=targregionnames{sub};
    end
end
scens=scens(selectscens);

    
function parseFields(params)

    flds=fieldnames(params);
    for qqq=1:length(flds)
        if flds{qqq} 
            assignin('caller',flds{qqq},params.(flds{qqq}));
        end
    end