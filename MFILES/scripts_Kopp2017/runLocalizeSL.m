% Example script for use of LocalizeSL
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-01-16 17:03:19 -0500

rootdir='~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootdir,'MFILES'));

selectedSites=[12 180];

%%%%

for ccc=1:2

    if ccc==1
        corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore.mat'));
        ccclab='K14';
    elseif ccc==2
        corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat'));
        ccclab='DP16';
    end

    for selectedSite=selectedSites(:)'

        % generate local samples

        [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);
        nameshort=sitenames{1}(1:3);

        % plot curves

        clf;
        [hp,hax,hax2]=PlotSLRProjection(sampslocrise,targyears);
        axes(hax(1));
        title([sitenames{1} ' - ' ccclab]);
        pdfwrite(['LSLproj_' ccclab '_' nameshort '_' num2str(siteids(1)) ]);

        % plot variance decomposition

        cols.colIS=[cols.colAIS cols.colGIS];
        cols.colLI=[cols.colIS cols.colGIC];
 
        subcomp={cols.colAIS,cols.colIS,cols.colLI, [cols.colLI cols.colTE], ...
        [cols.colLI cols.colTE cols.colLS], [cols.colLI cols.colTE cols.colLS cols.colOD], ...
        [cols.colLI cols.colTE cols.colLS cols.colOD cols.colGIA]};

        complabls={'AIS','GIS','GIC','TE','LWS','DSL','Geo'};

        clf;
        [hp,vars,fvars,hlg]=PlotSLRProjectionVariance(sampsloccomponents(:,[1 3 4]),targyears,cols,[2010 2100],1,2,1,subcomp,complabls,'rcbgmykrcbgm');
        title(sitenames{1});
        pdfwrite(['LSLprojvar_' ccclab '_' nameshort '_' num2str(siteids(1))]);

        % output quantiles of projections

        quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
        WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,['LSLproj_' ccclab '_' nameshort '_']);

        % output timing of height exceedances
        WriteTableSLRHeightExceedanceTiming(sampslocrise,[],siteids,sitenames,targyears,scens,1,['LSLheights_' ccclab '_' nameshort '_']);

        % output Monte Carlo samples

        WriteTableMC(sampsloccomponents,[],siteids,sitenames,targyears,scens,['LSLproj_MC_' ccclab '_' nameshort '_']);

        % output Monte Carlo samples without background trend,
        % to allow incorporation of alternative estimates of background trend

        WriteTableMC(sampsloccomponents,setdiff(1:size(sampsloccomponents{1},2),cols.colGIA),siteids,sitenames,targyears,scens,['LSLproj_MC_nobkgd_' ccclab '_' nameshort '_']);

        % output decomposition
        WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,['LSLproj_decomp_' ccclab '_' nameshort '_']);


    end
    % pull GSL samples
    [sampsGSLrise,sampsGSLcomponents,GSLsiteids,GSLsitenames,GSLtargyears,GSLscens,GSLcols] = LocalizeStoredProjections(0,corefile);
    WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,GSLtargyears,GSLcols,GSLscens,['GSLproj_decomp_' ccclab '_']);

end