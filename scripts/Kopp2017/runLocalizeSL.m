% Example script for use of LocalizeSL
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-03-25 20:55:45 -0400

rootdir='~/Documents/Github/LocalizeSL';
addpath(fullfile(rootdir,'MFILES'));

selectedSites=[24];


%%%%

for ccc=1:2

    if ccc==1
        corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore.mat'));
        ccclab='K14';
    elseif ccc==2
        corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat'));
        ccclab='DP16';
    end

    [sampsGSLrise,sampsGSLcomponents,~,~,targyears,~,colsGSL] = LocalizeStoredProjections(0,corefile);

    cols=colsGSL;
    cols.colIS=[cols.colAIS cols.colGIS];
    cols.colLI=[cols.colIS cols.colGIC];
    

    for selectedSite=selectedSites(:)'

        % generate local samples

        [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);
        nameshort=sitenames{1}(1:3);

        % separate out thermal expansion and DSL
        for nnn=1:length(sampsGSLcomponents)
            sampsloccomponents{nnn}(:,end+1,:)=sampsloccomponents{nnn}(:,cols.colOD,:)-sampsGSLcomponents{nnn}(:,colsGSL.colTE,:);
            sampsloccomponents{nnn}(:,cols.colTE,:)=sampsGSLcomponents{nnn}(:,colsGSL.colTE,:);
        end
        cols.colOD=size(sampsloccomponents{1},2);

        % plot curves

        clf;
        [hp,hax,hax2]=PlotSLRProjection(sampslocrise,targyears);
        axes(hax(1));
        title([sitenames{1} ' - ' ccclab]);
        %%pdfwrite(['LSLproj_' ccclab '_' nameshort '_' num2str(siteids(1)) ]);

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
        %pdfwrite(['LSLprojvar_' ccclab '_' nameshort '_' num2str(siteids(1))]);

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

        subcomp={cols.colAIS,cols.colGIS,cols.colGIC, cols.colTE, ...
        cols.colLS, cols.colOD, ...
        cols.colGIA,1:size(sampsloccomponents{1},2) };

        complabls={'AIS','GIS','GIC','TE','LWS','DSL','Geo','Total'};

        WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,['LSLproj_decomp_' ccclab '_' nameshort '_'],subcomp,complabls);

        clf;
        subplot(3,1,1);
        scensel=1; timesel=5;
        PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);
        subplot(3,1,2);
        scensel=1; timesel=10;
        PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
        xl=get(gca,'xlim');
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);
        subplot(3,1,3);
        scensel=4; timesel=10;
        PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
        set(gca,'xlim',xl);
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);
        %pdfwrite(['LSLproj_decomp_' ccclab '_' nameshort]);

    end
    % pull GSL samples
    [sampsGSLrise,sampsGSLcomponents,GSLsiteids,GSLsitenames,GSLtargyears,GSLscens,GSLcols] = LocalizeStoredProjections(0,corefile);
    WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,GSLtargyears,GSLcols,GSLscens,['GSLproj_decomp_' ccclab '_']);

end