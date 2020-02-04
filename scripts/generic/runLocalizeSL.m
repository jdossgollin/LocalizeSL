% Example script for use of LocalizeSL to generate common tables and plots
% for RSL projections.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-07-25 20:42:48 -0400

selectedSites=[12 180]; % specify the PSMSL or grid cell IDs of interest

rootdir='~/Dropbox/Code/LocalizeSL'; % change to directory containing LocalizeSL
corefiles={'SLRProjections170113GRIDDEDcore.mat','SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat','SLRProjections180124GRIDDEDcore_Tscens.mat','SLRProjections190726core_SEJ_full.mat','SLRProjections190726core_SEJ_full.mat','SLRProjections200204GRIDDEDcore_SROCC.mat'}; % specify corefiles to use
corefilelabs={'K14','DP16','R18','B19H','B19L','K14SROCC'}; % specify corefile labels
subcore={'','','','corefileH','corefileL',''}; % specify if corefile file contains multiple cores
quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999]; % quantiles to use

addpath(fullfile(rootdir,'MFILES'));

%%%%

for ccc=1:length(corefiles)

    corefile=load(fullfile(rootdir,['IFILES/' corefiles{ccc}]));
    if length(subcore{ccc})>0
        corefile=corefile.(subcore{ccc});
    end
    ccclab=corefilelabs{ccc};

    % pull GSL samples
    [sampsGSLrise,sampsGSLcomponents,GSLsiteids,GSLsitenames,targyears,GSLscens,GSLcols] = LocalizeStoredProjections(0,corefile);
    WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,targyears,GSLcols,GSLscens,['GSLproj_decomp_' ccclab '_']);

    for selectedSite=selectedSites(:)'

        % generate local samples

        [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);
        nameshort=sitenames{1}(1:3);

        % separate out thermal expansion and DSL
        for nnn=1:length(sampsGSLcomponents)
            sampsloccomponents{nnn}(:,end+1,:)=sampsloccomponents{nnn}(:,cols.colOD,:)-sampsGSLcomponents{nnn}(:,GSLcols.colTE,:);
            sampsloccomponents{nnn}(:,cols.colTE,:)=sampsGSLcomponents{nnn}(:,GSLcols.colTE,:);
        end
        cols.colOD=size(sampsloccomponents{1},2);

        % plot curves

        clf;
        if length(scens)<4
            subscens=1:length(scens);
        else
            subscens=[1 3 4];
        end
        [hp,hax,hax2]=PlotSLRProjection(sampslocrise,targyears,1,scens,subscens);
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

        if ismember(2010,targyears)
            clf;
            [hp,vars,fvars,hlg]=PlotSLRProjectionVariance(sampsloccomponents(:,subscens),targyears,cols,[2010 2100],1,min(2,length(subscens)),1,subcomp,complabls,'rcbgmykrcbgm');
            title(sitenames{1});
            pdfwrite(['LSLprojvar_' ccclab '_' nameshort '_' num2str(siteids(1))]);
        end

        % output quantiles of projections

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
        scensel=1; timesel=find(targyears==2050);
        PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);
        subplot(3,1,2);
        scensel=1; timesel=find(targyears==2100);
        PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
        xl=get(gca,'xlim');
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);

        if length(subscens)>1
            subplot(3,1,3);
            scensel=subscens(end); timesel=find(targyears==2100);
            PlotSLRDecomposition(sampsloccomponents,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
            set(gca,'xlim',xl);
            title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab]);
        end
        pdfwrite(['LSLproj_decomp_' ccclab '_' nameshort]);

    end

end