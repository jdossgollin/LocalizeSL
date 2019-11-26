% 2019-11-26 15:23:00 -0500

rootdir='~/Dropbox/Code/LocalizeSL'; % change to directory containing LocalizeSL
corefiles={'SLRProjections170113GRIDDEDcore','SLRProjections180124GRIDDEDcore_Tscens.mat','SLRProjections190726core_SEJ_full.mat','SLRProjections190726core_SEJ_full.mat'}; % specify corefiles to use
corefilelabs={'K14','R18','B19H','B19L'}; % specify corefile labels
subcore={'','','corefileH','corefileL'}; % specify if corefile file contains multiple cores
scenselect={[1],[2],[1],[1]};
quantlevs=[.01 .05 .17 .5 .83 .95 .99 .995 .999]; % quantiles to use

addpath(fullfile(rootdir,'MFILES'));

workdir='workdir-191126-proj';
if ~exist(workdir,'dir')
    mkdir(workdir);
end
cd(workdir);

selectedSites=[396 1636];
for sssss=1:length(selectedSites)
    selectedSite=selectedSites(sssss);

    clear storesampslocrise;

    for ccc=1:length(corefiles)

        corefile=load(fullfile(rootdir,['IFILES/' corefiles{ccc}]));
        if length(subcore{ccc})>0
            corefile=corefile.(subcore{ccc});
        end
        ccclab=corefilelabs{ccc};

        [sampsGSLrise,sampsGSLcomponents,GSLsiteids,GSLsitenames,targyears,GSLscens,GSLcols] = LocalizeStoredProjections(0,corefile,scenselect{ccc});

        [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,scenselect{ccc});
        nameshort=sitenames{1}(1:3);

        % separate out thermal expansion and DSL
        for nnn=1:length(sampsGSLcomponents)
            sampsloccomponents{nnn}(:,end+1,:)=sampsloccomponents{nnn}(:,cols.colOD,:)-sampsGSLcomponents{nnn}(:,GSLcols.colTE,:);
            sampsloccomponents{nnn}(:,cols.colTE,:)=sampsGSLcomponents{nnn}(:,GSLcols.colTE,:);
        end
        cols.colOD=size(sampsloccomponents{1},2);

        % output quantiles of projections

        WriteTableSLRProjection(sampsGSLrise,quantlevs,GSLsiteids,GSLsitenames,targyears,scens,['GSLproj_' ccclab '_']);
    

        WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,['LSLproj_' ccclab '_' nameshort '_']);
        subcomp={cols.colAIS,cols.colGIS,cols.colGIC, cols.colTE, ...
        cols.colLS, cols.colOD, ...
        cols.colGIA,1:size(sampsloccomponents{1},2) };

        complabls={'AIS','GIS','GIC','TE','LWS','DSL','Geo','Total'};

        WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,['LSLproj_decomp_' ccclab '_' nameshort '_'],subcomp,complabls);
        WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,targyears,GSLcols,GSLscens,['GSLproj_decomp_' ccclab '_']);

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
        WriteTableSLRHeightExceedanceTiming(sampslocrise,[1:10]*12*2.54,siteids,sitenames,targyears,scens,1,['LSLheights_ft_' ccclab '_' nameshort '_']);

        % output Monte Carlo samples

        WriteTableMC(sampsloccomponents,[],siteids,sitenames,targyears,scens,['LSLproj_MC_' ccclab '_' nameshort '_']);

        % output Monte Carlo samples without background trend,
        % to allow incorporation of alternative estimates of background trend

        %WriteTableMC(sampsloccomponentsRE,setdiff(1:size(sampsloccomponents{1},2),cols.colGIA),siteids,sitenames,targyears,scens,['LSLproj_MC_nobkgd_' ccclab '_reweight2_' nameshort '_']);

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
        title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab ]);

        if length(subscens)>1
            subplot(3,1,3);
            scensel=subscens(end); timesel=find(targyears==2100);
            PlotSLRDecomposition(sampsloccomponentsRE,[],[],scensel,timesel,targyears,cols,subcomp,complabls);
            set(gca,'xlim',xl);
            title([sitenames{1} ' - ' scens{scensel} ' - ' num2str(targyears(timesel)) ' - ' ccclab ]);
        end
        pdfwrite(['LSLproj_decomp_' ccclab '_' nameshort]);

        WriteTableDecomposition(sampsGSLcomponents,quantlevs,GSLsiteids,GSLsitenames,targyears,GSLcols,GSLscens,['GSLproj_decomp_' ccclab '_']);
        storesampslocrise{ccc}=sampslocrise{1};

    end

    %%%%


    %%%%

    newscenlabs={'RCP8.5-K14','2C-R18','RCP8.5-B19','2C-B19'};
    subscen2C=[2 4];
    subscen85=[1 3];

    clf;
    [hp,hax,hax2]=PlotSLRProjection(storesampslocrise,targyears,1,newscenlabs,subscen2C);
    axes(hax(1));
    title([sitenames{1} ' - 2C']);
    pdfwrite(['LSLproj_2C_' nameshort '_' num2str(siteids(1)) ]);

    clf;
    [hp,hax,hax2]=PlotSLRProjection(storesampslocrise,targyears,1,newscenlabs,subscen85);
    axes(hax(1));
    title([sitenames{1} ' - RCP 8.5']);
    pdfwrite(['LSLproj_RCP85_' nameshort '_' num2str(siteids(1)) ]);

    clear qvals;
    clear pooledqvals2C pooledqvals85 pooledall;
    for ccc=1:length(storesampslocrise)
        qvals(:,:,ccc)=quantile(storesampslocrise{ccc},quantlevs);
    end
    subt=find(targyears<=2050);
    sub=find(quantlevs<.5);
    pooledqvals2C(sub,:)=min(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=min(qvals(sub,:,subscen85),[],3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=min(comp,[],3);
    sub=find(quantlevs==.5);
    pooledqvals2C(sub,:)=mean(qvals(sub,:,subscen2C),3);
    pooledqvals85(sub,:)=mean(qvals(sub,:,subscen85),3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=mean(comp,3);
    sub=find(quantlevs>.5);
    pooledqvals2C(sub,:)=max(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=max(qvals(sub,:,subscen85),[],3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=max(comp,[],3);

    scalefactor=0.1; units='cm';
    fid=fopen(['sealevel_quantiles_' nameshort '.tsv'],'w');
    fprintf(fid,units);
    fprintf(fid,'\t%0.1f',quantlevs*100);
    fprintf(fid,'\n');
    fprintf(fid,['Early years (composite) \n']);
    subt=find(targyears<=2050);
    for ttt=subt
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledall(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['2C (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledqvals2C(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['RCP 8.5 (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledqvals85(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    for iii=subscen2C;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    for iii=subscen85;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);

    scalefactor=0.1/2.54/12; units='ft';
    fid=fopen(['sealevel_quantiles_ft_' nameshort '.tsv'],'w');
    fprintf(fid,units);
    fprintf(fid,'\t%0.1f',quantlevs*100);
    fprintf(fid,'\n');
    fprintf(fid,['Early years (composite) \n']);
    subt=find(targyears<=2050);
    for ttt=subt
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.1f',pooledall(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['2C (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals2C(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['RCP 8.5 (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals85(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    for iii=subscen2C;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    for iii=subscen85;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);


    %%% 
    % derivatives

    dtstep=4;

    difftargyears=(targyears(1:end-dtstep)+targyears((1+dtstep):end))/2;

    clear qvals;
    clear pooledqvals2C pooledqvals85;
    clear diffstoresampslocrise;
    for ccc=1:length(storesampslocrise)
        diffstoresampslocrise{ccc}=storesampslocrise{ccc}(:,(1+dtstep):end)-storesampslocrise{ccc}(:,1:end-dtstep);
        diffstoresampslocrise{ccc}=diffstoresampslocrise{ccc}/(targyears(1+dtstep)-targyears(1));
    %    subt=find(difftargyears==2100);
    %    diffstoresampslocrise{ccc}(:,subt+[0 1])=repmat(mean(diffstoresampslocrise{ccc}(:,subt+[0 1]),2),1,2);
        qvals(:,:,ccc)=quantile(diffstoresampslocrise{ccc},quantlevs);
    end
    %subt=find(difftargyears<=2040);
    sub=find(quantlevs<.5);
    pooledqvals2C(sub,:)=min(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=min(qvals(sub,:,subscen85),[],3);
    %clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    %pooledall(sub,:)=min(comp,[],3);

    sub=find(quantlevs==.5);
    pooledqvals2C(sub,:)=mean(qvals(sub,:,subscen2C),3);
    pooledqvals85(sub,:)=mean(qvals(sub,:,subscen85),3);
    %clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    %pooledall(sub,:)=mean(comp,3);

    sub=find(quantlevs>.5);
    pooledqvals2C(sub,:)=max(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=max(qvals(sub,:,subscen85),[],3);
    %clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    %pooledall(sub,:)=max(comp,[],3);


    scalefactor=1; units='mm/yr';
    fid=fopen(['sealevelrate_quantiles_' nameshort '.tsv'],'w');
    fprintf(fid,units);
    fprintf(fid,'\t%0.1f',quantlevs*100);
    fprintf(fid,'\n');
    fprintf(fid,['2C (composite) \n']);
    for ttt=1:length(difftargyears)
        fprintf(fid,'%0.0f',difftargyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals2C(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['RCP 8.5 (composite) \n']);
    for ttt=1:length(difftargyears)
        fprintf(fid,'%0.0f',difftargyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals85(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    for iii=subscen2C;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(difftargyears)
            fprintf(fid,'%0.0f',difftargyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    for iii=subscen85;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(difftargyears)
            fprintf(fid,'%0.0f',difftargyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);

    scalefactor=1/2.54; units='in/decade';
    fid=fopen(['sealevelrate_quantiles_in_' nameshort '.tsv'],'w');
    fprintf(fid,units);
    fprintf(fid,'\t%0.1f',quantlevs*100);
    fprintf(fid,'\n');
    fprintf(fid,['2C (composite) \n']);
    for ttt=1:length(difftargyears)
        fprintf(fid,'%0.0f',difftargyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals2C(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['RCP 8.5 (composite) \n']);
    for ttt=1:length(difftargyears)
        fprintf(fid,'%0.0f',difftargyears(ttt));
        fprintf(fid,'\t%0.1f',pooledqvals85(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    for iii=subscen2C;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(difftargyears)
            fprintf(fid,'%0.0f',difftargyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    for iii=subscen85;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(difftargyears)
            fprintf(fid,'%0.0f',difftargyears(ttt));
            fprintf(fid,'\t%0.1f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);

    % fuller set of quantiles

    quantlevs2=[.001 .005 .01:.01:.99 .995 .999];

    clear qvals;
    clear pooledqvals2C pooledqvals85 pooledall;
    for ccc=1:length(storesampslocrise)
        qvals(:,:,ccc)=quantile(storesampslocrise{ccc},quantlevs2);
    end
    subt=find(targyears<=2050);
    sub=find(quantlevs2<.5);
    pooledqvals2C(sub,:)=min(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=min(qvals(sub,:,subscen85),[],3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=min(comp,[],3);
    sub=find(quantlevs2==.5);
    pooledqvals2C(sub,:)=mean(qvals(sub,:,subscen2C),3);
    pooledqvals85(sub,:)=mean(qvals(sub,:,subscen85),3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=mean(comp,3);
    sub=find(quantlevs2>.5);
    pooledqvals2C(sub,:)=max(qvals(sub,:,subscen2C),[],3);
    pooledqvals85(sub,:)=max(qvals(sub,:,subscen85),[],3);
    clear comp; comp(:,:,1)=pooledqvals2C(sub,:); comp(:,:,2)=pooledqvals85(sub,:);
    pooledall(sub,:)=max(comp,[],3);

    scalefactor=0.1; units='cm';
    fid=fopen(['sealevel_quantiles_full_' nameshort '.tsv'],'w');
    fprintf(fid,units);
    fprintf(fid,'\t%0.1f',quantlevs2*100);
    fprintf(fid,'\n');
    fprintf(fid,['Early years (composite) \n']);
    subt=find(targyears<=2050);
    for ttt=subt
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledall(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['2C (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledqvals2C(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    fprintf(fid,['RCP 8.5 (composite) \n']);
    for ttt=1:length(targyears)
        fprintf(fid,'%0.0f',targyears(ttt));
        fprintf(fid,'\t%0.0f',pooledqvals85(:,ttt)*scalefactor);
        fprintf(fid,'\n');
    end
    for iii=subscen2C;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    for iii=subscen85;
        fprintf(fid,[newscenlabs{iii} '\n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',qvals(:,ttt,iii)*scalefactor);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);


    %%%%%%%%

    if sssss==1
        runLocalizeESL
    end
end