    % 2019-11-26 15:24:18 -0500

    extremesIFILES = fullfile(rootdir,'IFILES/extremes');
    datDir = fullfile(extremesIFILES, 'declustered');
    % import table of calibration parameters
    parmdat=importdata(fullfile(extremesIFILES,'GPDfits_withunc.tsv'),'\t',1);
    extracols=size(parmdat.textdata,2)-size(parmdat.data,2);
    lambdas=parmdat.data(:,find(strcmpi('lambda',parmdat.textdata(1,:)))-extracols);
    thresholds=parmdat.data(:,find(strcmpi('u',parmdat.textdata(1,:)))-extracols);
    scales=parmdat.data(:,find(strcmpi('scale',parmdat.textdata(1,:)))-extracols);;
    shapes=parmdat.data(:,find(strcmpi('shape',parmdat.textdata(1,:)))-extracols);
    AEP10pts=parmdat.data(:,find(strcmpi('AEP0.1',parmdat.textdata(1,:)))-extracols);
    Vscale=parmdat.data(:,find(strcmpi('Vscale',parmdat.textdata(1,:)))-extracols);
    Vshape=parmdat.data(:,find(strcmpi('Vshape',parmdat.textdata(1,:)))-extracols);
    Vscaleshape=parmdat.data(:,find(strcmpi('Vscaleshape',parmdat.textdata(1,:)))-extracols);
    psmslids=parmdat.data(:,1);
    tgids=parmdat.data(:,2);
    NOAAnames=parmdat.textdata(2:end,1);

    clear effcurve histcurve effcurve999 integratecurve;
    qqq=find(psmslids==selectedSite);
    wfile=fullfile(datDir,['maxtofit.dclist.' num2str(tgids(qqq)) '_xdat.dc.tsv']);
    historicaldata=importdata(wfile);

    longname=NOAAnames{qqq};
    shortname=longname(setdiff(1:length(longname),strfind(longname,' ')));
    shortname(strfind(shortname,','))='-';

    acov = [Vscale(qqq) Vscaleshape(qqq) ; Vscaleshape(qqq) Vshape(qqq)];
    parmsamps=lhsnorm([scales(qqq) shapes(qqq)],acov,1000);
    parmsamps(:,1)=max(eps,parmsamps(:,1));

    clear pm;
    pm.showuncertainty=1; pm.historicaldata=historicaldata;
    pm.showESLR=0; pm.showeffcurve=1; pm.show999=0; pm.showinteffcurve=0;
    pm.endyears=[2050 2100];

    clear effcurve;
    for iii=1:4

            samps=[zeros(size(storesampslocrise{iii},1),1) storesampslocrise{iii}]/1000; % add base year and convert to meters

            % truncate samples
            [s,si]=sort(storesampslocrise{iii}(:,find(targyears==2100)));
            [mi]=find(s>3000);
            if length(mi)>0
                subsi=(length(si)-mi(1)):mi(1);
                samps=samps(si(subsi),:);
            end
    %        samps=bsxfun(@min,samps,quantile(samps,.999)); % truncate samples viewed as physically implausible
            targyears2 = [2000 targyears];

        clf;

            [effcurve{iii},testz,histcurve,histcurvesamps]= ...
            SLRFloodNexpVsLevelCurves(samps,targyears2,thresholds(qqq), ...
            parmsamps(:,1),parmsamps(:,2),lambdas(qqq),longname,pm);     
    end

    historicalcolor='y';

    clear hl;
    clf;
    subplot(2,1,1);
    ct=sum(bsxfun(@gt,historicaldata,testz))/(length(historicaldata)/365.25);
    subct=find(diff(ct)<0);
    subct=intersect(subct,find(testz>thresholds(qqq)));
    plot(testz(subct),ct(subct),'s','Color',historicalcolor); hold on;
    hl(1)=plot(testz,histcurve,'k');
    plot(testz,quantile(histcurvesamps,[.17 .83],1),'color',[.6 .6 .6]);
    hold on; 
    set(gca,'yscale','log'); ylim([.002 10]);

    legstr={'RCP8.5-K14','2C-R18','RCP8.5-B19','2C-B19'};

    linespecs={'-','-','--','--'};
    colorspecs=[217 95 2 ; 27 158 119 ; 217 95 2 ; 27 158 119]/255;


    for iii=1:length(legstr)
        subyr=find(targyears==2050)+1;
        hl(end+1)=plot(testz,effcurve{iii}(subyr,:),'linestyle',linespecs{iii},'color',colorspecs(iii,:)); hold on;
    end
    ylabel({'Expected number of','annual exceedances'});
    xlabel('Extreme sea level (m)');
    title([longname ' - 2050']);
    u=legstr';
    u={'Historical',u{1:end}};
    hld=legend(hl,u,'location','northeast')
    set(hld,'fontsize',7)
    xlim([0 4.75]);

    subplot(2,1,2);
    ct=sum(bsxfun(@gt,historicaldata,testz))/(length(historicaldata)/365.25);
    subct=find(diff(ct)<0);
    subct=intersect(subct,find(testz>thresholds(qqq)));
    plot(testz(subct),ct(subct),'s','Color',historicalcolor); hold on;
    plot(testz,histcurve,'k');
    plot(testz,quantile(histcurvesamps,[.17 .83],1),'color',[.6 .6 .6]);
    hold on; 
    set(gca,'yscale','log'); ylim([.002 10]);


    for iii=1:length(legstr)
        subyr=find(targyears==2100)+1;
        hl(end+1)=plot(testz,effcurve{iii}(subyr,:),'linestyle',linespecs{iii},'color',colorspecs(iii,:)); hold on;
    end

    ylabel({'Expected number of','annual exceedances'});
    xlabel('Extreme sea level (m)');
    title([longname ' - 2100']);
    xlim([0 4.75]);
    pdfwrite(['ESL_' shortname]);

    fid=fopen(['ESL_' shortname '.tsv'],'w');
    fprintf(fid,'\t%0.3f',testz);
    fprintf(fid,'\n');
    fprintf(fid,'historical');
    fprintf(fid,'\t%0.3g',histcurve);
    fprintf(fid,'\n');
    fprintf(fid,'historical - 17th');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps,.17));
    fprintf(fid,'\n');
    fprintf(fid,'historical - 83rd');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps,.83));
    fprintf(fid,'\n');
    for doyr=2010:10:2150;
        for iii=1:length(effcurve)
            subyr=find(targyears==doyr)+1;
            fprintf(fid,[legstr{iii} ' - %0.0f'],doyr);
            fprintf(fid,'\t%0.3g',effcurve{iii}(subyr,:));
            fprintf(fid,'\n');
        end
    end

    fclose(fid);

    doAAA=[5.7 1 .1 .01];
    subdo=[];
    for ddd=1:length(doAAA)
        [m,mi]=min(abs(histcurve-doAAA(ddd)));
        subdo=[subdo mi];
    end

    fid=fopen(['ESL2_' shortname '.tsv'],'w');
    fprintf(fid,'\t%0.3f',testz(subdo));
    fprintf(fid,'\n');
    fprintf(fid,'historical');
    fprintf(fid,'\t%0.3g',histcurve(subdo));
    fprintf(fid,'\n');
    fprintf(fid,'historical - 17th');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps(:,subdo),.17));
    fprintf(fid,'\n');
    fprintf(fid,'historical - 83rd');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps(:,subdo),.83));
    fprintf(fid,'\n');
    for doyr=2010:10:2150;
        for iii=1:length(effcurve)
            subyr=find(targyears==doyr)+1;
            fprintf(fid,[legstr{iii} ' - %0.0f'],doyr);
            fprintf(fid,'\t%0.3g',effcurve{iii}(subyr,subdo));
            fprintf(fid,'\n');
        end
    end

    fclose(fid);

    %

    quantlevs=[.01 .05 .17 .5 .83 .95 .99 .995 .999]; % quantiles to use

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

    fid=fopen(['extremesealevel_quantiles_' shortname '.tsv'],'w');

    for ddd=1:length(doAAA)
        fprintf(fid,'Historical Nexp = %0.3f',doAAA(ddd));
        fprintf(fid,'\t%0.1f',quantlevs*100);
        fprintf(fid,'\n');
        fprintf(fid,['Early years (composite) \n']);
        subt=find(targyears<=2050);
        for ttt=subt
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledall(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        fprintf(fid,['2C (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals2C(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        fprintf(fid,['RCP 8.5 (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals85(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        for iii=subscen2C;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap'));
                fprintf(fid,'\n');
            end
        end
        for iii=subscen85;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap'));
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);

    %%%%%%

    doAAA=[.01];
    subdo=[];
    for ddd=1:length(doAAA)
        [m,mi]=min(abs(histcurve-doAAA(ddd)));
        subdo=[subdo mi];
    end

    fid=fopen(['extremesealevelAFs_quantiles_' shortname '.tsv'],'w');

    for ddd=1:length(doAAA)
        fprintf(fid,'Historical Nexp = %0.3f',doAAA(ddd));
        fprintf(fid,'\t%0.1f',quantlevs*100);
        fprintf(fid,'\n');
        fprintf(fid,['Early years (composite) \n']);
        subt=find(targyears<=2050);
        for ttt=subt
            fprintf(fid,'%0.0f',targyears(ttt));
            newlev=interp1(testz,histcurve,testz(subdo(ddd))-pooledall(:,ttt)/1000,'nearest','extrap');
            fprintf(fid,'\t%0.3f',newlev/doAAA(ddd));
            fprintf(fid,'\n');
        end
        fprintf(fid,['2C (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            newlev=interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals2C(:,ttt)/1000,'nearest','extrap');
            fprintf(fid,'\t%0.3f',newlev/doAAA(ddd));
            fprintf(fid,'\n');
        end
        fprintf(fid,['RCP 8.5 (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            newlev=interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals85(:,ttt)/1000,'nearest','extrap');
            fprintf(fid,'\t%0.3f',newlev/doAAA(ddd));
            fprintf(fid,'\n');
        end
        fprintf(fid,['Moderate (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            newlev=interp1(testz,histcurve,testz(subdo(ddd))-.5*(pooledqvals85(:,ttt)+pooledqvals2C(:,ttt))/1000,'nearest','extrap');
            fprintf(fid,'\t%0.3f',newlev/doAAA(ddd));
            fprintf(fid,'\n');
        end
        for iii=subscen2C;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap');
                fprintf(fid,'\t%0.3f',newlev/doAAA(ddd));
                fprintf(fid,'\n');
            end
        end
        for iii=subscen85;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap')/doAAA(ddd));
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);


    %%%%%%%

    dohts=[0.56 0.84 1.23]
    subdo=[];
    for ddd=1:length(dohts)
        [m,mi]=min(abs(testz-dohts(ddd)));
        subdo=[subdo mi];
    end

    fid=fopen(['ESL3_' shortname '.tsv'],'w');
    fprintf(fid,'\t%0.3f',testz(subdo));
    fprintf(fid,'\n');
    fprintf(fid,'historical');
    fprintf(fid,'\t%0.3g',histcurve(subdo));
    fprintf(fid,'\n');
    fprintf(fid,'historical - 17th');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps(:,subdo),.17));
    fprintf(fid,'\n');
    fprintf(fid,'historical - 83rd');
    fprintf(fid,'\t%0.3g',quantile(histcurvesamps(:,subdo),.83));
    fprintf(fid,'\n');
    for doyr=2010:10:2150;
        for iii=1:length(effcurve)
            subyr=find(targyears==doyr)+1;
            fprintf(fid,[legstr{iii} ' - %0.0f'],doyr);
            fprintf(fid,'\t%0.3g',effcurve{iii}(subyr,subdo));
            fprintf(fid,'\n');
        end
    end

    fclose(fid);

    %

    quantlevs=[.01 .05 .17 .5 .83 .95 .99 .995 .999]; % quantiles to use

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

    fid=fopen(['extremesealevel_quantiles2_' shortname '.tsv'],'w');

    for ddd=1:length(dohts)
        fprintf(fid,'Height = %0.2f m',dohts(ddd));
        fprintf(fid,'\t%0.1f',quantlevs*100);
        fprintf(fid,'\n');
        fprintf(fid,['Early years (composite) \n']);
        subt=find(targyears<=2050);
        for ttt=subt
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledall(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        fprintf(fid,['2C (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals2C(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        fprintf(fid,['RCP 8.5 (composite) \n']);
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-pooledqvals85(:,ttt)/1000,'nearest','extrap'));
            fprintf(fid,'\n');
        end
        for iii=subscen2C;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap'));
                fprintf(fid,'\n');
            end
        end
        for iii=subscen85;
            fprintf(fid,[newscenlabs{iii} '\n']);
            for ttt=1:length(targyears)
                fprintf(fid,'%0.0f',targyears(ttt));
                fprintf(fid,'\t%0.3f',interp1(testz,histcurve,testz(subdo(ddd))-qvals(:,ttt,iii)/1000,'nearest','extrap'));
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);