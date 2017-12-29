% 2017-11-02 01:33:13 -0400
% AIS rate plots

doscens=[1 3 4];
sslab={'RCP 8.5','RCP 4.5','RCP 2.6'};
limpairs=[2000 2110 ; 2090 2300];
displimpairs=[2000 2100 ; 2100 2300];
%colrs=[0 0 255 ; 0 255 0 ; 255 0 255 ; 255 0 0 ]/255;
colrs=[ 27 158 119 ;  217 95 2 ; 117 112 179 ;  231 41 138 ]/255;

ppplabs={'icesub1','icesub2'};
pppsubs={[2 3 4 5],[6 8 9 7]};

rrrlabs={'AIS','WAIS','EAIS'};
bbblabs={'','noBC-'};
preferredBCvals=[2 1];

figure;
for bbb=1:2
preferredBC=preferredBCvals(bbb);
for rrr=1:3
    if rrr==1
        doAIScol=colsGSL.colAIS;
    else
        doAIScol=colsGSL.colAIS(rrr-1);
    end
    for ppp=1:length(pppsubs)
    clf;
    for ss=1:3
    for xx=1:2
        limpair=limpairs(xx,:);
        subplot(3,2,(ss-1)*2+xx); box on;

        sss=doscens(ss);
        subt=find((targyears<=limpair(2)));
        subt2=find((targyearsGSL<=limpair(2)));;

        plotdat2.x=[targyears(subt(1:end-1))+targyears(subt(2:end))]'/2;

        wdat0=sum(sampsGSLcomponents{sss}(:,doAIScol,subt2),2);
        wdat=zeros(size(wdat0,1),1,1);
        wdat(:,:,[1:length(subt2)]+1)=wdat0;
        wdat=diff(wdat,[],3)/10;

        y=squeeze(quantile(wdat,[.01 .5 .99],1))';
        plotdat2.y=y(:,2); plotdat2.dy=[[y(:,2)-y(:,1)] [y(:,3)-y(:,2)]];
        [hl,hk]=PlotWithShadedErrors(plotdat2,[0 0 0],.95,'none','none',limpair);

        plotdat.x=[targyears(subt(1:end-1))+targyears(subt(2:end))]'/2;
        y=squeeze(quantile(wdat,[.05 .5 .95],1))';
        plotdat.y=y(:,2); plotdat.dy=[[y(:,2)-y(:,1)] [y(:,3)-y(:,2)]];

        [hl,hk]=PlotWithShadedErrors(plotdat,[0 0 0],.8,'-','none',limpair);
        hold on;

        plot(plotdat.x,squeeze(quantile(wdat,.999)),'k--');

        deltayears=.5*(targyears(subt(1:end-1))+targyears(subt(2:end)));

        % $$$ doset=find(ensembleset(:,3));
        % $$$ plot(deltayears,diff([zeros(1,length(doset)) ; doRD{ppp,sss}(subt2,doset)])/10,'b-');
        % $$$ doset=find(ensembleset(:,2));
        % $$$ plot(deltayears,diff([zeros(1,length(doset)) ; doRD{ppp,sss}(subt2,doset)])/10,'r-');

        doicesets=pppsubs{ppp};
        for uuu=1:length(doicesets)
            if rrr==1
                doRD = RDWAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:)+RDEAISsub{doicesets(uuu)}{sss}(subt2,:);
            elseif rrr==2
                doRD = RDWAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:);
            elseif rrr==3
                doRD = RDEAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:);
            end
            plot(deltayears,diff([zeros(1,size(doRD,2)) ; doRD/10]),'Color',colrs(uuu,:)); hold on;
        end

        xlim(displimpairs(xx,:));
        longticks(gca,2);
        ylabel('mm/yr');
        if ss==3
            xlabel('Year');
        end
        if xx==1
            title([sslab{ss}]);
        end

        y=squeeze(quantile(wdat,[.01 .5 .99],1))';
        plot(deltayears,y,'k:');
        if xx==1
            ylim([-10 45]);
        elseif xx==2
            ylim([-20 120]);
        end


    end
    end
    set(gcf,'PaperPosition',[.25 2.5 8 6]);
    pdfwrite(['fig-rate-' rrrlabs{rrr} '-' ppplabs{ppp} '-' bbblabs{bbb} 'oldandnew']);
end
end

% AIS level plots

for rrr=1:3
    if rrr==1
        doAIScol=colsGSL.colAIS;
    else
        doAIScol=colsGSL.colAIS(rrr-1);
    end
for ppp=1:length(pppsubs)
clf;
for ss=1:3
for xx=1:2
    limpair=limpairs(xx,:);
    subplot(3,2,(ss-1)*2+xx); box on;

    sss=doscens(ss);
    subt=find((targyears<=limpair(2)));
    subt2=find((targyearsGSL<=limpair(2)));;

    plotdat2.x=targyears(subt)';
    y=squeeze(quantile(sum(sampsGSLcomponents{sss}(:,doAIScol,subt2),2),[.01 .5 .99],1))';
    y=[[0 0 0];y]/10; plotdat2.y=y(:,2); plotdat2.dy=[[y(:,2)-y(:,1)] [y(:,3)-y(:,2)]];
    [hl,hk]=PlotWithShadedErrors(plotdat2,[0 0 0],.9,'none','none',limpair);

    plotdat.x=targyears(subt)';
    y=squeeze(quantile(sum(sampsGSLcomponents{sss}(:,doAIScol,subt2),2),[.05 .5 .95],1))';
    y=[[0 0 0];y]/10; plotdat.y=y(:,2); plotdat.dy=[[y(:,2)-y(:,1)] [y(:,3)-y(:,2)]];

    [hl,hk]=PlotWithShadedErrors(plotdat,[0 0 0],.8,'-','none',limpair);
    hold on;

    plot(targyears(subt),[0 ; squeeze(quantile(sum(sampsGSLcomponents{sss}(:,doAIScol,subt2),2),[.999],1))]/10,'k--');

    % $$$ doset=find(ensembleset(:,3));
    % $$$ plot(deltayears,diff([zeros(1,length(doset)) ; doRD{ppp,sss}(subt2,doset)])/10,'b-');
    % $$$ doset=find(ensembleset(:,2));
    % $$$ plot(deltayears,diff([zeros(1,length(doset)) ; doRD{ppp,sss}(subt2,doset)])/10,'r-');

    doicesets=pppsubs{ppp};
    for uuu=1:length(doicesets)
        if rrr==1
            doRD = RDWAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:)+RDEAISsub{doicesets(uuu)}{sss}(subt2,:);
        elseif rrr==2
            doRD = RDWAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:);
        elseif rrr==3
            doRD = RDEAISsub{doicesets(uuu)}{preferredBC,sss}(subt2,:);
        end
        plot([2000 targyearsGSL(subt2)],[zeros(1,size(doRD,2)) ; doRD/10],'Color',colrs(uuu,:)); hold on;
    end

    xlim(displimpairs(xx,:));
    longticks(gca,2);
    ylabel('cm');
    if ss==3
        xlabel('Year');
    end
    if xx==1
        title([sslab{ss}]);
    end

    y=squeeze(quantile(sum(sampsGSLcomponents{sss}(:,doAIScol,subt2),2),[.01 .5 .99],1))';
    y=[0 0 0 ; y]/10;
    plot(targyears(subt),y(subt,:),'k:');

    % if xx==1
    %     ylim([-10 45]);
    % elseif xx==2
    %     ylim([-20 120]);
    % end


    
end
end
set(gcf,'PaperPosition',[.25 2.5 8 6]);
pdfwrite(['fig-level-' rrrlabs{rrr} '-' ppplabs{ppp} '-' bbblabs{bbb} 'oldandnew']);
end
end
end