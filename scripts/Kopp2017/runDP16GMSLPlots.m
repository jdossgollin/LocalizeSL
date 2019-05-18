% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-10-25 18:36:45 -0400

% first run runDP16subsets; we will use the 'full' subset

colrs=[ 27 158 119 ; 217 95 2 ; 117 112 179 ]/255;
doscens=[4 3 1];
endyears=[2100 2300];
ylims=[0 2.5 ; 0 16];
titls={'K14','','DP16',''};
scenlabs={'RCP 2.6','RCP 4.5','RCP 8.5'};
letrs='abcd';
clf;
for yyy=1:2
    if yyy==1
        doGMSL=sampsGSLrise;
    else
        doGMSL=sampsGSLrise2sub{1};
    end
for xxx=1:length(endyears)
    hax(xxx+2*(yyy-1))=subplot(4,1,xxx+2*(yyy-1));
    if xxx==2
        plot([2000 2100],[1 1]*ylims(1,2),'color',[.5 .5 .5]); hold on;
        plot([2100 2100],[0 1]*ylims(1,2),'color',[.5 .5 .5]);        
    end
    clear hp;
    for aaa=1:length(doscens)
        hp(aaa)=plot([2000 targyearsGSL],[0 median(doGMSL{doscens(aaa)})]/1000,'r','linew',2,'color',colrs(aaa,:)); hold on;
    end

    xlim([2000 endyears(xxx)]);
    if (xxx+2*(yyy-1))==1
        legend(hp(end:-1:1),scenlabs(end:-1:1),'Location','Northeast');
    end
    ha=hax(xxx+2*(yyy-1));
    longticks(gca,2); ylabel('GMSL (m)');
    if xxx==2
        xlabel('Year');
    end
    title(titls{xxx+2*(yyy-1)});
    pos=get(ha,'Position'); pos0=pos;
    pos(3)=pos(3)*.9;
    set(ha,'Position',pos)

    pos2=[pos(1)+pos(3) pos(2) (pos0(3)-pos(3)) pos(4)];
    ha2=axes('position',pos2);
    set(ha2,'xlim',[.5 3.5],'box','off','color','none','ycolor','none','xcolor','none','xlim',[0 4]);

    axes(ha2);
    minq=0; maxq=0;
    for aaa=1:length(doscens)
        q=quantile(doGMSL{doscens(aaa)}(:,find(targyearsGSL==endyears(xxx))),[.05 .95])/1000;
        hr=rectangle('position',[.75+(aaa-1) q(1) .5 q(2)-q(1)],'edgecolor',colrs(aaa,:),'facecolor',colrs(aaa,:));
        minq=min(min(q),minq); maxq=max(max(q),max(q));
    end
    set(ha,'ylim',ylims(xxx,:)); set(ha2,'ylim',ylims(xxx,:));
    hax2(xxx+2*(yyy-1))=ha2;
    u=get(ha,'ytick');

    if sum(abs(u-round(u)))>0
        ytl=get(ha,'yticklabel');
        for ppp=1:length(ytl)
            ytl{ppp}=sprintf('%0.1f',u(ppp));
        end
        set(ha,'yticklabel',ytl);
    end


    axes(ha);
    ht=text(2000+[endyears(xxx)-2000]*.025,ylims(xxx,2)*.85,letrs(xxx+2*(yyy-1)));
    set(ht,'fontweight','bold','fontsize',12,'verticalalignment','bottom');
end
end
movev([hax(3:4) hax2(3:4)],-.03);
set(gcf,'PaperPosition',[.25 .5 8 8]);
pdfwrite('fig3-GMSLoldandnew');

%%% filtered projections

titls={'K14','DP16'};

clear q;
targlevels=[500 2000]; targwins=[100 100];
clear subA;
for yyy=1:2
    if yyy==1
        doGMSL=sampsGSLrise;
    else
        doGMSL=sampsGSLrise2sub{1};
    end
    doGMSL2=[doGMSL{1} ; doGMSL{3} ; doGMSL{4}];
    for sss=1:length(targlevels)
      sub=find(abs(doGMSL2(:,find(targyearsGSL==2100))-targlevels(sss))<targwins(sss));
        q(:,:,yyy,sss)=quantile(doGMSL2(sub,:),[.5 .05 .95])/10;
    end
end

clf;
plotdat.x=[2000 targyearsGSL]';
for yyy=1:2
    subplot(3,1,yyy);

    for sss=1:2
        plotdat.y=[0 q(1,:,yyy,sss)]';
        plotdat.dy=[0 q(1,:,yyy,sss)-q(2,:,yyy,sss) ; 0 q(3,:,yyy,sss)-q(1,:,yyy,sss)]';
        [hl(sss),hk]=PlotWithShadedErrors(plotdat,colrs(sss,:),.8,'-','--',[2000 2100],[0 250]);
    end
    if yyy==1
        legend(hl([2 1]),{'200 cm','50 cm'},'Location','Northwest');
    end
    for sss=1:2
        plotdat.y=[0 q(1,:,yyy,sss)]';
        plotdat.dy=[0 q(1,:,yyy,sss)-q(2,:,yyy,sss) ; 0 q(3,:,yyy,sss)-q(1,:,yyy,sss)]';
        PlotWithShadedErrors(plotdat,colrs(sss,:),'none','none','--',[2000 2100],[0 250]);
    end
    box on;
    longticks(gca,2);
    ylabel('GMSL (cm)');
    if yyy==2
        xlabel('Year');
    end
    title(titls{yyy});
end
pdfwrite('condscens');

% now invert to get learning under different pathways
for yyy=1:2
    if yyy==1
        doGMSL=sampsGSLrise;
    else
        doGMSL=sampsGSLrise2sub{1};
    end
    doGMSL2=[doGMSL{1} ; doGMSL{3} ; doGMSL{4}];
    targt=find(targyearsGSL==2100);
    for sss=1:2;
        for ttt=1:targt
            bnds=10*q(2:3,ttt,yyy,sss);
            sub=find((doGMSL2(:,ttt)>=bnds(1)).*(doGMSL2(:,ttt)<=bnds(2)));
            q2(:,ttt,yyy,sss)=quantile(doGMSL2(sub,targt),[.5 .05 .95])/10;
        end
    q20(:,yyy,sss)=quantile(doGMSL2(:,targt),[.5 .05 .95])/10;
    end
end

clear plotdat;
clf;
plotdat.x=[2000 targyearsGSL(1:targt)]';
for yyy=1:2
    subplot(3,1,yyy);
    for sss=1:2
        plotdat.y=[q20(1,yyy,sss) q2(1,:,yyy,sss)]';
        plotdat.dy=[[q20(1,yyy,sss) q2(1,:,yyy,sss)]-[q20(2,yyy,sss) q2(2,:,yyy,sss)] ; [q20(3,yyy,sss) q2(3,:,yyy,sss)]-[q20(1,yyy,sss) q2(1,:,yyy,sss)]]';
        [hl(sss),hk]=PlotWithShadedErrors(plotdat,colrs(sss,:),.8,'-','--',[2000 2100],[0 250]);
    end
    if yyy==1
        legend(hl([2 1]),{'200 cm','50 cm'},'Location','Northwest');
    end

    for sss=1:2
        plotdat.y=[q20(1,yyy,sss) q2(1,:,yyy,sss)]';
        plotdat.dy=[[q20(1,yyy,sss) q2(1,:,yyy,sss)]-[q20(2,yyy,sss) q2(2,:,yyy,sss)] ; [q20(3,yyy,sss) q2(3,:,yyy,sss)]-[q20(1,yyy,sss) q2(1,:,yyy,sss)]]';
        PlotWithShadedErrors(plotdat,colrs(sss,:),'none','-','--',[2000 2100],[0 250]);
    end
    box on;
    longticks(gca,2);
    ylabel('GMSL (cm)');
    if yyy==2
        xlabel('Year');
    end
    title(titls{yyy});
end

% plot to show learning


titls={'K14','DP16'};


clear plotdat;
for yyy=1:2
    clf;
    subplot(3,1,1);
    plotdat.x=[2000 targyearsGSL]';
    for sss=1:2
        plotdat.y=[0 q(1,:,yyy,sss)]';
        plotdat.dy=[0 q(1,:,yyy,sss)-q(2,:,yyy,sss) ; 0 q(3,:,yyy,sss)-q(1,:,yyy,sss)]';
        [hl(sss),hk]=PlotWithShadedErrors(plotdat,colrs(sss,:),.8,'-','--',[2000 2100],[0 250]);
    end
    if yyy==1
        legend(hl([2 1]),{'200 cm','50 cm'},'Location','Northwest');
    end
    for sss=1:2
        plotdat.y=[0 q(1,:,yyy,sss)]';
        plotdat.dy=[0 q(1,:,yyy,sss)-q(2,:,yyy,sss) ; 0 q(3,:,yyy,sss)-q(1,:,yyy,sss)]';
        PlotWithShadedErrors(plotdat,colrs(sss,:),'none','-','--',[2000 2100],[0 250]);
    end
    box on;
    longticks(gca,2);
    ylabel('GMSL (cm)');
    title(titls{yyy});

    subplot(3,1,2)
    plotdat.x=[2000 targyearsGSL(1:targt)]';
    for sss=1:2
        plotdat.y=[q20(1,yyy,sss) q2(1,:,yyy,sss)]';
        plotdat.dy=[[q20(1,yyy,sss) q2(1,:,yyy,sss)]-[q20(2,yyy,sss) q2(2,:,yyy,sss)] ; [q20(3,yyy,sss) q2(3,:,yyy,sss)]-[q20(1,yyy,sss) q2(1,:,yyy,sss)]]';
        [hl(sss),hk]=PlotWithShadedErrors(plotdat,colrs(sss,:),.8,'-','--',[2000 2100],[0 250]);
    end
 
    for sss=1:2
        plotdat.y=[q20(1,yyy,sss) q2(1,:,yyy,sss)]';
        plotdat.dy=[[q20(1,yyy,sss) q2(1,:,yyy,sss)]-[q20(2,yyy,sss) q2(2,:,yyy,sss)] ; [q20(3,yyy,sss) q2(3,:,yyy,sss)]-[q20(1,yyy,sss) q2(1,:,yyy,sss)]]';
        PlotWithShadedErrors(plotdat,colrs(sss,:),'none','-','--',[2000 2100],[0 250]);
    end
    box on;
    longticks(gca,2);
    ylabel('Proj. 2100 GMSL (cm)');
    if yyy==2
        xlabel('Year');
    end
    pdfwrite(['learning-' titls{yyy}])
end