function [hp,hax,hax2]=PlotSLRProjection(sampslocrise,targyears,sitesel,scens,subscens,endyears)

% PlotSLRProjection(sampslocrise,targyears,[sitesel],[scens],[subscens],[endyears])
%
% Plot sea-level rise projections.
%
% sampslocrise, targyears, and scens are outputted by LocalizeStoredProjections.
% sitesel is the row id of the site of interest in sampslocrise (default: 1)
% subscens can be used to select scenarios by sequential id (default: [1 3 4])
% endyears is used to select the end years of the plots (default: [2100 2200])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-01-16 16:59:49 -0500


defval('sitesel',1);
defval('scens',{'RCP85','RCP60','RCP45','RCP26'});
defval('subscens',[4 3 1]);
colrs=[ 27 158 119 ; 217 95 2 ; 117 112 179 ]/255;
defval('endyears',[2100 2200]);
defval('startyear',[2000]);

scenlab={};
for kk=subscens
	scenlab={scenlab{:},[upper(scens{kk}(1:3)) ' ' scens{kk}(4) '.' scens{kk}(5)]};
end

for nn=1:length(endyears)
    for vvvv=1:length(subscens)
        qvals(vvvv,:,nn)=quantile(sampslocrise{sitesel,subscens(vvvv)}(:,find(targyears==endyears(nn))),[.05 .95]);
        qvals(vvvv,:,nn)=quantile(sampslocrise{sitesel,subscens(vvvv)}(:,find(targyears==endyears(nn))),[.05 .95]);
    end
    ylims(nn,:)=[floor(min(min(squeeze(qvals(:,1,nn))*2/1000)))/2 ceil(max(max(squeeze(qvals(:,2,nn))*2/1000)))/2];
end
ylims(:,1)=min(0,ylims(:,1));

clf;

for xxx=1:length(endyears)
    hax(xxx)=subplot(2,1,xxx);
    plot([startyear endyears(xxx)],[0 0],'color',[.5 .5 .5]); hold on
    if xxx>1
        plot([startyear endyears(xxx-1)],[1 1]*ylims(xxx-1,2),'color',[.5 .5 .5]); hold on;
        plot([endyears(xxx-1) endyears(xxx-1)],[0 1]*ylims(xxx-1,2),'color',[.5 .5 .5]);        
    end
    clear hp;
    for aaa=1:length(subscens)
        hp(aaa)=plot([startyear targyears],[0 median(sampslocrise{sitesel,subscens(aaa)})]/1000,'linew',2,'color',colrs(aaa,:)); hold on;
    end

    xlim([startyear endyears(xxx)]);
    if (xxx)==1
        legend(hp(end:-1:1),scenlab(end:-1:1),'Location','Northeast');
    end
    ha=hax(xxx);
    longticks(gca,2); ylabel('RSL (m)');
    if xxx==length(endyears)
        xlabel('Year');
    end
    pos=get(ha,'Position'); pos0=pos;
    pos(3)=pos(3)*.9;
    set(ha,'Position',pos)

    pos2=[pos(1)+pos(3) pos(2) (pos0(3)-pos(3)) pos(4)];
    ha2=axes('position',pos2);
    set(ha2,'xlim',[.5 3.5],'box','off','color','none','ycolor','none','xcolor','none','xlim',[0 4]);

    axes(ha2);
    minq=0; maxq=0;
    for aaa=1:length(subscens)
        q=quantile(sampslocrise{sitesel,subscens(aaa)}(:,find(targyears==endyears(xxx))),[.05 .95])/1000;
        hr=rectangle('position',[.75+(aaa-1) q(1) .5 q(2)-q(1)],'edgecolor',colrs(aaa,:),'facecolor',colrs(aaa,:));
        minq=min(min(q),minq); maxq=max(max(q),max(q));
    end
    set(ha,'ylim',ylims(xxx,:)); set(ha2,'ylim',ylims(xxx,:));
    hax2(xxx)=ha2;
    u=get(ha,'ytick');

    if sum(abs(u-round(u)))>0
        ytl=get(ha,'yticklabel');
        for ppp=1:length(ytl)
            ytl{ppp}=sprintf('%0.1f',u(ppp));
        end
        set(ha,'yticklabel',ytl);
    end

    axes(ha);
end

