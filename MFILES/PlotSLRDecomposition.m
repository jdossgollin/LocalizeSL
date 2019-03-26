function PlotSLRDecomposition(sampsloccomponents,quantlevs,sitesel,scensel,timesel,targyears,cols,subcomp,labls)
 
% PlotSLRDecomposition(sampsloccomponents,quantlevs,sitesel,sitenames,targyears,cols,[subcomp],[labls])
%
% Plot with bars decompositing into major components.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-03-25 20:39:47 -0400

defval('quantlevs',[.16 .84 ; .05 .95]);

if cols.colOD==cols.colTE
    cols.colTE2=[];
else
    cols.colTE2=cols.colTE;
end

defval('sitesel',1);
defval('scensel',1);
defval('timesel',10);
defval('subcomp',{cols.colAIS,cols.colGIS,[cols.colTE2 cols.colOD],[cols.colGIC],[cols.colLS],[cols.colGIA]});
defval('labls',{'AIS','GIS','Ocean','GIC','LWS','Bkgd'}); 

jjj=sitesel;
kkk=scensel;
ttt=timesel;

for sss=1:length(subcomp)
    u=squeeze(sum(sampsloccomponents{jjj,kkk}(:,subcomp{sss},ttt),2));
    for www=1:size(quantlevs,1)
        ht=.25/www;
        qq=quantile(u/10,quantlevs(www,:));
        patch(qq([1 2 2 1]),sss+[ht ht -ht -ht],'b'); hold on;
    end
    plot(median(u/10)*[1 1],sss+[.25 -.25],'Color',[.8 .8 .8]);
end
xlabel('cm');
set(gca,'ytick',1:length(subcomp),'yticklabel',labls,'ydir','reverse');