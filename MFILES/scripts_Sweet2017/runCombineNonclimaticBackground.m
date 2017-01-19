% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Jan 11 15:07:11 EST 2017

scenimport=ImportScenarioOutput('../workdir-161101/LocalScenarios-NoBkgd.tsv');
scenimport2=ImportScenarioOutput('../workdir-170109-Low/LocalScenarios-NoBkgd.tsv');
bkgd=ImportBkgdRateOutput('bkgdrate.tsv');
scenmerge=scenimport;
targyears=scenmerge.yrs;

for www=1:length(scenmerge.scen)
    if scenmerge.scen(www)==0.3
        sub=find(scenimport2.psmslid==scenmerge.psmslid(www));
        sub=intersect(sub,find(strcmpi(scenmerge.lev{www},scenimport2.lev)));
        scenmerge.projs(www,:)=scenimport2.projs(sub,:);
    end
    scenmerge.appliedbkgd(www)=0;
end

scenwithbkgd=scenmerge;
for www=1:length(scenwithbkgd.scen)
    scenwithbkgd.appliedbkgd(www)=0;
    if scenimport.psmslid(www)>0
        sub=find(bkgd.psmslid==scenimport.psmslid(www));
        dobkgd=bkgd.bkgdrate(sub);
        dobkgd_1s=bkgd.bkgdrate_std(sub);
        scenwithbkgd.appliedbkgd(www)=dobkgd;
        scenwithbkgd.projs(www,:)=scenwithbkgd.projs(www,:)+(scenwithbkgd.yrs-scenwithbkgd.yrs(1))*.1*dobkgd;
        
        if strcmpi(scenimport.lev{www},'LOW')
            scenwithbkgd.projs(www,:)=scenwithbkgd.projs(www,:)-(scenwithbkgd.yrs-scenwithbkgd.yrs(1))*.1*dobkgd_1s;
            scenwithbkgd.appliedbkgd(www)=dobkgd-dobkgd_1s;
        elseif strcmpi(scenimport.lev{www},'HIGH');
            scenwithbkgd.projs(www,:)=scenwithbkgd.projs(www,:)+(scenwithbkgd.yrs-scenwithbkgd.yrs(1))*.1*dobkgd_1s;
            scenwithbkgd.appliedbkgd(www)=dobkgd+dobkgd_1s;
        end
        
    end
    
end

%%%%


%%%%

for sss=1:2
    if sss==1
        filesuffix='-Climatic';
        longlabel='Climatic with Uncorrelated Low';
        doscen=scenmerge;
    elseif sss==2
        filesuffix='-ClimaticPlusBackground';
        longlabel='Climatic with Uncorrelated Low and Background Added';
        doscen=scenwithbkgd;
    end
    

    fid=fopen(['LocalScenarios' filesuffix '.tsv'],'w');

    fprintf(fid,['Local scenarios (cm) - ' longlabel '\n']);
    today=date;
    fprintf(fid,['Produced by Robert Kopp on ' today '\n\n']);
    fprintf(fid,'Site\tID\tLatitude\tLongitude\tScenario\tBackground');
    fprintf(fid,'\t%0.0f',targyears);
    fprintf(fid,'\n');

    for www=1:length(doscen.sitelab)
        fprintf(fid,doscen.sitelab{www});
        fprintf(fid,'\t%0.0f',doscen.psmslid(www));
        fprintf(fid,'\t%0.2f',[doscen.lat(www) doscen.lon(www)]);
        fprintf(fid,['\t%0.1f - ' doscen.lev{www}],doscen.scen(www));
        fprintf(fid,'\t%0.2f',doscen.appliedbkgd(www));
        fprintf(fid,'\t%0.0f',doscen.projs(www,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

%%%%


doyear=2100;
dot=find(targyears==doyear);
crange=[-80 80];
cmap=brewermap(16,'RdYlBu');
cmap=cmap(end:-1:1,:);

for sss=1:2
    if sss==1
        filesuffix='-Climatic';
        longlabel='Climatic with Uncorrelated Low';
        doscen=scenmerge;
    elseif sss==2
        filesuffix='-ClimaticPlusBackground';
        longlabel='Climatic with Uncorrelated Low and Background Added';
        doscen=scenwithbkgd;
    end
    
    condtargs=unique(doscen.scen)*1000;
    condsets=unique(doscen.lev);

    for qqq=1:length(condtargs)
        for rrr=1:length(condsets)
            labsets=condsets{rrr};
            if strcmpi('MED',labsets)
                labsets='';
            else
                labsets=[' (' labsets ')'];
            end
            titl=(sprintf(['Adder for %0.0f cm Scenario - %0.0f' labsets],[condtargs(qqq)/10 doyear]));
            disp(titl);
            doset=find(strcmpi(condsets{rrr},doscen.lev));
            doset=intersect(doset,find(doscen.scen==condtargs(qqq)/1000));
            
            
            clf;
            %worldmap('North America');
            worldmap([-20 73],[120 -25]);

            setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
            ax=gca;
            stateColor=[.8 .8 .8];
            geoshow('landareas.shp','facecolor',stateColor);
            states=shaperead('usastatehi', 'UseGeoCoords', true);
            geoshow(states,'FaceColor','none');
            u=doscen.projs(doset,dot)-condtargs(qqq)/10;
            
            scatterm(doscen.lat(doset),doscen.lon(doset),10,u,'s','filled');

            pos0=get(gca,'position');
            hcb=colorbar('SouthOutside');

            ht=title(titl);
            colormap(cmap);
            setm(gca,'mlinevisible','off','grid','off')
            caxis(crange);
            pdfwrite(['LocalScenarioAdder' '-' condsets{rrr} '-' num2str(qqq) filesuffix]);               
        end
    end
end
