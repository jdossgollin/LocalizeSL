function ConditionalDistributionsPlotSeaLevelComponents(p,condtargs,projLOC,projLOChi,projLOClo,projLOCcomp,projLOCcomphi,projLOCcomplo,colsCOMP,colsCOMPlab,filesuffix,cmap)

% ConditionalDistributionsPlotSeaLevelComponents(p,condtargs,projLOC,projLOChi,projLOClo,projLOCcomp,projLOCcomphi,projLOCcomplo,colsCOMP,colsCOMPlab,filesuffix)
%
% Generate output plots and tables for conditional scenarios.
%
% INPUT
% -----
% p: core sea-level structure 
% condtargs: target heights (in mm) upon which conditioned
% projLOC: median LSL projection for each scenario
% projLOChi: high LSL projection for each scenario
% projLOClo: low LSL projection for each scenario
% projLOCcomp: median LSL projection for each scenario
% projLOChicomp: high LSL projection for each scenario
% projLOClocomp: low LSL projection for each scenario
% colsCOMP: columns of core files used for compribution breakdown
% colsCOMPlab: labels for compribution breakdown
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Jan 09 13:46:51 EST 2017
%

defval('filesuffix','');
defval('cmap','parula');

disp('Local Scenarios component table');
fid=fopen(['LocalScenariosComponents' filesuffix '.tsv'],'w');

fprintf(fid,'Local scenarios (cm)\n');
today=date;
fprintf(fid,['Produced by Robert Kopp on ' today '\n\n']);
fprintf(fid,'Site\tID\tLatitude\tLongitude\tScenario\tComponent\t2000');
fprintf(fid,'\t%0.0f',p.targyears);
fprintf(fid,'\n');

for www=1:length(p.targregions)
    selectedSite=p.targregions(www);
    % generate local samples

    for qqq=1:size(projLOCcomp,1)
        fprintf(fid,p.targregionnames{www});
        fprintf(fid,'\t%0.0f',p.targregions(www));
        fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
        fprintf(fid,'\t%0.1f - MED',condtargs(1,qqq)/1000);
        fprintf(fid,'\tTotal');
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projLOC(qqq,:,www)/10);
        fprintf(fid,'\n');
        
        for ttt=1:length(colsCOMP)
            fprintf(fid,p.targregionnames{www});
            fprintf(fid,'\t%0.0f',p.targregions(www));
            fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
            fprintf(fid,'\t%0.1f - MED',condtargs(1,qqq)/1000);
            fprintf(fid,['\t' colsCOMPlab{ttt}]);
            fprintf(fid,'\t0');
            fprintf(fid,'\t%0.1f',projLOCcomp(qqq,:,ttt,www)/10);
            fprintf(fid,'\n');
        end
        
        fprintf(fid,p.targregionnames{www});
        fprintf(fid,'\t%0.0f',p.targregions(www));
        fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
        fprintf(fid,'\t%0.1f - LOW',condtargs(1,qqq)/1000);
        fprintf(fid,'\tTotal');
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projLOClo(qqq,:,www)/10);
        fprintf(fid,'\n');
        
        for ttt=1:length(colsCOMP)
            fprintf(fid,p.targregionnames{www});
            fprintf(fid,'\t%0.0f',p.targregions(www));
            fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
            fprintf(fid,'\t%0.1f - LOW',condtargs(1,qqq)/1000);
            fprintf(fid,['\t' colsCOMPlab{ttt}]);
            fprintf(fid,'\t0');
            fprintf(fid,'\t%0.1f',projLOCcomplo(qqq,:,ttt,www)/10);
            fprintf(fid,'\n');
        end
        
        fprintf(fid,p.targregionnames{www});
        fprintf(fid,'\t%0.0f',p.targregions(www));
        fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
        fprintf(fid,'\t%0.1f - HIGH',condtargs(1,qqq)/1000);
        fprintf(fid,'\tTotal');
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projLOChi(qqq,:,www)/10);
        fprintf(fid,'\n');
        
        for ttt=1:length(colsCOMP)
            fprintf(fid,p.targregionnames{www});
            fprintf(fid,'\t%0.0f',p.targregions(www));
            fprintf(fid,'\t%0.2f',p.targsitecoords(www,:));
            fprintf(fid,'\t%0.1f - HI',condtargs(1,qqq)/1000);
            fprintf(fid,['\t' colsCOMPlab{ttt}]);
            fprintf(fid,'\t0');
            fprintf(fid,'\t%0.1f',projLOCcomphi(qqq,:,ttt,www)/10);
            fprintf(fid,'\n');
        end

        

    end
    
end
fclose(fid);

%%%

doyear=2100;
dot=find(p.targyears==doyear);
for ttt=1:length(colsCOMP)
    curlab=colsCOMPlab{ttt};
    for qqq=1:size(projLOC,1)
        disp(sprintf([curlab ' Component for %0.0f cm Scenario - %0.0f'],[condtargs(1,qqq)/10 doyear]));
        clf;
        %worldmap('North America');
        worldmap([-20 73],[120 -25]);
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        stateColor=[.8 .8 .8];
        geoshow('landareas.shp','facecolor',stateColor);
        
        u=(squeeze(projLOCcomp(qqq,dot,ttt,:)))/10;
        ulo=(squeeze(projLOCcomplo(qqq,dot,ttt,:)))/10;
        uhi=(squeeze(projLOCcomphi(qqq,dot,ttt,:)))/10;

        scatterm(p.targsitecoords(:,1),p.targsitecoords(:,2),10,u,'s','filled');

        pos0=get(gca,'position');
        hcb=colorbar('SouthOutside');

        ht=title(sprintf([curlab ' Component for %0.0f cm Scenario - %0.0f'],[condtargs(1,qqq)/10 doyear]));
        colormap(cmap);
        setm(gca,'mlinevisible','off','grid','off')
        crange=quantile([u(:)' ulo(:)' uhi(:)'],[.01 .99]);
        caxis(crange);
        pdfwrite(['LocalScenarioComponent-' curlab num2str(qqq) filesuffix]);    
        
        clf;
        %worldmap('North America');
        worldmap([-20 73],[120 -25]);
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        stateColor=[.8 .8 .8];
        geoshow('landareas.shp','facecolor',stateColor);
        scatterm(p.targsitecoords(:,1),p.targsitecoords(:,2),10,ulo,'s','filled');

        pos0=get(gca,'position');
        hcb=colorbar('SouthOutside');

        ht=title(sprintf([curlab ' Component for %0.0f cm Scenario - %0.0f (Low)'],[condtargs(1,qqq)/10 doyear]));
        colormap(cmap);
        setm(gca,'mlinevisible','off','grid','off')
        caxis(crange);
        pdfwrite(['LocalScenarioComponent-' curlab '-Low' num2str(qqq) filesuffix]);    
        
        clf;
        %worldmap('North America');
        worldmap([-20 73],[120 -25]);
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        stateColor=[.8 .8 .8];
        geoshow('landareas.shp','facecolor',stateColor);
        scatterm(p.targsitecoords(:,1),p.targsitecoords(:,2),10,uhi,'s','filled');

        pos0=get(gca,'position');
        hcb=colorbar('SouthOutside');

        ht=title(sprintf([curlab ' Component for %0.0f cm Scenario - %0.0f (High)'],[condtargs(1,qqq)/10 doyear]));
        colormap(cmap);
        setm(gca,'mlinevisible','off','grid','off')
        caxis(crange);
        pdfwrite(['LocalScenarioComponent-' curlab '-Hi' num2str(qqq) filesuffix]);    
        
    end
end

%%%%%

% background rate


clf;
%worldmap('North America');
worldmap([-20 73],[120 -25]);
setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
ax=gca;
stateColor=[.8 .8 .8];
geoshow('landareas.shp','facecolor',stateColor);
scatterm(p.targsitecoords(:,1),p.targsitecoords(:,2),10,p.rateprojs,'s','filled');
pos0=get(gca,'position');
hcb=colorbar('SouthOutside');
ht=title('Background rate (mm/yr)');
colormap(cmap);
setm(gca,'mlinevisible','off','grid','off')
pdfwrite(['LocalScenarioComponent-Bkgd' filesuffix]);    

clf;
%worldmap('North America');
worldmap([-20 73],[120 -25]);
setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
ax=gca;
stateColor=[.8 .8 .8];
geoshow('landareas.shp','facecolor',stateColor);
scatterm(p.targsitecoords(:,1),p.targsitecoords(:,2),10,p.rateprojssd,'s','filled');
pos0=get(gca,'position');
hcb=colorbar('SouthOutside');
ht=title('Background rate - standard deviation (mm/yr)');
colormap(cmap);
setm(gca,'mlinevisible','off','grid','off')
pdfwrite(['LocalScenarioComponent-Bkgd-Std' filesuffix]);    
