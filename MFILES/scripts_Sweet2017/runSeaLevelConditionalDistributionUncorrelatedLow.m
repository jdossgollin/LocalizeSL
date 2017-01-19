% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Jan 09 12:25:45 EST 2017

% this one does the uncorrelated low only

rootdir='~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootdir,'MFILES'));

savefilecore=fullfile(rootdir,'IFILES/SLRProjections161027GRIDDEDcore.mat');
p=load(savefilecore);

condtargyrs=[2100 2050 2030];
condtargs=[30 50 100 150 200 250 ;
           15 NaN NaN NaN NaN NaN ;
           9 NaN NaN NaN NaN NaN] * 10;
condtargwins=[20 20 20 50 50 150 ;
              20 20 20 20 20 20 ;
              10 10 10 10 10 10];

% turn off correlation of TE and OD
p.OceanDynTECorr=p.OceanDynTECorr*0;
condtargs=condtargs(:,1);

disp('Conditionalizing GSL...');
[projGSL,condsubscen]=ConditionalDistributionsGSL(p,condtargyrs,condtargs,condtargwins);

% only do condsubscen(1)

disp('Conditionalizing LSL...');
projLOC=ConditionalDistributionsLSL(p,condsubscen);

% now output

fid=fopen('conditionalScenariosN.tsv','w'); fprintf(fid,'scenario\tN\n');
for ppp=1:length(condsubscen)
    fprintf(fid,'%0.0f',ppp);
    fprintf(fid,'\t%0.0f',length(condsubscen{ppp}));
    fprintf(fid,'\n');
end
fclose(fid);

crange0=[-80 80];
cmap=brewermap(16,'RdYlBu');
cmap=cmap(end:-1:1,:);

subrateyrs=find(projGSL.targyearrates<2100);



filesuffix='';
ConditionalDistributionsPlotSeaLevel(p,condtargs,projGSL.proj,projGSL.projhi,projGSL.projlo,projLOC.projLOC,projLOC.projLOChi,projLOC.projLOClo,projGSL.targyearrates(subrateyrs),projGSL.projrate(:,subrateyrs),projGSL.projratehi(:,subrateyrs),projGSL.projratelo(:,subrateyrs),projLOC.projLOCrate(:,subrateyrs,:),projLOC.projLOCratehi(:,subrateyrs,:),projLOC.projLOCratelo(:,subrateyrs,:),filesuffix,crange0,cmap);

filesuffix='-NoBkgd';
ConditionalDistributionsPlotSeaLevel(p,condtargs,projGSL.proj,projGSL.projhi,projGSL.projlo,projLOC.projLOC0,projLOC.projLOC0hi,projLOC.projLOC0lo,projGSL.targyearrates(subrateyrs),projGSL.projrate(:,subrateyrs),projGSL.projratehi(:,subrateyrs),projGSL.projratelo(:,subrateyrs),projLOC.projLOC0rate(:,subrateyrs,:),projLOC.projLOC0ratehi(:,subrateyrs,:),projLOC.projLOC0ratelo(:,subrateyrs,:),filesuffix,crange0,cmap);

filesuffix='';

ConditionalDistributionsPlotGSLComponents(p,condtargs,projGSL.proj,projGSL.projhi,projGSL.projlo,projGSL.projCONT,projGSL.projCONThi,projGSL.projCONTlo,projGSL.colsCONT,projGSL.colsCONTlab);
   
ConditionalDistributionsPlotSeaLevelComponents(p,condtargs,projLOC.projLOC,projLOC.projLOChi,projLOC.projLOClo,projLOC.projLOCcomp,projLOC.projLOCcomphi,projLOC.projLOCcomplo,projLOC.colsCOMP,projLOC.colsCOMPlab,filesuffix);
