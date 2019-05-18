% 2019-02-10 19:04:17 -0500

rootdir='~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootdir,'MFILES'));

baseline=[275 0 200 75]/360; % in mm/yr 

dat=importdata('../USUK_2100H.csv');
data=bsxfun(@plus,dat.data,baseline*100);

dat=importdata('../USUK_2100L.csv');
data2=bsxfun(@plus,dat.data,baseline*100);

dat=importdata('../USUK_2300H.csv');
dataB=bsxfun(@plus,dat.data,baseline*300);

dat=importdata('../USUK_2300L.csv');
dataB2=bsxfun(@plus,dat.data,baseline*300);

dat=importdata('../USUK_2050H.csv');
dataA=bsxfun(@plus,dat.data,baseline*50);

dat=importdata('../USUK_2050L.csv');
dataA2=bsxfun(@plus,dat.data,baseline*50);

corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore.mat'));
corefile.parent='SLRProjections170113GRIDDEDcore';
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile);

corefile2=corefile;
subt=find(ismember(targyears,[2050 2100 2300]));
subscens=[1];
corefile2.OceanDynMean=corefile2.OceanDynMean(:,:,subscens);
corefile2.OceanDynN=corefile2.OceanDynN(:,:,subscens);
corefile2.OceanDynStd=corefile2.OceanDynStd(:,:,subscens);
corefile2.OceanDynTECorr=corefile2.OceanDynTECorr(:,:,subscens);
corefile2.ThermExpMean=corefile2.ThermExpMean(:,subscens);
corefile2.ThermExpN=corefile2.ThermExpN(:,subscens);
corefile2.ThermExpStd=corefile2.ThermExpStd(:,subscens);
corefile2.samps=corefile2.samps(:,:,subt,subscens);
corefile2.targyears=corefile2.targyears(subt);
corefile2.scens={'rcp85+H'};

[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile2,1);

[s,si]=sort(data(:,1));
subdat=si(round(linspace(1,size(data,1),size(corefile.samps,1))));

subt=find(targyears==2100);
corefile2.samps(:,cols.colAIS(2),subt,1)=data(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=data(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=data(subdat,3);

subt=find(targyears==2050);
corefile2.samps(:,cols.colAIS(2),subt,1)=dataA(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=dataA(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=dataA(subdat,3);


subt=find(targyears==2300);
corefile2.samps(:,cols.colAIS(2),subt,1)=dataB(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=dataB(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=dataB(subdat,3);

corefileH=corefile2;

%%%%

corefile=load(fullfile(rootdir,'IFILES/SLRProjections180124GRIDDEDcore_Tscens.mat'));
corefile.parent='SLRProjections180124GRIDDEDcore_Tscens';
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile,1);

corefile2=corefile;
subt=find(ismember(targyears,[2050 2100 2300]));
subscens=[2];
corefile2.OceanDynMean=corefile2.OceanDynMean(:,:,subscens);
corefile2.OceanDynN=corefile2.OceanDynN(:,:,subscens);
corefile2.OceanDynStd=corefile2.OceanDynStd(:,:,subscens);
corefile2.OceanDynTECorr=corefile2.OceanDynTECorr(:,:,subscens);
corefile2.ThermExpMean=corefile2.ThermExpMean(:,subscens);
corefile2.ThermExpN=corefile2.ThermExpN(:,subscens);
corefile2.ThermExpStd=corefile2.ThermExpStd(:,subscens);
corefile2.samps=corefile2.samps(:,:,subt,subscens);
corefile2.targyears=corefile2.targyears(subt);
corefile2.scens={'2p0degree+L'};

[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile2,1);

[s,si]=sort(data2(:,1));
subdat=si(round(linspace(1,size(data2,1),size(corefile.samps,1))));

subt=find(targyears==2100);
corefile2.samps(:,cols.colAIS(2),subt,1)=data2(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=data2(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=data2(subdat,3);

subt=find(targyears==2050);
corefile2.samps(:,cols.colAIS(2),subt,1)=dataA2(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=dataA2(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=dataA2(subdat,3);


subt=find(targyears==2300);
corefile2.samps(:,cols.colAIS(2),subt,1)=dataB2(subdat,2);
corefile2.samps(:,cols.colAIS(1),subt,1)=dataB2(subdat,4);
corefile2.samps(:,cols.colGIS,subt,1)=dataB2(subdat,3);

corefileL=corefile2;
clear corefile corefile2;

save SLRProjections190301core_SEJ corefileH corefileL

%

clear;
rootdir='~/Dropbox/Code/LocalizeSL';
corefiles=load(fullfile(rootdir,'IFILES/SLRProjections190301core_SEJ.mat'));

selectedSite=0;
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefiles.corefileH,1);
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefiles.corefileL,1);
