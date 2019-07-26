% 2019-07-25 19:59:55 -0400

rootdir='~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootdir,'MFILES'));

% import SEJ data

baseline=[275 0 200 75]/360; % in mm/yr 

dat=importdata('../USUK_2100H.csv');
data=bsxfun(@plus,dat.data,baseline*100);
[s,si]=sort(data(:,1));
data=data(si,:);

dat=importdata('../USUK_2100L.csv');
data2=bsxfun(@plus,dat.data,baseline*100);
[s,si]=sort(data2(:,1));
data2=data2(si,:);

dat=importdata('../USUK_2200H.csv');
dataC=bsxfun(@plus,dat.data,baseline*200);
[s,si]=sort(dataC(:,1));
dataC=dataC(si,:);

dat=importdata('../USUK_2200L.csv');
dataC2=bsxfun(@plus,dat.data,baseline*200);
[s,si]=sort(dataC2(:,1));
dataC2=dataC2(si,:);

dat=importdata('../USUK_2300H.csv');
dataB=bsxfun(@plus,dat.data,baseline*300);
[s,si]=sort(dataB(:,1));
dataB=dataB(si,:);

dat=importdata('../USUK_2300L.csv');
dataB2=bsxfun(@plus,dat.data,baseline*300);
[s,si]=sort(dataB2(:,1));
dataB2=dataB2(si,:);

dat=importdata('../USUK_2050H.csv');
dataA=bsxfun(@plus,dat.data,baseline*50);
[s,si]=sort(dataA(:,1));
dataA=dataA(si,:);

dat=importdata('../USUK_2050L.csv');
dataA2=bsxfun(@plus,dat.data,baseline*50);
[s,si]=sort(dataA2(:,1));
dataA2=dataA2(si,:);

%%%%

% derive correlation matrix from DP16

targyears2=[2050 2100 2200 2300];
corefileA=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat'));
ISdp16h=squeeze(sum(corefileA.samps(:,[corefileA.colAIS corefileA.colGIS],:,1),2));
ISdp16m=squeeze(sum(corefileA.samps(:,[corefileA.colAIS corefileA.colGIS],:,3),2));
ISdp16l=squeeze(sum(corefileA.samps(:,[corefileA.colAIS corefileA.colGIS],:,4),2));
subt=find(ismember(corefileA.targyears,targyears2));
ISdp16=[ISdp16h ; ISdp16m ; ISdp16l];
ISdp16=ISdp16(:,subt);
CV2=corr(ISdp16);
clear corefileA;

samps=normcdf(mvnrnd([0 0 0 0],CV2,10000));
sampsi=ceil(samps*size(data,1));

clear dataH dataL dataBaseline;
dataH(:,:,1) = dataA(sampsi(:,1),:);
dataH(:,:,2) = data(sampsi(:,2),:);
dataH(:,:,3) = dataC(sampsi(:,3),:);
dataH(:,:,4) = dataB(sampsi(:,4),:);

dataL(:,:,1) = dataA2(sampsi(:,1),:);
dataL(:,:,2) = data2(sampsi(:,2),:);
dataL(:,:,3) = dataC2(sampsi(:,3),:);
dataL(:,:,4) = dataB2(sampsi(:,4),:);

dataBaseline(:,:,1) = baseline*50;
dataBaseline(:,:,2) = baseline*100;
dataBaseline(:,:,3) = baseline*200;
dataBaseline(:,:,4) = baseline*300;

dataH = permute(dataH,[1 3 2]);
dataL = permute(dataL,[1 3 2]);
dataBaseline = permute(dataBaseline,[1 3 2]);

%%%

% check correlations

dataIS=[dataH(:,:,2)+dataH(:,:,3)+dataH(:,:,4) ; dataL(:,:,2)+ dataL(:,:,3) + dataL(:,:,4)];
for ttt=1:(length(targyears2)-1)
    for sss=(ttt+1):length(targyears2)
        clf;
        plot(dataIS(:,ttt),dataIS(:,sss),'.');
        hold on;
        plot(ISdp16(1:100:end,ttt),ISdp16(1:100:end,sss),'rd');
        xlabel(num2str(targyears2(ttt)));
        ylabel(num2str(targyears2(sss)));
        legend('B19','K17+DP16');
        pdfwrite(['DP16B19_' num2str(targyears2(sss)) 'vs' num2str(targyears2(ttt))]);
    end
end

corr(dataIS)

%%%

% this is difference between 1990-2000 rate and 2000-2010 rate, from Mouginot et al 2019 and Rignot et al. 2019
baselineanom1990=[-262 -73 -145 -44]/360;

% this is difference between 2010-201x rate and 2000-2010 rate, from Mouginot et al 2019 and Rignot et al. 2019
baselineanom2013=[185 -31 100 116]/360;
baselinelev2013=baselineanom2013*3;
slope2000=baselineanom1990/2;

t=[0 10 13.5 50 100 200 300];
tex=10:10:300;
clear dataHex dataLex;
for sss=1:size(dataH,3)
    debaselined = [zeros(size(dataH,1),1) zeros(size(dataH,1),1) ones(size(dataH,1),1)*baselinelev2013(sss) dataH(:,:,sss)-dataBaseline(:,:,sss)];
    slope=(debaselined(:,[end-1 end])-debaselined(:,[end-2 end-1]))/(t(end)-t(end-1));
    enslope=slope(:,2)+.5*(slope(:,2)-slope(:,1));
    debaselined2 = [ones(size(dataH,1),1)*baselineanom1990(sss) debaselined enslope];
    dataHex(:,:,sss)=spline(t,debaselined2,tex)+bsxfun(@times,baseline(sss),tex);

    debaselined = [zeros(size(dataL,1),1) zeros(size(dataL,1),1) ones(size(dataL,1),1)*baselinelev2013(sss) dataL(:,:,sss)-dataBaseline(:,:,sss)];
    slope=(debaselined(:,[end-1 end])-debaselined(:,[end-2 end-1]))/(t(end)-t(end-1));
    enslope=slope(:,2)+.5*(slope(:,2)-slope(:,1));
    debaselined2 = [ones(size(dataH,1),1)*baselineanom1990(sss) debaselined enslope];
    dataLex(:,:,sss)=spline(t,debaselined2,tex)+bsxfun(@times,baseline(sss),tex);
end

%%%%


labl={'','EAIS','GrIS','WAIS'};
quantlevs=[.05 .17 .5 .83 .95];
for sss=1:size(dataH,3)
    clf;
    subplot(2,1,1);
    plot(tex,quantile(dataHex(:,:,sss),quantlevs));
    title([labl{sss} ' - H']); ylabel('mm');
    subplot(2,1,2);
    plot(.5*tex(1:end-1)+.5*tex(2:end),quantile(diff(dataHex(:,:,sss),[],2)/10,quantlevs));
    ylabel('mm/yr');
    pdfwrite(['fitH' num2str(sss)]);

    fid=fopen(['fitH' num2str(sss) '.tsv'],'w');
    fprintf(fid,[labl{sss} ' - H\n']);
    fprintf(fid,'year\tmm\n');
    fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
    for ttt=1:length(tex)
        fprintf(fid,'%0.0f',2000+tex(ttt));
        fprintf(fid,'\t%0.0f',quantile(dataHex(:,ttt,sss),quantlevs));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'year\tmm/yr\n');
    fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
    for ttt=1:length(tex)-1
        fprintf(fid,'%0.0f',2000+tex(ttt));
        fprintf(fid,'\t%0.1f',quantile(diff(dataHex(:,ttt+[0 1],sss),[],2)/10,quantlevs));
        fprintf(fid,'\n');
    end
    fclose(fid);

    clf;
    subplot(2,1,1);
    plot(tex,quantile(dataLex(:,:,sss),quantlevs));
    title([labl{sss} ' - L']); ylabel('mm');
    subplot(2,1,2);
    plot(.5*tex(1:end-1)+.5*tex(2:end),quantile(diff(dataLex(:,:,sss),[],2)/10,quantlevs));
    ylabel('mm/yr');

    pdfwrite(['fitL' num2str(sss)]);

    fid=fopen(['fitL' num2str(sss) '.tsv'],'w');
    fprintf(fid,[labl{sss} ' - H\n']);
    fprintf(fid,'year\tmm\n');
    fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
    for ttt=1:length(tex)
        fprintf(fid,'%0.0f',2000+tex(ttt));
        fprintf(fid,'\t%0.0f',quantile(dataLex(:,ttt,sss),quantlevs));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'year\tmm/yr\n');
    fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
    for ttt=1:length(tex)-1
        fprintf(fid,'%0.0f',2000+tex(ttt));
        fprintf(fid,'\t%0.1f',quantile(diff(dataLex(:,ttt+[0 1],sss),[],2)/10,quantlevs));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

clf;
subplot(2,1,1);
plot(tex,quantile(sum(dataHex(:,:,2:4),3),quantlevs));
title('total - H'); ylabel('mm');
subplot(2,1,2);
plot(.5*tex(1:end-1)+.5*tex(2:end),quantile(diff(sum(dataHex(:,:,2:4),3),[],2)/10,quantlevs));
ylabel('mm/yr');

pdfwrite(['fitH_total']);

fid=fopen(['fitH_total.tsv'],'w');
fprintf(fid,['total - H\n']);
fprintf(fid,'year\tmm\n');
fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
for ttt=1:length(tex)
    fprintf(fid,'%0.0f',2000+tex(ttt));
    fprintf(fid,'\t%0.0f',quantile(sum(dataHex(:,ttt,2:4),3),quantlevs));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'year\tmm/yr\n');
fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
for ttt=1:length(tex)-1
    fprintf(fid,'%0.0f',2000+tex(ttt));
    fprintf(fid,'\t%0.1f',quantile(diff(sum(dataHex(:,ttt+[0 1],2:4),3),[],2)/10,quantlevs));
    fprintf(fid,'\n');
end
fclose(fid);

clf
subplot(2,1,1);
plot(tex,quantile(sum(dataLex(:,:,2:4),3),quantlevs));
title('total - L'); ylabel('mm');
subplot(2,1,2);
plot(.5*tex(1:end-1)+.5*tex(2:end),quantile(diff(sum(dataLex(:,:,2:4),3),[],2)/10,quantlevs));
ylabel('mm/yr');

pdfwrite(['fitL_total']);

fid=fopen(['fitL_total.tsv'],'w');
fprintf(fid,[labl{sss} ' - H\n']);
fprintf(fid,'year\tmm\n');
fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
for ttt=1:length(tex)
    fprintf(fid,'%0.0f',2000+tex(ttt));
    fprintf(fid,'\t%0.0f',quantile(sum(dataLex(:,ttt,2:4),3),quantlevs));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'year\tmm/yr\n');
fprintf(fid,'\t%0.2f',quantlevs); fprintf(fid,'\n');
for ttt=1:length(tex)-1
    fprintf(fid,'%0.0f',2000+tex(ttt));
    fprintf(fid,'\t%0.1f',quantile(diff(sum(dataLex(:,ttt+[0 1],2:4),3),[],2)/10,quantlevs));
    fprintf(fid,'\n');
end
fclose(fid);

%%%

corefile=load(fullfile(rootdir,'IFILES/SLRProjections170113GRIDDEDcore.mat'));
corefile.parent='SLRProjections170113GRIDDEDcore';
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile);

corefile2=corefile;
subt=1:length(targyears);
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
subdat=si(round(linspace(1,size(dataHex,1),size(corefile.samps,1))));

corefile2.samps(:,cols.colAIS(2),:,1)=dataHex(subdat,:,2);
corefile2.samps(:,cols.colAIS(1),:,1)=dataHex(subdat,:,4);
corefile2.samps(:,cols.colGIS,:,1)=dataHex(subdat,:,3);

corefileH=corefile2;

%%%%

corefile=load(fullfile(rootdir,'IFILES/SLRProjections180124GRIDDEDcore_Tscens.mat'));
corefile.parent='SLRProjections180124GRIDDEDcore_Tscens';
[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile,1);

corefile2=corefile;
subt=1:length(targyears);
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

[s,si]=sort(data(:,1));
subdat=si(round(linspace(1,size(dataLex,1),size(corefile.samps,1))));

corefile2.samps(:,cols.colAIS(2),:,1)=dataLex(subdat,:,2);
corefile2.samps(:,cols.colAIS(1),:,1)=dataLex(subdat,:,4);
corefile2.samps(:,cols.colGIS,:,1)=dataLex(subdat,:,3);

corefileL=corefile2;
clear corefile corefile2;

save SLRProjections190725core_SEJ_full corefileH corefileL

%

% clear;
% rootdir='~/Dropbox/Code/LocalizeSL';
% corefiles=load(fullfile(rootdir,'IFILES/SLRProjections190301core_SEJ.mat'));

% selectedSite=0;
% [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefiles.corefileH,1);
% [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefiles.corefileL,1);
