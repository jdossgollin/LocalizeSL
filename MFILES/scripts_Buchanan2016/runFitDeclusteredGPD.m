% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Apr 22 10:50:57 EDT 2016

rootpath='~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootpath,'MFILES'));

IFILESpath=fullfile(rootpath,'IFILES/extremes');
NOAAidfile = fullfile(IFILESpath,'USStationsLinearSeaLevelTrends.csv');
PSMSLidfile = fullfile(IFILESpath,'PSMSLIDs.csv');
datDir = fullfile(IFILESpath,'declustered');

% open file with list of tide gauges and 99th percentile thresholds

q99dat=importdata(fullfile(datDir,'q99.tsv'));
tgids=q99dat(:,1);
q99s=q99dat(:,2);

% open NOAA file so we can label sites
NOAAdat=importdata(NOAAidfile);
NOAAidlist=NOAAdat.textdata(2:end,1);
NOAAnames=NOAAdat.textdata(2:end,2);
NOAAlat=NOAAdat.data(:,9);
NOAAlong=NOAAdat.data(:,10);

% idlist is a bunch of strings - convert to numbers
for qqq=1:length(NOAAidlist)
    NOAAids(qqq)=str2num(NOAAidlist{qqq});
end

% open PSMSL file to match to NOAA
dat1=importdata(PSMSLidfile);
PSMSLsitenames=dat1.textdata(2:end,1);
PSMSLids=dat1.data(:,1);
PSMSLlat=dat1.data(:,2);
PSMSLlong=dat1.data(:,3);

distance_tolerance = .05;

clear NOAAsiteid NOAA_PSMSLid NOAA_PSMSLname mdist;
for qqq=1:length(NOAAnames)
    [mdist(qqq),mi]=min((NOAAlat(qqq)-PSMSLlat).^2+(NOAAlong(qqq)-PSMSLlong).^2);
    if mdist(qqq)<distance_tolerance
        NOAA_PSMSLid(qqq)=PSMSLids(mi);
        NOAA_PSMSLname{qqq}=PSMSLsitenames{mi};
    else
        NOAA_PSMSLid(qqq)=NaN;
        NOAA_PSMSLname{qqq}='not found';
    end
end

% run through list and fit GPD to exceedances of threshold

clear sitenames lambdas parmfit parmfit_gumbel N10 N01;
for qqq=1:length(tgids)
    wfile=fullfile(datDir,['maxtofit.dclist.' num2str(tgids(qqq)) '_xdat.dc.tsv']);
    
    % make sure file exists
    
    if exist(wfile,'file')  
        
        disp(tgids(qqq));
        
        % open file
        
        data=importdata(wfile);
        threshold=q99s(qqq);
        
        % identify exceedances
        data2=data(data>threshold)-threshold;
        
        % calculate Poission and GPD parameters
        lambdas(qqq)=length(data2)/(length(data)/365.25);
        parmfit(qqq,:) = gpfit(data2);
        parmfit_gumbel(qqq,:) = gpfit_forceshape(data2,[],[],0);
        
        
        % calculate 10% and 1% AEP heights
        if parmfit(qqq,1)==0
            z0 = @(N0) -parmfit(qqq,2)*log(N0)/log(lambdas(qqq)) + threshold;
        else
            z0 = @(N0) ((N0/lambdas(qqq))^(-parmfit(qqq,1))-1) * parmfit(qqq,2)/parmfit(qqq,1) + threshold;
        end
        N10(qqq)=z0(0.1);
        N01(qqq)=z0(0.01);
 
        % calculate 10% and 1% AEP heights for gumbel
        z0_gumbel = @(N0) -parmfit_gumbel(qqq,2)*log(N0)/log(lambdas(qqq)) + threshold;
        N10_gumbel(qqq)=z0_gumbel(0.1);
        N01_gumbel(qqq)=z0_gumbel(0.01);

    else
        % oops, couldn't find data file!
        
        lambdas(qqq)=NaN;
        parmfit(qqq,:)=[NaN NaN];
        N10(qqq)=NaN;
        N01(qqq)=NaN;
        
    end
end

% output table

fid=fopen('GPDfits.tsv','w');
fprintf(fid,'Site\tSite (PSMSL Name)\tPSMSL ID\tNOAA Station ID\tlambda\tu\tshape\tscale\tAEP0.1\tAEP0.010\n');
for qqq=1:length(tgids)
    sub=find(NOAAids==tgids(qqq));
    if (~isnan(lambdas(qqq)))&&(length(sub)==1)
        fprintf(fid,NOAAnames{sub});
        fprintf(fid,['\t' NOAA_PSMSLname{sub}]);
        fprintf(fid,'\t%0.0f',NOAA_PSMSLid(sub));
        fprintf(fid,'\t%0.0f',tgids(qqq));
        fprintf(fid,'\t%0.4f',[lambdas(qqq) q99s(qqq) parmfit(qqq,:) N10(qqq) N01(qqq)]);
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('GPDfits_gumbel.tsv','w');
fprintf(fid,'Site\tSite (PSMSL Name)\tPSMSL ID\tNOAA Station ID\tlambda\tu\tshape\tscale\tAEP0.1\tAEP0.010\n');
for qqq=1:length(tgids)
    sub=find(NOAAids==tgids(qqq));
    if (~isnan(lambdas(qqq)))&&(length(sub)==1)
        fprintf(fid,NOAAnames{sub});
        fprintf(fid,['\t' NOAA_PSMSLname{sub}]);
        fprintf(fid,'\t%0.0f',NOAA_PSMSLid(sub));
        fprintf(fid,'\t%0.0f',tgids(qqq));
        fprintf(fid,'\t%0.4f',[lambdas(qqq) q99s(qqq)]);
        fprintf(fid,'\t%0.4f',[parmfit_gumbel(qqq,:) N10_gumbel(qqq) N01_gumbel(qqq)]);
        fprintf(fid,'\n');
    end
end
fclose(fid);

% plot and table with uncertainty

% load data
clear acov;
fid=fopen('GPDfits_withunc.tsv','w');
fprintf(fid,'Site\tSite (PSMSL Name)\tPSMSL ID\tNOAA Station ID\tlambda\tu\tshape\t5\t95\tscale\t5\t95\tAEP1.0\t5\t95\tAEP0.1\t5\t95\tAEP0.01\t5\t95\tAEP0.001\t5\t95\tAEP 1.0 expected\tAEP0.1 expected\tAEP0.01 expected\tAEP0.001 expected\tVshape\tVscale\tVscaleshape\n');
for qqq=1:length(tgids)
    subNOAA=find(NOAAids==tgids(qqq));
    if (~isnan(lambdas(qqq)))&&(length(subNOAA)==1)
        disp(tgids(qqq));

        wfile=fullfile(datDir,['maxtofit.dclist.' num2str(tgids(qqq)) '_xdat.dc.tsv']);

        data=importdata(wfile);
        threshold=q99s(qqq);
        data2=data(data>threshold)-threshold;

        testht=0:.01:5;
        % construct empirical distribution
        ct=sum(bsxfun(@gt,data,testht))/(length(data)/365.25);
        sub=find(diff(ct)<0);;

        % calculate curve for best fit values
        y=GPDLogNExceedances(testht,lambdas(qqq),parmfit(qqq,1),parmfit(qqq,2),-q99s(qqq));

        % NOW HERE WE CALCULATE THE CONFIDENCE INTERVALS

        [nlogl,acov(:,:,qqq)]=gplike(parmfit(qqq,:),data2);
	% gplike calculates the negative log likelihood of the parameters parmfit(qqq,:) for data 2
	% acov is the covariance of the parameters, from which we will now sample a lot of pairs

        % sample lots of pairs
        Nsamps=1000; quantlevs=[.05 .5 .95];
        sampparam=lhsnorm(parmfit(qqq,:),acov(:,:,qqq),Nsamps); % we’re sampling with a latin hypercube (lhsnorm) instead of randomly (mvnrnd) for efficiency
        sampparam(:,2)=max(.0001,sampparam(:,2));

        % now calculate the curve for each of the parameter samples
        clear ysamps;
        for iii=1:size(sampparam,1)
            ysamps(iii,:)=GPDLogNExceedances(testht,lambdas(qqq),sampparam(iii,1),sampparam(iii,2),-q99s(qqq));
        end
        
        qsampparam=quantile(sampparam,quantlevs);
        yq=quantile(ysamps,quantlevs);

        % calculate uncertainty on 10-year and 100-year flood levels
        sampz0 = @(N0) ((N0/lambdas(qqq)).^(-sampparam(:,1))-1) .* sampparam(:,2)./sampparam(:,1) + threshold;
        qlevel1 = quantile(sampz0(1),quantlevs);
        qlevel10 = quantile(sampz0(.1),quantlevs);
        qlevel100 = quantile(sampz0(.01),quantlevs);
        qlevel500 = quantile(sampz0(.002),quantlevs);
        subd=find(diff(mean(exp(ysamps)))<0);
        expect1 = interp1(mean(exp(ysamps(:,subd))),testht(subd),1)+threshold;
        expect10 = interp1(mean(exp(ysamps(:,subd))),testht(subd),.1)+threshold;
        expect100 = interp1(mean(exp(ysamps(:,subd))),testht(subd),.01)+threshold;
        expect500 = interp1(mean(exp(ysamps(:,subd))),testht(subd),.002)+threshold;
        
        ymean=log(mean(exp(ysamps)));
        
        % plot empirical and fit
        clf;
        plot(testht(sub),ct(sub),'bs');
        hold on;
        plot(testht+threshold,exp(y),'r');
        set(gca,'ysc','log');
        ylim([5e-3 1]);
        plot(testht+threshold,exp(yq),'r--');
        title(NOAAnames{subNOAA});
        xlabel('Flood level (m)');
      
        % added expected curve
        plot(testht+threshold,exp(ymean),'k','linew',2);
     
        fprintf(fid,NOAAnames{subNOAA});
        fprintf(fid,['\t' NOAA_PSMSLname{subNOAA}]);
        fprintf(fid,'\t%0.0f',NOAA_PSMSLid(subNOAA));
        fprintf(fid,'\t%0.0f',tgids(qqq));
        fprintf(fid,'\t%0.4f',[lambdas(qqq) q99s(qqq)]);
        fprintf(fid,'\t%0.4f',[parmfit(qqq,1) qsampparam([1 3],1)']);
        fprintf(fid,'\t%0.4f',[parmfit(qqq,2) qsampparam([1 3],2)']);
        fprintf(fid,'\t%0.4f',[qlevel1([2 1 3])]);
        fprintf(fid,'\t%0.4f',[qlevel10([2 1 3])]);
        fprintf(fid,'\t%0.4f',[qlevel100([2 1 3])]);
        fprintf(fid,'\t%0.4f',[qlevel500([2 1 3])]);
        fprintf(fid,'\t%0.4f',[expect1 expect10 expect100 expect500])
        %fprintf(fid,'\t%0.4f',[parmfitCE(qqq,1) parmfitCE(qqq,2) CEz10 CEz100]);
        fprintf(fid,'\t%0.4g',[acov(1,1,qqq) acov(2,2,qqq) acov(1,2,qqq)]);
        %fprintf(fid,'\t%0.4f',[parmfit_gumbel(qqq,2) N10_gumbel(qqq) N01_gumbel(qqq)]);
        fprintf(fid,'\n');
        
    end
end
fclose(fid);
