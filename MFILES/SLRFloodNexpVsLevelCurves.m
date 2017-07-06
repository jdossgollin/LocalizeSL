function [effcurve,testz,histcurve,histcurvesamps,effcurveESLR,effcurve999,integratecurve]=SLRFloodNexpVsLevelCurves(samps,targyears,threshold,scale,shape,lambda,sitelab,params)

% [effcurve,testz,histcurve,histcurvesamps,effcurveESLR,effcurve999,integratecurve]=
%    SLRFloodNexpVsLevelCurves(samps,targyears,threshold,scale,shape,lambda,
%    [sitelab],[params])
%
% Output a figure plotting expected number of floods vs flood height for
% historic GPD, shifted GPDs, and AADLLs.
%
% INPUTS:
%
%     samps: matrix of SLR samples (rows: samples, cols: time)
%     targyears: year identifiers for columns on samps
%     threshold: GPD threshold
%     scale: GPD scale factor
%     shape: GPD shape factor
%     lambda: mean number of exceedances per year 
%     sitelab: full name of site
%     params: an optional structure with fields setting several parameters
%         - startyear: start year for AADLL curves (default: 2020) 
%         - endyears: end years for AADLLs; instantaneous effective flood risk
%                     for these years are also plotted
%                     (default: [2050 2100])
%         - testz: vector of heights used for calculating the curve
%                  (default: [0:.01:10])
%         - doplot: generate plot (default: 1)
%
% OUTPUTS:
%
%     effcurve: expected number of floods of different heights (rows: years as
%               in targyears; cols: heights as specified in testz)
%     testz: heights for effcurve
%     histcurve: historical curve for testz (expected)
%     histcurvesamps: historical curve for testz (samples from GPD parameters)
%     effcurveESLR: expected number of floods under expected SLR
%     effcurve999: expected number of floods under 99.9th percentile SLR
%     integratecurve: function (curve,t1,t2) to integrate effcurve
%
% EXAMPLE:
%
%     sitelab='Boston';
%     selectedSite = 235; %PSMSL ID for Boston
%     threshold = 0.6214; % GPD threshold
%     scale = 0.1109; % GPD scale
%     shape = 0.0739; % GPD shape
%     lambda = 2.5699; % Poisson Lambda
%     
%     [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = ...
%         LocalizeStoredProjections(selectedSite,corefile,1);
%     samps=[zeros(size(sampslocrise{1,1},1),1) ...
%            sampslocrise{1,1}]/1000; % add base yr, convert to meters
%     samps=bsxfun(@min,samps,quantile(samps,.999));
%     % truncate samples viewed as physically implausible
%     targyears = [2000 targyears]; % add base year
%     
%     clf;
%     [effcurve,testz,histcurve,histcurvesamps,effcurveESLR,effcurve999,integratecurve]= ...
%         SLRFloodNexpVsLevelCurves(samps,targyears,threshold, ...
%                                   scale,shape,lambda,sitelab);     
%     pdfwrite([sitelab '_returncurves']);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-07-06 19:39:35 -0400

defval('sitelab',[]);
defval('testz',0:.01:10);
defval('startyear',2020);
defval('endyears',[2050 2100]);
defval('doplot',1);
defval('effcurve',[]);
defval('showuncertainty',0);
defval('historicaldata',[]);
defval('showESLR',1);
defval('showeffcurve',1);
defval('show999',1);
defval('showinteffcurve',1);
defval('historicalcolor','y');

if exist('params')
    parseFields(params)
end

ESLR = mean(samps);
SLR999=quantile(samps,.999);

if size(shape,1)>1
    logN = @(z) GPDELogNExceedances(z-threshold,lambda,shape,scale,-threshold); 
else
    logN = @(z) GPDLogNExceedances(z-threshold,lambda,shape,scale,-threshold); 
end

if length(effcurve)==0
    for ttt=1:length(targyears)
        effcurve(ttt,:) = real(mean(exp(logN(bsxfun(@minus,testz,samps(:,ttt)))),1));
    end
end

for ttt=1:length(targyears)
    effcurveESLR(ttt,:)=exp(logN(testz-ESLR(ttt)));
    effcurve999(ttt,:)=exp(logN(testz-SLR999(ttt)));
end

if size(shape,1)>1
    for iii=1:size(shape,1)
        logNs = @(z) GPDLogNExceedances(z-threshold,lambda,shape(iii),scale(iii),-threshold);
        histcurvesamps(iii,:) = real(exp(logNs(testz)));
    end
else
    histcurvesamps = effcurve(1,:);
end

histcurve=exp(logN(testz));
integratecurve=@(curve,t1,t2) trapz(t1:t2,interp1(targyears,curve,t1:t2))./(t2-t1);


%%%

if doplot

    hp1(1)=subplot(2,1,1); hold on;

    iii=1;
    if showuncertainty
        if size(shape,1)>1
            greyc=[.7 .7 .7];
            plot(testz,quantile(histcurvesamps,.5),'-','Color',greyc);
            hold on;
            plot(testz,quantile(histcurvesamps,.167),'--','Color',greyc);
            plot(testz,quantile(histcurvesamps,.833),'--','Color',greyc);
        end
    end
    
    if length(historicaldata)>0
        ct=sum(bsxfun(@gt,historicaldata,testz))/(length(historicaldata)/365.25);
        subct=find(diff(ct)<0);
        subct=intersect(subct,find(testz>threshold));
        plot(testz(subct),ct(subct),'s','Color',historicalcolor);
    end
    
    hl(iii)=plot(testz,histcurve,'k-');
    legstr{iii}='N';
    iii=iii+1; hold on;

    colrs1='brgcmy';
    colrs2='cmybrg';

    for qqq=1:length(endyears)
        t=find(targyears==endyears(qqq));
        
        if showESLR
            hl(iii)=plot(testz,effcurveESLR(t,:),[colrs1(qqq) '--']);
            legstr{iii}=['N+E(SL_{' num2str(endyears(qqq)) '})']; iii=iii+1;
        end
        
        if showeffcurve
            hl(iii)=plot(testz,effcurve(t,:), [colrs1(qqq) '-']);
            legstr{iii} = ['N_e(' num2str(endyears(qqq)) ')'];  iii=iii+1;
        end
        
        if show999
            hl(iii)=plot(testz,effcurve999(t,:),[colrs1(qqq) ':'] );
            legstr{iii}=['N+SL_{99.9}(' num2str(endyears(qqq)) ')']; iii=iii+1;
        end
        
    end

    for qqq=1:length(endyears)
        
        if showinteffcurve
            hl(iii)=plot(testz,integratecurve(effcurve,startyear,endyears(qqq)),[colrs2(qqq) '-']);
            legstr{iii}=['N_e(' num2str(startyear) ',' num2str(endyears(qqq)) ')'];
            iii=iii+1;
        end
        
    end

    hld=legend(hl,legstr,'Location','Northeast');
    set(hld,'fontsize',7);

    set(gca,'ysc','log');
    ylim([1e-4 10]);

    title(sitelab);
    xlabel('meters'); ylabel('expected events/year');

end

function parseFields(params)

flds=fieldnames(params);
for qqq=1:length(flds)
    if flds{qqq} 
        assignin('caller',flds{qqq},params.(flds{qqq}));
    end
end
