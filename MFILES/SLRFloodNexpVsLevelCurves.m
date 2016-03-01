function [effcurve,testz]=SLRFloodNexpVsLevelCurves(samps,targyears,threshold,scale,shape,lambda,sitelab,shortname,params)

% [effcurve,testz]=SLRFloodNexpVsLevelCurves(samps,targyears,threshold,scale,shape,lambda,sitelab,shortname,params)
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
%     shortname: short name of site (for filename; no spaces)
%     params: an optional structure with fields setting several parameters
%         - startyear: start year for AADLL curves (default: 2020) 
%         - endyears: end years for AADLLs; instantaneous effective flood risk
%                     for these years are also plotted
%                     (default: [2050 2100])
%         - testz: vector of heights used for calculating the curve
%                  (default: [0:.01:10])
%         - doplot: generate plot (default: 1)
%         - effcurve: if passed, will skip calculation and use as specified
%
% OUTPUTS:
%
%     effcurve: expected number of floods of different heights (rows: years as
%               in targyears; cols: heights as specified in testz)
%     testz: heights for effcurve
%
% EXAMPLE:
%
%    shortname='NYC';
%    selectedSite = 12; %PSMSL ID for New York City
%    threshold = 0.5148;
%    scale = 0.1285; % GPD scale
%    shape = 0.1879; % GPD shape
%    lambda = 2.8116; % Poisson Lambda
%    
%    [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = ...
%             LocalizeStoredProjections(selectedSite,corefile,1);
%    samps=[zeros(size(sampslocrise{1,1},1),1) ...
%           sampslocrise{1,1}]/1000; % add base year and convert to meters
%    samps=bsxfun(@min,samps,quantile(samps,.999));  
%                      % truncate samples viewed as physically implausible
%    targyears = [2000 targyears]; % add base year
%      
%    clf;
%    [effcurve,testz]=SLRFloodNexpVsLevelCurves(samps,targyears,threshold, ...
%                     scale,shape,lambda,sitenames{1},shortname);     
%    pdfwrite([shortname '_returncurves']);
%
% Note that, to calculate am AADLL curve for a design life starting at t1 and
% ending at t2 from the effcurve matrix, simply run:
%
%    mean(interp1(targyears,effcurve,t1:t2))
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Mar 01 15:16:33 EST 2016

defval('siteshortname','SL');
defval('sitelab',[]);
defval('testz',0:.01:10);
defval('startyear',2020);
defval('endyears',[2050 2100]);
defval('doplot',1);
defval('effcurve',[]);
defval('showuncertainty',0);
defval('historicaldata',[]);

if exist('params')
    parseFields(params)
end

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

if showuncertainty
    if size(shape,1)>1
        for iii=1:size(shape,1)
            logNs = @(z) GPDLogNExceedances(z-threshold,lambda,shape(iii),scale(iii),-threshold);
            curve0s(iii,:) = real(exp(logNs(testz)));
        end
    else
        curve0s = effcurve(1,:);
    end
end

integrateeffcurve=@(t1,t2) (1/(t2-t1)) * sum(interp1(targyears,effcurve,t1:t2));

ESLR = mean(samps);
SLR999=quantile(samps,.999);

%%%

if doplot

    hp1(1)=subplot(2,1,1);

    iii=1;
    if showuncertainty
        if size(shape,1)>1
            greyc=[.7 .7 .7];
            plot(testz,quantile(curve0s,.5),'-','Color',greyc);
            hold on;
            plot(testz,quantile(curve0s,.167),'--','Color',greyc);
            plot(testz,quantile(curve0s,.833),'--','Color',greyc);
        end
    end
    
    if length(historicaldata)>0
        ct=sum(bsxfun(@gt,historicaldata,testz))/(length(historicaldata)/365.25);
        subct=find(diff(ct)<0);
        subct=intersect(subct,find(testz>threshold));
        plot(testz(subct),ct(subct),'ys');
    end
    
    hl(iii)=plot(testz,exp(logN(testz)),'k-');
    legstr{iii}='N';
    iii=iii+1; hold on;

    colrs1='brg';
    colrs2='cmy';

    for qqq=1:length(endyears)
        t=find(targyears==endyears(qqq));
        
        hl(iii)=plot(testz,exp(logN(testz-ESLR(t))),[colrs1(qqq) '--']);
        legstr{iii}=['N+E(SL_{' num2str(endyears(qqq)) '})']; iii=iii+1;

        hl(iii)=plot(testz,effcurve(t,:), [colrs1(qqq) '-']);
        legstr{iii} = ['N_e(' num2str(endyears(qqq)) ')'];  iii=iii+1;
        
        hl(iii)=plot(testz,exp(logN(testz-SLR999(t))),[colrs1(qqq) ':'] );
        legstr{iii}=['N+SL_{99.9}(' num2str(endyears(qqq)) ')']; iii=iii+1;
    end

    for qqq=1:length(endyears)
        
        hl(iii)=plot(testz,integrateeffcurve(startyear,endyears(qqq)),[colrs2(qqq) '-']);
        legstr{iii}=['N_e(' num2str(startyear) ',' num2str(endyears(qqq)) ')'];
        iii=iii+1;
        
    end

    hld=legend(hl,legstr,'Location','Northeast');
    set(hld,'fontsize',7);

    set(gca,'ysc','log');
    ylim([1e-4 10]);

    title(sitelab);
    xlabel('meters'); ylabel('effective annual expected');

end

function parseFields(params)

flds=fieldnames(params);
for qqq=1:length(flds)
    if flds{qqq} 
        assignin('caller',flds{qqq},params.(flds{qqq}));
    end
end
