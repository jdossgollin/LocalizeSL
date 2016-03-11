function [A,z0]=SLRAllowanceFromCurve(testz,histcurve,effcurve,N0)

% [A,z0]=SLRAllowanceFromCurve(testz,histcurve,effcurve,N0)
%
% Calculate sea-level rise allowances a la Hunter et al. (2012, 2013)
% and Buchanan et al. (2016) from an expected flood frequency curve
% returned by SLRFloodNexpVsLevelCurves.
%
%  INPUTS
%
%     testz: vector of flood heights used by SLRFloodNexpVsLevelCurves
%     histcurve: vector of historic frequencies at heights testz
%     effcurve: matrix (rows = years, cols = heights) of expected
%               numbers of floods
%     N0: original expected number (e.g., 0.01 for 1% flood)
%
%  OUTPUTS
%
%     A: sea-level rise allowance
%     z0: Height with initial expected number N0
%
% EXAMPLE:
%
%       sitelab='Boston';
%       selectedSite = 235; %PSMSL ID for Boston
%       threshold = 0.6214; % GPD threshold
%       scale = 0.1109; % GPD scale
%       shape = 0.0739; % GPD shape
%       lambda = 2.5699; % Poisson Lambda
%  
%       [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = ...
%                LocalizeStoredProjections(selectedSite,corefile,1);
%       samps=[zeros(size(sampslocrise{1,1},1),1) ...
%              sampslocrise{1,1}]/1000; % add base yr, convert to meters
%       samps=bsxfun(@min,samps,quantile(samps,.999));
%                   % truncate samples viewed as physically implausible
%       targyears = [2000 targyears]; % add base year
%  
%       pm.doplot=0;
%       [effcurve,testz,histcurve,histcurvesamps,effcurveESLR,effcurve999,integratecurve]= ...
%           SLRFloodNexpVsLevelCurves(samps,targyears,threshold, ...
%           scale,shape,lambda,sitelab,pm);     
%  
%       N0s = [0.1 0.01 .002]; % we will calculate allowances for the current 10%, 1%, and .2% levels
%       beta = 0.95; % 95% confidence in SLR PDF, 5% extra wt to worst case
%   
%       weffcurve=beta*effcurve+(1-beta)*effcurve999;
%       for www=1:length(N0s)
%           [Ainst(:,www),z0]=SLRAllowanceFromCurve(testz,histcurve,weffcurve,N0s(www));
%  
%           % now create a time series of design-life allowances
%           % from 2020:2100
%  
%           t1 = 2020; t2 = 2100;
%           ADLendyears=targyears(find((targyears>t1).*(targyears<=t2)));
%           for ttt=1:length(ADLendyears)
%               ADL(ttt,www)=SLRAllowanceFromCurve(testz,histcurve, ...
%                  integratecurve(weffcurve,t1,ADLendyears(ttt)),N0s(www));
%           end
%       end
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Mar 11 15:57:27 EST 2016

for qqq=1:length(N0)
    [m,mi]=min(abs(histcurve-N0(qqq)));
    z0(qqq)=testz(mi);
    [m,mi]=min(abs(effcurve-N0(qqq)),[],2);
    A(:,qqq)=testz(mi)-z0(qqq);
end