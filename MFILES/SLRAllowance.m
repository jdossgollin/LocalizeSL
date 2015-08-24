function [Ainst,z0,ESLR,poissonL,expectedN]=SLRAllowance(samps,N0,threshold,scale,shape,poissonL,beta,MHHW)

% [Ainst,z0,ESLR,poissonL,expectedN]=SLRAllowance(samps,N0,threshold,scale,shape,[poissonL],[beta],[MHHW])
%
% Calculate instantaneous sea-level rise allowances a la
% Hunter et al. (2012, 2013) assuming samples of local sea-level
% rise as specified and a GPD distribution of extremes.
%
%  INPUTS
%
%     samps: matrix of SLR samples (rows: samples, cols: time)
%     N0: original expected number (e.g., 0.01 for 1% flood)
%     threshold: GPD threshold
%     scale: GPD scale factor
%     shape: GPD shape factor
%     poissonL: mean number of exceedances per year
%               (default: top 1% of days, so 3.6525)
%               alternatively, provide a pair of [N ht]
%               and will calculate poissonL
%               (e.g., for 1% flood of 2000 mm, [0.01 2000])
%     beta: degree of confidence (1-beta is weight put on 99.9th percentile of samps)
%     MHHW: If specified, assume exceedances below threshold
%           follow a Gumbel distribution, with poissonL exceedances
%           at threshold and MHHW exceedances at z = 0 (e.g., 365.25
%           for 1 flood per day at z = 0).
%
%  OUTPUTS
%
%     Ainst: instantaneous sea-level rise allowance
%     z0: Height with initial expected number N0
%     ESLR: expected sea-level rise
%     poissonL: mean number of exceedances per year
%     expectedN: expectedN with allowances (for debugging)
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
%             LocalizeStoredProjections(selectedS[Ainst,z0,ESLR,poissonL,expectedN]=SLRAllowance(samps,N0,threshold,scale,shape,poissonL,beta,MHHW)ite,corefile,1);
%    samps=[zeros(size(sampslocrise{1,1},1),1) ...
%           sampslocrise{1,1}]/1000; % add base yr, convert to meters
%    samps=bsxfun(@min,samps,quantile(samps,.999));  
%                % truncate samples viewed as physically implausible
%    targyears = [2000 targ[Ainst,z0,ESLR,poissonL,expectedN]=SLRAllowance(samps,N0,threshold,scale,shape,poissonL,beta,MHHW)years]; % add base year
%
%
%    N0 = 0.01; % we will calculate allowances for the current 1% level
%    beta = 0.95; % 95% confidence in SLR PDF, 5% extra wt to worst case
%    MHHW = 365.25; % 365.25 floods/year for z = 0
%
%    [Ainst,z0,ESLR]=SLRAllowance(samps,N0, ...
%    threshold,scale,shape,lambda,beta,MHHW);
%
%    % now create a time series of integrated allowances
%    % from 2020:2100
%
%    t1 = 2020; t2 = 2100;
%    Aint = cumsum(interp1(targyears,Ainst(:,j),t1:t2))./([t1:t2]-t1+1);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 24 17:57:23 EDT 2015

defval('N0',.01);
%defval('poissonL',365.25*24*.01); % top 1% of hourlies
defval('poissonL',3.6525); % top 1% of days
defval('beta',1);
defval('MHHW',[]);

% if we're given a tuple instead of lambda, get ready to calibrate it
docalib = 0;
if length(poissonL)>1
    docalib = 1;
    Nref = poissonL(1);
    htref = poissonL(2);
    poissonL = 1;
end

% if we have a vector, assume it is a vector of samples from a single time point

if prod(size(samps))==length(samps)
    samps=samps(:);
end

% define N and calculate z0 that goes with N0

if length(MHHW)>0
    MHHW=[-threshold MHHW(1)];
end
if shape==0
    logN = @(z) GPDLogNExceedances(z,poissonL,shape,scale,MHHW);
    if docalib
       scaler = log(Nref)-logN(htref-threshold);
       poissonL = exp(scaler+log(poissonL));
       logN = @(z) GPDLogNExceedances(z,poissonL,shape,scale,MHHW);
    end
    z0 = -scale*log(N0)/log(poissonL) + threshold;
else
    logN = @(z) GPDLogNExceedances(z,poissonL,shape,scale,MHHW);
    if docalib
       scaler = log(Nref)-logN(htref-threshold);
       poissonL = exp(scaler+log(poissonL));
       logN = @(z) GPDLogNExceedances(z,poissonL,shape,scale,MHHW);
    end
    z0 = ((N0/poissonL)^(-shape)-1) * scale/shape + threshold;
end

% calculate expected sea-level rise
ESLR=mean(samps,1);
SLR999=quantile(samps,.999);

% find the allowance that minimizes the difference between the expected number
% of floods with a given allowance and A0

for ttt=1:size(samps,2)
    if shape>=0
        mx=max(samps(:));
    else
        mx = -.999*scale/shape;
    end

    %    mx2 = min(max(samps(:,ttt)),max(mx,mx-(z0-threshold-max(samps(:,ttt)))));
    mx2 = min(max(samps(:,ttt)),min(mx-(z0-threshold-max(samps(:,ttt)))));
    EN = @(A) checkExpectedN(logN,z0-threshold + A - samps(:,ttt));
    EN999 = @(A) checkExpectedN(logN,z0-threshold+A-SLR999(ttt));
    Ainst(ttt) = fminbnd(@(A) abs( log(beta*EN(A)+(1-beta)*EN999(A))-log(N0)),0,mx2,optimset('maxfunevals',20000));
    expectedN(ttt) = (beta*EN(Ainst(ttt))+(1-beta)*EN999(Ainst(ttt)));
    if (abs(expectedN(ttt)-N0)/N0)>.01
            Ainst(ttt)=NaN;
    end
end

end

function EN=checkExpectedN(logN,y)
    EN = mean(exp(logN(max(0,y))));
    if imag(EN)>(1e-6*real(EN))
        EN=1e10;
    else
        EN=real(EN);
    end
end

