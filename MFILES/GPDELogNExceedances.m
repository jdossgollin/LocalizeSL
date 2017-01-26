function y = GPDELogNExceedances(z,lambda,shape,scale,MHHW,maxinterpN)

% y = GPDELogNExceedances(z,lambda,shape,scale,[MHHW],[maxinterpN])
%
% Calculate log of the expected number of exceedances of z from a
% Poisson-Generalized Pareto Distribution with Poisson mean
% lambda and the specified shape and scale factors. Assumes
% the exceedance threshold has already been removed from z.
%
% shape and scale are samples from the distribution.
%
% This function ensures that values returned stay within the
% support of the GPD.
%
% For values of z below zero, the function will treat as though
% z = 0 unless MHHW is specified. If MHHW is specified
% value, exceedances below zero will be assumed to fall on a
% Gumbel distribution between lambda exceedances at zero and a
% specified value at z = MHHW(1) < 0. If MHHW(2) exists, it is
% the value of exceedances at z = MHHW(1); otherwise, will default
% to 365.25.
%
% Developed for Buchanan et al. (2016).
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Mar 01 15:47:11 EST 2016

defval('MHHW',[]);
defval('maxinterpN',1000);

if size(shape)==1
     y = GPDLogNExceedances(z,lambda,shape,scale,MHHW);
else
    
    minz=.99999*abs(scale./shape);
    logN = @(z) log(lambda.*bsxfun(@power,1+bsxfun(@times,(shape+eps)./scale,max(0,z)),-1./(shape+eps)));

    % interpolate so we can deal with z of different dimensions
    testz=unique([0  z(:)']);
    if length(testz)>maxinterpN
        testz=linspace(testz(1)-eps,testz(end)+eps,maxinterpN);
    elseif length(testz)==1
        testz=[testz-eps testz+eps];
    end
    testz0=testz;
    testz=bsxfun(@times,(shape>=0),testz) + bsxfun(@times,(shape<0),bsxfun(@min,testz,minz));
    logNref=logN(testz);
    logENref = log(mean(exp(logNref)));
    logN = @(z) interp1(testz0,logENref,z,'linear');

    z1=max(0,z);
    y=logN(z1);

    % for those points below threshold to MHHW, put on a Gumbel
    if length(MHHW)>=1
        z=max(z,MHHW(1));
        if length(MHHW)>=2
            MHHWfreq=MHHW(2);
        else
            MHHWfreq=365.25/2;
        end  
        y = (z<0).*(log(lambda)+(log(MHHWfreq)-log(lambda)).*z/MHHW(1))+(z>=0).*y;
    end
end
