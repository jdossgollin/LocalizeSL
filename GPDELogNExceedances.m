function y = GPDLogNExceedances(z,lambda,shape,scale,MHHW)

% y = GPDLogNExceedances(z,lambda,shape,scale,[MHHW])
%
% Calculate log of the number of exceedances of z from a
% Poisson-Generalized Pareto Distribution with Poisson mean
% lambda and the specified shape and scale factor. Assumes
% the exceedance threshold has already been removed from z.
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
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 29 23:06:12 EST 2016

defval('MHHW',[]);

% $$$ if shape<0
% $$$     z=min(z,.99999*-scale/shape);
% $$$     logN = @(z) log(lambda.*(1+(shape).*z/scale).^(-1./(shape)));
% $$$ elseif shape==0
% $$$         logN = @(z) log(lambda)-z./scale;
% $$$ else
% $$$     logN = @(z) log(lambda.*(1+(shape).*z./scale).^(-1./(shape)));
% $$$ end

minz=.99999*abs(scale./shape);
if size(shape,1)>1
    z=bsxfun(@times,(shape>=0),z) + bsxfun(@times,(shape<0),bsxfun(@min,z,minz));
    logN = @(z) log(lambda.*bsxfun(@power,1+bsxfun(@times,(shape+eps)./scale,max(0,z)),-1./(shape+eps)));
else
    z=(shape>=0).*z + (shape<0).*min(z,minz);
    logN = @(z) log(lambda.*(1+(shape+eps).*max(0,z)./scale).^(-1./(shape+eps)));
end
y=logN(max(0,z));
keyboard
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
