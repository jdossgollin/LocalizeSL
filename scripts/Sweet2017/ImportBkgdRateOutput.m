function y=ImportBkgdRateOutput(datfile)
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Jan 11 14:40:38 EST 2017

dat=importdata(datfile);
y.sitelab=dat.textdata{:,1};
y.psmslid=dat.data(:,1);
y.lat=dat.data(:,2);
y.lon=dat.data(:,3);
y.bkgdrate=dat.data(:,4);
y.bkgdrate_std=dat.data(:,5);
