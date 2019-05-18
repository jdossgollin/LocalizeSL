% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Mar 16 11:59:39 EDT 2016

rootdir='~/Dropbox/Code/LocalizeSL';
corefile=fullfile(rootdir,'IFILES/SLRProjections140523core.mat');
addpath(fullfile(rootdir,'MFILES'));

extremesIFILES = fullfile(rootdir,'IFILES/extremes');
datDir = fullfile(extremesIFILES, 'declustered');


% import table of calibration parameters
parmdat=importdata(fullfile(extremesIFILES,'GPDfits_withunc.tsv'),'\t',1);
extracols=size(parmdat.textdata,2)-size(parmdat.data,2);
lambdas=parmdat.data(:,find(strcmpi('lambda',parmdat.textdata(1,:)))-extracols);
thresholds=parmdat.data(:,find(strcmpi('u',parmdat.textdata(1,:)))-extracols);
scales=parmdat.data(:,find(strcmpi('scale',parmdat.textdata(1,:)))-extracols);;
shapes=parmdat.data(:,find(strcmpi('shape',parmdat.textdata(1,:)))-extracols);
AEP10pts=parmdat.data(:,find(strcmpi('AEP0.1',parmdat.textdata(1,:)))-extracols);
Vscale=parmdat.data(:,find(strcmpi('Vscale',parmdat.textdata(1,:)))-extracols);
Vshape=parmdat.data(:,find(strcmpi('Vshape',parmdat.textdata(1,:)))-extracols);
Vscaleshape=parmdat.data(:,find(strcmpi('Vscaleshape',parmdat.textdata(1,:)))-extracols);
psmslids=parmdat.data(:,1);
tgids=parmdat.data(:,2);
NOAAnames=parmdat.textdata(2:end,1);

% find instantaneous allowances with specified GPD fits
clear effcurve histcurve effcurve999 integratecurve;
for qqq=1:length(psmslids)
    selectedSite = psmslids(qqq); % use PSMSL ID here to select site
    wfile=fullfile(datDir,['maxtofit.dclist.' num2str(tgids(qqq)) '_xdat.dc.tsv']);
    historicaldata=importdata(wfile);

try
    [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,1);
    %longname=sitenames{1};
    %shortname=longname(setdiff(1:length(longname),strfind(longname,' ')));
    
    longname=NOAAnames{qqq};
    shortname=longname(setdiff(1:length(longname),strfind(longname,' ')));
    shortname(strfind(shortname,','))='-';
    
    samps=[zeros(size(sampslocrise{1,1},1),1) sampslocrise{1,1}]/1000; % add base year and convert to meters
    samps=bsxfun(@min,samps,quantile(samps,.999)); % truncate samples viewed as physically implausible
    targyears = [2000 targyears];
    
    acov = [Vscale(qqq) Vscaleshape(qqq) ; Vscaleshape(qqq) Vshape(qqq)];
    parmsamps=lhsnorm([scales(qqq) shapes(qqq)],acov,1000);
    parmsamps(:,1)=max(eps,parmsamps(:,1));
       
    clf;
    clear pm;
    pm.showuncertainty=1; pm.historicaldata=historicaldata;
         [effcurve{qqq},testz,histcurve{qqq},histcurvesamps,effcurveESLR,effcurve999{qqq},integratecurve]= ...
         SLRFloodNexpVsLevelCurves(samps,targyears,thresholds(qqq), ...
         parmsamps(:,1),parmsamps(:,2),lambdas(qqq),longname,pm);     
      
    pdfwrite([shortname '_returncurves']);

    clf; 
    clear pm; pm.betas=[0 .5 .9 1]; [Ainst,ALDC,ADLfromstart,ADLfp,ADLLDCfromstart,ADLendyears,z0,hp]=SLRAllowancePlot(samps,targyears,effcurve{qqq},testz,histcurve{qqq},effcurve999{qqq},integratecurve,longname,pm);
    pdfwrite([shortname '_Allowances']);

[A,ADL,z0]=SLRAllowanceWriteTable([shortname '_Allowances'],targyears,effcurve{qqq},testz,histcurve{qqq},effcurve999{qqq},integratecurve,[.01 .1 .002],[1 .9 .5 .0],[2100],longname)
end
end    

save effcurves targyears effcurve testz histcurve effcurve999 integratecurve