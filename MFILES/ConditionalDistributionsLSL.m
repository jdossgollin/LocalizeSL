function [projLOC,projLOChi,projLOClo,projLOCrate,projLOCratehi,projLOCratelo]=ConditionalDistributionsLSL(p,condsubscen,substitutep)

% [projLOC,projLOChi,projLOClo,projLOCrate,projLOCratehi,projLOCratelo]=ConditionalDistributionsLSL(p,condsubscen,substitutep)
%
% Generate local sea-level scenarios by conditionalizing probabilistic projections.
%
% INPUT
% -----
% p: core sea-level structure 
% condsubscen: cell array of sets of sample indices for each scenario
% substitutep: subtitutions to make in p
%
% OUTPUT
% ------
% projLOC: median LSL projection for each scenario
% projLOChi: high LSL projection for each scenario
% projLOClo: low LSL projection for each scenario
% projLOCrate: median rate projections for each scenario
% projLOCratehi: high rate projections for each scenario
% projLOCratelo: low rate projections for each scenario
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 10:00:57 EDT 2016

defval('difftimestep',20);
defval('Nslice',20);

Nslice=20;
fullindex=1:length(p.targregions);
cnt=1;
slicesub{1}=1:min(Nslice,length(fullindex));
lastmax=max(slicesub{end});
while lastmax<max(fullindex)
    slicesub{end+1}=((lastmax+1):min(Nslice+lastmax,length(fullindex)));
    lastmax=max(slicesub{end});
end

targyears=p.targyears;
targregions=p.targregions;


for sss=1:length(slicesub)
    slicedp=slicep(p,slicesub{sss});
    clear wwprojLOC wwprojLOChi wwprojLOClo wwprojLOCrate wwprojLOCratehi wwprojLOCratelo;
    dtstep=difftimestep;
    parfor www=1:length(slicedp.targregions)
        selectedSite=targregions(slicesub{sss}(www));
        [wsamps] = LocalizeStoredProjections(selectedSite,slicedp,[1 3 4],substitutep);
        wlocalsamps=[wsamps{1}; wsamps{2} ; wsamps{3}];
        [wsampsdrise,wtargyeardrates]=SampleSLRates(wsamps,targyears,dtstep);
        wlocalsamprates=[wsampsdrise{1} ; wsampsdrise{2} ;wsampsdrise{3}];
        
        wprojLOC = zeros(length(condsubscen),length(targyears));
        wprojLOChi = wprojLOC;
        wprojLOClo = wprojLOC;
        
        wprojLOCrate = zeros(length(condsubscen),length(targyears)-1);
        wprojLOCratehi = wprojLOCrate;
        wprojLOCratelo = wprojLOCrate;

        for qqq=1:length(condsubscen)
            wprojLOC(qqq,:)=quantile(wlocalsamps(condsubscen{qqq},:),.5);
            wprojLOChi(qqq,:)=quantile(wlocalsamps(condsubscen{qqq},:),.833);
            wprojLOClo(qqq,:)=quantile(wlocalsamps(condsubscen{qqq},:),.167);
            wprojLOCrate(qqq,:)=quantile(wlocalsamprates(condsubscen{qqq},:),.5);
            wprojLOCratehi(qqq,:)=quantile(wlocalsamprates(condsubscen{qqq},:),.833);
            wprojLOCratelo(qqq,:)=quantile(wlocalsamprates(condsubscen{qqq},:),.167);
        end
        wwprojLOC(:,:,www)=wprojLOC;
        wwprojLOChi(:,:,www)=wprojLOChi;
        wwprojLOClo(:,:,www)=wprojLOClo;
        wwprojLOCrate(:,:,www)=wprojLOCrate;
        wwprojLOCratehi(:,:,www)=wprojLOCratehi;
        wwprojLOCratelo(:,:,www)=wprojLOCratelo;

    end
    
    projLOC(:,:,slicesub{sss})=wwprojLOC;
    projLOChi(:,:,slicesub{sss})=wwprojLOChi;
    projLOClo(:,:,slicesub{sss})=wwprojLOClo;
    projLOCrate(:,:,slicesub{sss})=wwprojLOCrate;
    projLOCratehi(:,:,slicesub{sss})=wwprojLOCratehi;
    projLOCratelo(:,:,slicesub{sss})=wwprojLOCratelo;
end
