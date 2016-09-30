function [projLOCcomp,projLOCcomphi,projLOCcomplo,colsCOMP,colsCOMPlab]=ConditionalDistributionsLSLComponents(p,condsubscen,substitutep)

% [projLOCcomp,projLOCcomphi,projLOCcomplo,colsCOMP,colsCOMPlab]=ConditionalDistributionsLSLComponents(p,condsubscen,substitutep)
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
% projLOCcomp: median LSL projection for each scenario
% projLOChicomp: high LSL projection for each scenario
% projLOClocomp: low LSL projection for each scenario
% colsCOMP: columns of core files used for compribution breakdown
% colsCOMPlab: labels for compribution breakdown
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 17:41:10 EDT 2016

defval('Nslice',20);

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

colsCOMP={p.colGIC,p.colTE,p.colGIS,p.colAIS,p.colTE}; colsCOMPlab={'GIC','TE','GIS','AIS','Oc'};


for sss=1:length(slicesub)
    slicedp=slicep(p,slicesub{sss});
    clear wwprojLOC wwprojLOChi wwprojLOClo
    parfor www=1:length(slicedp.targregions)
        selectedSite=targregions(slicesub{sss}(www));
        [~,wsampscomp] = LocalizeStoredProjections(selectedSite,slicedp,[1 3 4],substitutep);
        wlocalsampscomp=[wsampscomp{1}; wsampscomp{2} ; wsampscomp{3}];
        wprojLOC = zeros(length(condsubscen),length(targyears),length(colsCOMP));
        wprojLOChi = wprojLOC;
        wprojLOClo = wprojLOC;
        for ttt=1:length(colsCOMP)
            wlocalsamps=squeeze(sum(wlocalsampscomp(:,colsCOMP{ttt},:),2));
            for qqq=1:length(condsubscen)
                wprojLOC(qqq,:,ttt)=quantile(wlocalsamps(condsubscen{qqq},:),.5);
                wprojLOChi(qqq,:,ttt)=quantile(wlocalsamps(condsubscen{qqq},:),.833);
                wprojLOClo(qqq,:,ttt)=quantile(wlocalsamps(condsubscen{qqq},:),.167);
            end
        end      
        wwprojLOC(:,:,:,www)=wprojLOC;
        wwprojLOChi(:,:,:,www)=wprojLOChi;
        wwprojLOClo(:,:,:,www)=wprojLOClo;

    end
    
    projLOCcomp(:,:,:,slicesub{sss})=wwprojLOC;
    projLOCcomphi(:,:,:,slicesub{sss})=wwprojLOChi;
    projLOCcomplo(:,:,:,slicesub{sss})=wwprojLOClo;
end
