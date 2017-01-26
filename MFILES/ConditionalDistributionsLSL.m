function projections=ConditionalDistributionsLSL(p,condsubscen,substitutep)

% [projections]=ConditionalDistributionsLSL(p,condsubscen,substitutep)
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
% projections: projection structure
%
% Fields of projection structure
% ------------------------------
% projLOC: median LSL projection for each scenario
% projLOChi: high LSL projection for each scenario
% projLOClo: low LSL projection for each scenario
% projLOCrate: median rate projections for each scenario
% projLOCratehi: high rate projections for each scenario
% projLOCratelo: low rate projections for each scenario
% projLOC0: median LSL projection for each scenario, excluding background
% projLOC0hi: high LSL projection for each scenario, excluding background
% projLOC0lo: low LSL projection for each scenario, excluding background
% projLOC0rate: median rate projections for each scenario, excluding background
% projLOC0ratehi: high rate projections for each scenario, excluding background
% projLOC0ratelo: low rate projections for each scenario, excluding background
% projLOCcomp: median LSL component projection for each scenario
% projLOCcomphi: high LSL component projection for each scenario
% projLOCcomplo: low LSL component projection for each scenario
% colsCOMP: component columns
% colsCOMPlab: component labels
%
% Developed for Sweet et al. (2017).
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Nov 02 00:19:24 EDT 2016

defval('difftimestep',20);
defval('Nslice',20);
defval('substitutep',[]);

Nbkgdsamps=17;
docomponents=1;

if length(substitutep)==0
    clear substitutep; substitutep.filler=0;
end


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

colsCOMP={p.colGIC,p.colGIS,p.colAIS,p.colTE}; colsCOMPlab={'GIC','GIS','AIS','Oc'};

% turn off background rates
substitutep.rateprojs=p.rateprojs*0;
substitutep.rateprojssd=p.rateprojssd*0;

for sss=1:length(slicesub)
    slicedp=slicep(p,slicesub{sss});
    clear wwprojLOC wwprojLOChi wwprojLOClo wwprojLOCrate wwprojLOCratehi wwprojLOCratelo;
    clear wwprojLOC0 wwprojLOC0hi wwprojLOC0lo wwprojLOC0rate wwprojLOC0ratehi wwprojLOC0ratelo;
    clear wwprojLOCc wwprojLOCchi wwprojLOCclo;
    
    dtstep=difftimestep;
    parfor www=1:length(slicedp.targregions)
        selectedSite=targregions(slicesub{sss}(www));
        if docomponents
            [wsamps,wsampscomp] = LocalizeStoredProjections(selectedSite,slicedp,[1 3 4],substitutep);
        else
            [wsamps] = LocalizeStoredProjections(selectedSite,slicedp,[1 3 4],substitutep);
        end
        
        wlocalsamps=[wsamps{1}; wsamps{2} ; wsamps{3}];
        [wsampsdrise,wtargyeardrates]=SampleSLRates(wsamps,targyears,dtstep);
        wlocalsamprates=[wsampsdrise{1} ; wsampsdrise{2} ;wsampsdrise{3}];
        
        wprojLOC0 = zeros(length(condsubscen),length(targyears));
        wprojLOC0hi = wprojLOC0;
        wprojLOC0lo = wprojLOC0;
        
        wprojLOC = zeros(length(condsubscen),length(targyears));
        wprojLOChi = wprojLOC;
        wprojLOClo = wprojLOC;
                
        wprojLOCrate = zeros(length(condsubscen),length(targyears)-1);
        wprojLOCratehi = wprojLOCrate;
        wprojLOCratelo = wprojLOCrate;
        
        wprojLOC0rate = zeros(length(condsubscen),length(targyears)-1);
        wprojLOC0ratehi = wprojLOC0rate;
        wprojLOC0ratelo = wprojLOC0rate;

        u=norminv(linspace(0,1,Nbkgdsamps+2)); u=u(2:end-1);
        wbkgdratesamps=reshape(slicedp.rateprojs(www)+u*slicedp.rateprojssd(www),1,1,[]);
        wbkgdlevels=bsxfun(@times,slicedp.targyears-2000,wbkgdratesamps);
        
        if docomponents
            wlocalsampscomp=[wsampscomp{1}; wsampscomp{2} ; wsampscomp{3}];
            wprojLOCc = zeros(length(condsubscen),length(targyears),length(colsCOMP));
            wprojLOCchi = wprojLOCc;
            wprojLOCclo = wprojLOCc;
        end
        
        for qqq=1:length(condsubscen)
            wlsamps0 = wlocalsamps(condsubscen{qqq},:);
            wlsamps=bsxfun(@plus,repmat(wlsamps0,1,1,Nbkgdsamps),wbkgdlevels);
            wlsamps=reshape(permute(wlsamps,[3 1 2]),Nbkgdsamps*size(wlsamps0,1),[]);
            
            wlsamprates0 = wlocalsamprates(condsubscen{qqq},:);
            wlsamprates = bsxfun(@plus,repmat(wlsamprates0,1,1,Nbkgdsamps),repmat(wbkgdratesamps,size(wlsamprates0,1),size(wlsamprates0,2)));
            wlsamprates=reshape(permute(wlsamprates,[3 1 2]),Nbkgdsamps*size(wlsamprates0,1),[]);
            
            wprojLOC(qqq,:)=quantile(wlsamps,.5);
            wprojLOChi(qqq,:)=quantile(wlsamps,.833);
            wprojLOClo(qqq,:)=quantile(wlsamps,.167);
            wprojLOCrate(qqq,:)=quantile(wlsamprates,.5);
            wprojLOCratehi(qqq,:)=quantile(wlsamprates,.833);
            wprojLOCratelo(qqq,:)=quantile(wlsamprates,.167);
            
            wprojLOC0(qqq,:)=quantile(wlsamps0,.5);
            wprojLOC0hi(qqq,:)=quantile(wlsamps0,.833);
            wprojLOC0lo(qqq,:)=quantile(wlsamps0,.167);
            wprojLOC0rate(qqq,:)=quantile(wlsamprates0,.5);
            wprojLOC0ratehi(qqq,:)=quantile(wlsamprates0,.833);
            wprojLOC0ratelo(qqq,:)=quantile(wlsamprates0,.167);

        end
        wwprojLOC(:,:,www)=wprojLOC;
        wwprojLOChi(:,:,www)=wprojLOChi;
        wwprojLOClo(:,:,www)=wprojLOClo;
        wwprojLOCrate(:,:,www)=wprojLOCrate;
        wwprojLOCratehi(:,:,www)=wprojLOCratehi;
        wwprojLOCratelo(:,:,www)=wprojLOCratelo;
        
        wwprojLOC0(:,:,www)=wprojLOC0;
        wwprojLOC0hi(:,:,www)=wprojLOC0hi;
        wwprojLOC0lo(:,:,www)=wprojLOC0lo;
        wwprojLOC0rate(:,:,www)=wprojLOC0rate;
        wwprojLOC0ratehi(:,:,www)=wprojLOC0ratehi;
        wwprojLOC0ratelo(:,:,www)=wprojLOC0ratelo;
        
        if docomponents
            for ttt=1:length(colsCOMP)
                wlocalsampsc=squeeze(sum(wlocalsampscomp(:,colsCOMP{ttt},:),2));
                for qqq=1:length(condsubscen)
                    wprojLOCc(qqq,:,ttt)=quantile(wlocalsampsc(condsubscen{qqq},:),.5);
                    wprojLOCchi(qqq,:,ttt)=quantile(wlocalsampsc(condsubscen{qqq},:),.833);
                    wprojLOCclo(qqq,:,ttt)=quantile(wlocalsampsc(condsubscen{qqq},:),.167);
                end
            end      

            wwprojLOCc(:,:,:,www)=wprojLOCc;
            wwprojLOCchi(:,:,:,www)=wprojLOCchi;
            wwprojLOCclo(:,:,:,www)=wprojLOCclo;
        end


    end

    projections.projLOC0(:,:,slicesub{sss})=wwprojLOC0;
    projections.projLOC0hi(:,:,slicesub{sss})=wwprojLOC0hi;
    projections.projLOC0lo(:,:,slicesub{sss})=wwprojLOC0lo;
    projections.projLOC0rate(:,:,slicesub{sss})=wwprojLOC0rate;
    projections.projLOC0ratehi(:,:,slicesub{sss})=wwprojLOC0ratehi;
    projections.projLOC0ratelo(:,:,slicesub{sss})=wwprojLOC0ratelo;
    
    projections.projLOC(:,:,slicesub{sss})=wwprojLOC;
    projections.projLOChi(:,:,slicesub{sss})=wwprojLOChi;
    projections.projLOClo(:,:,slicesub{sss})=wwprojLOClo;
    projections.projLOCrate(:,:,slicesub{sss})=wwprojLOCrate;
    projections.projLOCratehi(:,:,slicesub{sss})=wwprojLOCratehi;
    projections.projLOCratelo(:,:,slicesub{sss})=wwprojLOCratelo;
    
    if docomponents
        projections.projLOCcomp(:,:,:,slicesub{sss})=wwprojLOCc;
        projections.projLOCcomphi(:,:,:,slicesub{sss})=wwprojLOCchi;
        projections.projLOCcomplo(:,:,:,slicesub{sss})=wwprojLOCclo;

    end
    
end
projections.colsCOMP=colsCOMP; projections.colsCOMPlab=colsCOMPlab;
