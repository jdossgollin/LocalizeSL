function [proj,condsubscen,projhi,projlo,targyearrates,projrate,projratehi,projratelo,projCONT,projCONThi,projCONTlo,colsCONT,colsCONTlab]=ConditionalDistributionsGSL(p,condtargyrs,condtargs,condtargwins,substitutep)

% [proj,condsubscen,projhi,projlo,projrate,projratehi,projratelo,projCONT,projCONThi,projCONTlo,colsCONT,colsCONTlab]=GSLConditionalDistributions(p,condtargyrs,condtargs,condtargwins,substitutep)
%
% Generate global sea-level scenarios by conditionalizing probabilistic projections.
%
% INPUT
% -----
% p: core sea-level structure 
% condtargyrs: years on which to condition GSL
% condtargs: target heights (in mm) upon which to condition;
%            rows correpsond to years and columns to targets
% condtargwins: tolerance (in mm) for deviation from condtargs;
%               same dimensions as condtargs
% substitutep: subtitutions to make in p
%
% OUTPUT
% ------
% proj: median GSL projection for each scenario
% condsubscen: cell array of sets of sample indices for each scenario
% projhi: high GSL projection for each scenario
% projlo: low GSL projection for each scenario
% targyearrates: years for rates
% projrate: median rate projections for each scenario
% projratehi: high rate projections for each scenario
% projratelo: low rate projections for each scenario
% projCONT: median contribution projections for each scenario
% projCONThi: high contribution projections for each scenario
% projCONTlo: low contribution projections for each scenario
% colsCONT: columns of core files used for contribution breakdown
% colsCONTlab: labels for contribution breakdown
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 11:10:19 EDT 2016

defval('condtargyrs',[2100 2050 2030]);
defval('condtargs',[30 50 100 150 200 250 ;
           15 NaN NaN NaN NaN NaN ;
           9 NaN NaN NaN NaN NaN] * 10);
defval('condtargwins',[20 20 20 20 50 50 ;
              10 10 10 10 10 10 ;
              5 5 5 5 5 5]);
defval('difftimestep',20);
defval('substitutep',[]);
if length(substitutep)==0
    substitutep.filler=0;
end


%%%%

[sampsGSLrise,sampsGSLcomponents,siteidsGSL,sitenamesGSL,targyearsGSL,scensGSL,colsGSL] = LocalizeStoredProjections(0,p,[1 3 4],substitutep);
[sampsdGSLrise,targyearGSLrates]=SampleSLRates(sampsGSLrise,targyearsGSL,difftimestep);

targyears=targyearsGSL; targyearrates=targyearGSLrates;
pooledGSL=[sampsGSLrise{1} ; sampsGSLrise{2} ;sampsGSLrise{3}];
pooledGSLcont=[sampsGSLcomponents{1} ; sampsGSLcomponents{2} ;sampsGSLcomponents{3}];
pooledGSLrate=[sampsdGSLrise{1} ; sampsdGSLrise{2} ;sampsdGSLrise{3}];


%%%%%

clear proj projhi projlo projrate projCONT projCONThi projCONTlo  projratehi projratelo;
colsCONT={colsGSL.colGIC,colsGSL.colTE,colsGSL.colGIS,colsGSL.colAIS,[1:23]}; colsCONTlab={'GIC','TE','GIS','AIS','all'};


clear condsubscen;
for qqq=1:size(condtargs,2)
    sub=1:size(pooledGSL,1);
    for rrr=1:length(condtargyrs)
       if ~isnan(condtargs(rrr,qqq)) sub=intersect(sub,find(abs(pooledGSL(:,find(targyears==condtargyrs(rrr)))-condtargs(rrr,qqq))<condtargwins(rrr,qqq)));
       end
       
    end
    
 
    condsubscen{qqq}=sub;
     proj(qqq,:)=quantile(pooledGSL(sub,:),.5);
     projhi(qqq,:)=quantile(pooledGSL(sub,:),.833);
     projlo(qqq,:)=quantile(pooledGSL(sub,:),.167);
     projrate(qqq,:)=quantile(pooledGSLrate(sub,:),.5);
     projratehi(qqq,:)=quantile(pooledGSLrate(sub,:),.833);
     projratelo(qqq,:)=quantile(pooledGSLrate(sub,:),.167);
     for www=1:length(colsCONT)
         clear w; w{1}=squeeze(sum(pooledGSLcont(sub,colsCONT{www},:),2));
         [sampsdCONT]=SampleSLRates(w,targyearsGSL,difftimestep);
         projCONT(qqq,:,www)=quantile(w{1},.5);
         projCONThi(qqq,:,www)=quantile(w{1},.833);
         projCONTlo(qqq,:,www)=quantile(w{1},.167);
     end     
end