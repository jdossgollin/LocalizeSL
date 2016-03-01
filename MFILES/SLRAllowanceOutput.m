function [Ainst,ALDC,z0,lambda,hp,params]=SLRAllowanceOutput(samps,targyears,threshold,scale,shape,lambda,sitelab,shortname,params)

% [Ainst,ALDC,z0,lambda,hp,params]=SLRAllowanceOutput(samps,targyears,threshold,scale,shape,lambda,[sitelab],[shortname],[params])
%
% Output a figure plotting sea level and various flavors of SLR allowances.
%
% INPUTS:
%
%     samps: matrix of SLR samples (rows: samples, cols: time)
%     targyears: year identifiers for columns on samps
%     threshold: GPD threshold
%     scale: GPD scale factor
%     shape: GPD shape factor
%     lambda: mean number of exceedances per year 
%             alternatively, provide a pair of [N ht] and will calculate lambda
%             (e.g., for 1% flood of 2000 mm, [0.01 2000])
%     sitelab: full name of site
%     shortname: short name of site (for filename; no spaces)
%     params: an optional structure with fields setting several parameters
%         - startyear: start year for integrated allowances (default: 2020) 
%         - endyear: end year for allowances (default: 2100)
%         - intperiod: integration period (in yrs) for allowances
%                      with variable start (default: 30)
%         - betas: beta values for Limited Degree of Confidence (LDC) allowances
%                  (default: [0 .5 .67 .9 .95 .99 1])
%         - N0s: expected numbers for allowance calculation; first one specified
%                will be used for the (LDC allowances) (default: [.01 .1 .002])
%         - doplot: which plots to generate; set to 0 for none
%                   (default: 1:6)
%                   1 - SLR; 2 - Instantaneous; 3 - Integrated, fixed start
%                   4 - Integrated, fixed part; 5 - LDC, Instantaneous
%                   6 - LDC, Integrated
%         - Ainst and ALDC: if passed, will skip calculation of these and use
%                           as specified (if Ainst specified, will not calculate
%                           z0 or lambda)
%
% OUTPUTS:
%
%     Ainst: Instantaneous allowances (for beta = 1); rows correspond to years
%            and columns to N0s
%     ALDC: Instantaneous allowances with LDC for N0 = N0s(1); rows correspond 
%           to years and columns to betas
%     z0: height of flood level for each N0 without sea-level rise
%     lambda: mean number of exceedances per year for Poisson-GPD
%     hp: subplot handles
%     params: parameter values
%
% EXAMPLE:
%
%    shortname='NYC';
%    selectedSite = 12; %PSMSL ID for New York City
%    threshold = 0.5148;
%    scale = 0.1285; % GPD scale
%    shape = 0.1879; % GPD shape
%    AEP10pt = 1.111; % height of 10% AEP flood
%    
%    [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = ...
%             LocalizeStoredProjections(selectedSite,corefile,1);
%    samps=[zeros(size(sampslocrise{1,1},1),1) ...
%           sampslocrise{1,1}]/1000; % add base year and convert to meters
%    samps=bsxfun(@min,samps,quantile(samps,.999));  
%                      % truncate samples viewed as physically implausible
%    targyears = [2000 targyears]; % add base year
%      
%    clf;
%    [Ainst,ALDC,z0,lambda,hp,params]=SLRAllowanceOutput(samps, ...
%            targyears,threshold,scale,shape,[0.1 AEP10pt], ...
%            sitenames{1},shortname);
%
% Note that, to calculate a time series of integrated sea-level rise allowances
% starting at t1 and for years running up to t2 from a vector of instantaneous
% allowaces Ainst(:,j), simply run:
%
%    cumsum(interp1(targyears,Ainst(:,j),t1:t2))./([t1:t2]-t1+1);
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 24 22:16:41 EDT 2015

    defval('siteshortname','SL');
    defval('sitelab',[]);
    defval('startyear',2020);
    defval('endyear',2100);
    defval('intperiod',30);
    defval('betas',[0 .5 .67 .9 .95 .99 1]);
    defval('N0s',[.01 .1 .002]);
    defval('doplot',1:6);
    defval('Ainst',[]);
    defval('ALDC',[]);

    if exist('params')
        parseFields(params)
    end


    params.N0s=N0s;
    params.betas=betas;
    params.intperiod=30;
    params.startyear=startyear;
    params.endyear=endyear;
    params.doplot=doplot;

    logN = @(z) GPDLogNExceedances(z-threshold,lambda,shape,scale,-threshold); 

    if length(Ainst)==0
        for rrr=1:length(N0s)

            [Ainst(:,rrr),z0(rrr),ESLR,lambda,expectedN(:,rrr)]=SLRAllowance(samps,N0s(rrr),threshold,scale,shape,lambda,[],365.25/2);

        end
    else
        z0=[]; lambda=[];
    end

    for rrr=1:length(N0s)
        N0legstr{rrr}=[sprintf('%0.1f',100*N0s(rrr)) '%'];
    end

    if length(ALDC)==0
        for bbb=1:length(betas)
            ALDC(:,bbb)=SLRAllowance(samps,N0s(1),threshold,scale,shape,lambda,betas(bbb),365.25/2);
        end
    end

    for bbb=1:length(betas)
        betalegstr{bbb}=sprintf('\\beta = %0.2f',betas(bbb));
    end


    %%%
    if sum(doplot)>0
        shownbetaleg=0;
        shownNleg=0;
        for counter=1:length(doplot)
            curplot=doplot(counter);
            hp(counter)=subplot(3,2,counter);

            if curplot==1
                
                % sea level plot

                clear projs projs2 projs3;
                projs.x=targyears';
                projs.y=quantile(samps,.5)';
                projs.dy=abs(bsxfun(@minus,quantile(samps,[.167 .833])',projs.y));
                projs2=projs;
                projs2.dy=abs(bsxfun(@minus,quantile(samps,[.05 .95])',projs.y));
                projs3=projs;
                projs3.dy=abs(bsxfun(@minus,quantile(samps,[.01 .99])',projs.y));

                PlotWithShadedErrors(projs3,[1 0 0 ; 0 0 1],.95,'none',':',[2000 2100]);
                hl=PlotWithShadedErrors(projs2,[1 0 0 ;0 0 1],.9,'none','-.',[2000 2100]);
                hl=PlotWithShadedErrors(projs,[1 0 0 ;0 0 1],.7,'-','--',[2000 2100]);
 
                sub=find(projs.x<=2100);
                
               ylim(.2*[floor(5*min(projs3.y(sub)-projs3.dy(sub,1))) ceil(5*max(projs3.y(sub)+projs3.dy(sub,2)))]);
                box on; longticks(gca,2);
                ylabel('Sea level (m)');
                title(sitelab);

            elseif curplot==2
                
                % instantaneous allowances

                plot(targyears,Ainst);
                xlim([targyears(1) endyear]);
                ylabel('m');
                title(['Instantaneous']);

                
                if ~shownNleg
                    legend(N0legstr,'Location','Northwest');
                    shownNleg = 1;
                end
                
            elseif curplot == 3
                % integrated allowances starting in startyeart

                integrateallowanceseries=@(t1,t2,j) cumsum(interp1(targyears,Ainst(:,j),t1:t2))./([t1:t2]-t1+1);
                for jjj=1:length(N0s)
                    plot(2020:2100,integrateallowanceseries(startyear,endyear,jjj)); hold on;
                end                

                ylabel('m');
                title(['Integrated starting ' num2str(startyear)]);
               
                if ~shownNleg
                    legend(N0legstr,'Location','Northwest');
                    shownNleg = 1;
                end
                
            elseif curplot == 4
                % integrated for intperiod years, starting arbitrarily
                M = zeros((endyear-intperiod-startyear+1),(endyear-startyear)+1);
                for ii=1:size(M,1)
                    M(ii,ii:(ii+(intperiod-1)))=1/intperiod;
                end
                wseries=@(t1,t2,j) (interp1(targyears(sub),Ainst(sub,j),t1:t2));

                for jjj=1:length(N0s)
                    plot(startyear:(endyear-intperiod),M*wseries(startyear,endyear,jjj)'); hold on;
                end                
                ylabel('m');
                title(['Integrated ' num2str(intperiod) ' years']);


                if ~shownNleg
                    legend(N0legstr,'Location','Northwest');
                    shownNleg = 1;
                end

            elseif curplot == 5
                %  LDC instantaneous

                plot(targyears,ALDC);
                xlim([targyears(1) endyear]);
                ylabel('m');
                title(['Instantaneous (' sprintf('%0.1f',N0s(1)*100) '%)']);

                if ~shownbetaleg
                    legend(betalegstr,'Location','Northwest');
                    shownbetaleg = 1;
                end
                
            elseif curplot == 6
                % LDC integrated

                integrateallowanceseries=@(t1,t2,b) cumsum(interp1(targyears,ALDC(:,b),t1:t2))./([t1:t2]-t1+1);
                for bbb=1:length(betas)
                    plot(startyear:endyear,integrateallowanceseries(startyear,endyear,bbb)); hold on;
                end                
                ylabel('m');
                title(['Integrated starting ' num2str(startyear) ' (' sprintf('%0.1f',N0s(1)*100) '%)']);

                if ~shownbetaleg
                    legend(betalegstr,'Location','Northwest');
                    shownbetaleg = 1;
                end
            
            end
        end
        
    else
        hp = [];
    end
    
function parseFields(params)

    flds=fieldnames(params);
    for qqq=1:length(flds)
        if flds{qqq} 
            assignin('caller',flds{qqq},params.(flds{qqq}));
        end
    end