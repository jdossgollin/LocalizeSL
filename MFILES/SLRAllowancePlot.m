function [Ainst,ALDC,ADLfromstart,ADLfp,ADLLDCfromstart,ADLendyears,z0,hp]=SLRAllowancePlot(samps,targyears,effcurve,testz,histcurve,effcurve999,integratecurve,sitelab,params)

% [Ainst,ALDC,ADLfromstart,ADLfp,ADLLDCfromstart,ADLendyears,z0,hp]=
%    SLRAllowancePlot(samps,targyears,effcurve,testz,histcurve,effcurve999,
%    integratecurve,[sitelab],[params])
%
% Output a figure plotting sea level and various flavors of SLR allowances.
% It takes a variety of parameters outputted by SLRFloodNexpVsLevelCurves.
%
% INPUTS:
%
%     samps: matrix of SLR samples (rows: samples, cols: time)
%     targyears: year identifiers for columns on samps
%     effcurve: expected number of floods of different heights (rows: years as
%               in targyears; cols: heights as specified in testz)
%     testz: heights for effcurve
%     histcurve: historical curve for testz (expected)
%     effcurve999: expected number of floods under 99.9th percentile SLR
%     integratecurve: function (curve,t1,t2) to integrate effcurve
%     sitelab: full name of site
%     params: an optional structure with fields setting several parameters
%         - startyear: start year for design-life allowances (default: 2020) 
%         - endyear: end year for allowances (default: 2100)
%         - intperiod: integration period (in yrs) for allowances
%                      with variable start (default: 30)
%         - betas: beta values for Limited Degree of Confidence (LDC) allowances
%                  (default: [0 .5 .67 .9 .95 .99 1])
%         - N0s: expected numbers for allowance calculation; first one specified
%                will be used for the (LDC allowances) (default: [.01 .1 .002])
%         - doplot: which plots to generate; set to 0 for none
%                   (default: 1:6)
%                   1 - SLR; 2 - Instantaneous; 3 - Design-life, fixed start
%                   4 - Design-life, fixed part; 5 - LDC, Instantaneous
%                   6 - LDC, Design-life
%
% OUTPUTS:
%
%     Ainst: Instantaneous allowances (for beta = 1); rows correspond to years
%            and columns to N0s
%     ALDC: Instantaneous allowances with LDC for N0 = N0s(1); rows correspond 
%           to years and columns to betas
%     ADLfromstart: Design-life allowance with fixed starting point
%     ADLfp: Design-life allowance with fixed length
%     ADLLDCfromstart: LDC design-life allowance with fixed starting point
%     ADLendyears: End years for the fixed period design-life allowances
%     z0: height of flood level for each N0 without sea-level rise
%     hp: subplot handles
%
% EXAMPLE:
%
%       sitelab='Boston';
%       filename=['allowance-' sitelab];
%       selectedSite = 235; %PSMSL ID for Boston
%       threshold = 0.6214; % GPD threshold
%       scale = 0.1109; % GPD scale
%       shape = 0.0739; % GPD shape
%       lambda = 2.5699; % Poisson Lambda
%  
%       [sampslocrise,~,siteids,sitenames,targyears,scens,cols] = ...
%                LocalizeStoredProjections(selectedSite,corefile,1);
%       samps=[zeros(size(sampslocrise{1,1},1),1) ...
%              sampslocrise{1,1}]/1000; % add base yr, convert to meters
%       samps=bsxfun(@min,samps,quantile(samps,.999));
%                   % truncate samples viewed as physically implausible
%       targyears = [2000 targyears]; % add base year
%  
%       pm.doplot=0;
%       [effcurve,testz,histcurve,histcurvesamps,effcurveESLR,effcurve999,integratecurve]= ...
%           SLRFloodNexpVsLevelCurves(samps,targyears,threshold, ...
%           scale,shape,lambda,sitelab,pm);
%      
%       clf;
%       [Ainst,ALDC,ADLfromstart,ADLfp,ADLLDCfromstart,ADLendyears,z0,hp]=...
%           SLRAllowancePlot(samps,targyears,effcurve,testz,histcurve,effcurve999, ...
%           integratecurve,sitelab)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Mar 15 17:54:45 EDT 2016

    defval('sitelab','');
    defval('startyear',2020);
    defval('endyear',2100);
    defval('intperiod',30);
    defval('betas',[0 .5 .67 .9 .95 .99 1]);
    defval('N0s',[.01 .1 .002]);
    defval('doplot',1:6);

    if exist('params')
        parseFields(params)
    end


    params.N0s=N0s;
    params.betas=betas;
    params.intperiod=30;
    params.startyear=startyear;
    params.endyear=endyear;
    params.doplot=doplot;
    
    for rrr=1:length(N0s)
        [Ainst(:,rrr),z0(rrr)]=SLRAllowanceFromCurve(testz,histcurve,effcurve,N0s(rrr));
        N0legstr{rrr}=[sprintf('%0.1f',100*N0s(rrr)) '%'];
        
        ADLendyears=[(startyear+1) (startyear+10):10:endyear];
        for ttt=1:length(ADLendyears)
            ADLfromstart(ttt,rrr)=SLRAllowanceFromCurve(testz,histcurve,integratecurve(effcurve,startyear,ADLendyears(ttt)),N0s(rrr));
           end       
        ADLfpyearsubstart=find((targyears+intperiod)<=endyear);
        for ttt=1:length(ADLfpyearsubstart)
            ADLfp(ttt,rrr)=SLRAllowanceFromCurve(testz,histcurve,integratecurve(effcurve,targyears(ADLfpyearsubstart(ttt)),targyears(ADLfpyearsubstart(ttt))+intperiod),N0s(rrr));
        end
        

    end

    for bbb=1:length(betas)
        if betas(bbb)>=0
            weightedcurve=betas(bbb)*effcurve+(1-betas(bbb))*effcurve999;
        else
            weightedcurve=bsxfun(@plus,abs(betas(bbb))*effcurve,(1-abs(betas(bbb)))*histcurve);
        end
        
        [ALDC(:,bbb),z0(rrr)]=SLRAllowanceFromCurve(testz,histcurve,weightedcurve,N0s(1));
        if betas(bbb)>=0
            betalegstr{bbb}=sprintf('\\beta = %0.2f',betas(bbb));
        else
            quote='''';
            betalegstr{bbb}=sprintf(['\\beta' quote ' = %0.2f'],-betas(bbb));
        end
        
        
        for ttt=1:length(ADLendyears)
            ADLLDCfromstart(ttt,bbb)=SLRAllowanceFromCurve(testz,histcurve,integratecurve(weightedcurve,startyear,ADLendyears(ttt)),N0s(1));
            end

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
                % design-life allowances starting in startyear

                plot(ADLendyears,ADLfromstart);
                xlim([startyear endyear]);
                ylabel('m');

                title(['Design life starting ' num2str(startyear)]);
                
                if ~shownNleg
                    legend(N0legstr,'Location','Northwest');
                    shownNleg = 1;
                end
                
            elseif curplot == 4
                plot(targyears(ADLfpyearsubstart),ADLfp);
                xlim([startyear endyear-intperiod]);
                ylabel('m');
                title(['Design life ' num2str(intperiod) ' years']);


                if ~shownNleg
                    legend(N0legstr,'Location','Northwest');
                    shownNleg = 1;
                end

            elseif curplot == 5
                %  LDC instantaneous
                
                clear hl1 hl2;
                sub=find(betas>=0);
                hl1=plot(targyears,ALDC(:,sub));
                hold on;
                sub=find(betas<0);
                hl2=plot(targyears,ALDC(:,sub),'--');
                
                xlim([targyears(1) endyear]);
                ylabel('m');
                title(['Instantaneous (' sprintf('%0.1f',N0s(1)*100) '%)']);

                if ~shownbetaleg
                    hl=legend([hl1 ; hl2],betalegstr,'Location','Northwest');
                    shownbetaleg = 1;
                end
                
            elseif curplot == 6
                % LDC Design-life

                sub=find(betas>=0);
                plot(ADLendyears,ADLLDCfromstart(:,sub));
                hold on;
                sub=find(betas<0);
                plot(ADLendyears,ADLLDCfromstart(:,sub),'--');
                
                xlim([startyear endyear]);
                ylabel('m');
                title(['Design life starting ' num2str(startyear) ' (' sprintf('%0.1f',N0s(1)*100) '%)']);

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