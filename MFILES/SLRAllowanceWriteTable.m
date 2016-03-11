function [A,ADL,z0]=SLRAllowanceWriteTable(filename,targyears,effcurve,testz,histcurve,effcurve999,integratecurve,N0s,betas,endyear,sitename)

% [A,ADL,z0]=SLRAllowanceOutput(filename,effcurve,testz,histcurve,effcurve999,
%            [integratecurve],[N0s],[betas],[endyear],[sitename])
%
% Writes a table of instantaneous and design-life SLR allowances with
% various specified values of N0 and beta. It takes a variety of
% parameters outputted by SLRFloodNexpVsLevelCurves.
%
% INPUTS:
%
%     filename: output filename (automatically appends '.tsv')
%               (default: 'allowances')
%     targyears: year identifiers for columns on samps
%     effcurve: expected number of floods of different heights (rows: years as
%               in targyears; cols: heights as specified in testz)
%     testz: heights for effcurve
%     histcurve: historical curve for testz (expected)
%     effcurve999: expected number of floods under 99.9th percentile SLR
%     integratecurve: function (curve,t1,t2) to integrate effcurve
%                     if blank, will not do design-life allowances
%     N0s: expected numbers for allowance calculation
%          (default: [.01 .1 .002])
%     betas: beta values for Limited Degree of Confidence (LDC) allowances
%            (default: [1 .99 .95 .9 .67 .5 0])
%     endyear: end year for allowances (default: max(targyears))
%     sitename: full name of site
%
% OUTPUTS:
%
%     A: Instantaneous allowances; rows correspond to years,
%        columns to N0s, third dimension to betas
%     ADL: Design-life allowances; dimensions are N0s, betas,
%          start years, and end years
%     z0: height of flood level for each N0 without sea-level rise
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
%       [A,ADL,z0]=SLRAllowanceWriteTable(filename,targyears, ...
%            effcurve,testz,histcurve,effcurve999,integratecurve, ...
%            [],[],2100,sitelab);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Mar 11 16:04:42 EST 2016

    defval('filename','allowances');
    defval('sitename','');
    defval('endyear',max(targyears));
    defval('betas',[1 .99 .95 .9 .67 .5 0]);
    defval('N0s',[.01 .1 .002]);
    defval('integratecurve',[]);
    subyears=find(targyears<=endyear);
    
    if length(integratecurve)==0
        doADL=0;
        ADL=[];
    else
        doADL=1;
        ADL=NaN*ones(length(N0s),length(betas),length(subyears),length(subyears));
    end
    
    for bbb=1:length(betas)
        for rrr=1:length(N0s)
            [A(:,rrr,bbb),z0(rrr)]=SLRAllowanceFromCurve(testz,histcurve,betas(bbb)*effcurve+(1-betas(bbb))*effcurve999,N0s(rrr));
            if doADL
                for ttt=1:length(subyears)
                    for sss=(ttt+1):length(subyears)
                        [ADL(rrr,bbb,ttt,sss),z0(rrr)]=SLRAllowanceFromCurve(testz,histcurve,integratecurve(betas(bbb)*effcurve+(1-betas(bbb))*effcurve999,targyears(subyears(ttt)),targyears(subyears(sss))),N0s(rrr));
                    end
                    
                end
            end
            
            
        end
    end
    
    %%%
    
    
    fid = fopen([filename '.tsv'],'w');
    if length(sitename)>0; fprintf(fid,[sitename '\n']); end
    fprintf(fid,'\nInstantaneous sea-level rise allowances\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'beta\tN0');
    fprintf(fid,'\t%0.0f',targyears(subyears));
    fprintf(fid,'\n');

    for bbb=1:length(betas)
        for rrr=1:length(N0s)
            fprintf(fid,'%0.3g',betas(bbb));
            fprintf(fid,'\t%0.3g',N0s(rrr));
            fprintf(fid,'\t%0.2f',A(subyears,rrr,bbb));
            fprintf(fid,'\n');
        end
    end
    
    if doADL
        for bbb=1:length(betas)
            for rrr=1:length(N0s)
                fprintf(fid,'\n\nDesign-life sea-level rise allowances\n');
                fprintf(fid,'Beta:\t%0.3g\tN0:\t%0.3g\n',[betas(bbb) N0s(rrr)]);
                fprintf(fid,'\n');
                
                fprintf(fid,'\t%0.0f',targyears(subyears));
                fprintf(fid,'\n');

                for sss=1:length(subyears)
                    fprintf(fid,'%0.0f',targyears(subyears(sss)));
                    fprintf(fid,'\t%0.2f',squeeze(ADL(rrr,bbb,:,sss)));
                    fprintf(fid,'\n');
                end
                
            end
        end
    end
    
    
    fclose(fid);
    

    