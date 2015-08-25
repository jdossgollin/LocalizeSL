function [A,z0,lambda]=SLRAllowanceWriteTable(filename,samps,targyears,threshold,scale,shape,lambda,N0s,betas,endyear,sitename)

% [A,z0,lambda]=SLRAllowanceOutput(filename,samps,targyears,threshold
%               scale,shape,lambda,[N0s],[betas],[endyear],[sitename])
%
% Writes a table of instantaneous SLR allowances with various specified
% values of N0 and beta. It outputs annual values, so to calculated
% integrated allowance, you simply take the mean over the years of
% relevance.
%
% INPUTS:
%
%     filename: output filename (automatically appends '.tsv')
%               (default: 'allowances')
%     samps: matrix of SLR samples (rows: samples, cols: time)
%     targyears: year identifiers for columns on samps
%     threshold: GPD threshold
%     scale: GPD scale factor
%     shape: GPD shape factor
%     lambda: mean number of exceedances per year 
%             alternatively, provide a pair of [N ht] and will calculate lambda
%             (e.g., for 1% flood of 2.0 m, [0.01 2.0])
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
%     z0: height of flood level for each N0 without sea-level rise
%     lambda: mean number of exceedances per year for Poisson-GPD
%
% EXAMPLE:
%
%    filename='NYC_allowances';
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
%    [A,z0,lambda]=SLRAllowanceWriteTable(filename,samps, ...
%            targyears,threshold,scale,shape,[0.1 AEP10pt],[],[],[],...
%            'New York City');
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Aug 24 22:15:01 EDT 2015

    defval('filename','allowances');
    defval('sitename',[]);
    defval('endyear',max(targyears));
    defval('betas',[1 .99 .95 .9 .67 .5 0]);
    defval('N0s',[.01 .1 .002]);

    logN = @(z) GPDLogNExceedances(z-threshold,lambda,shape,scale,-threshold); 

    for bbb=1:length(betas)
        for rrr=1:length(N0s)
            [A(:,rrr,bbb),z0(rrr),~,lambda]=SLRAllowance(samps,N0s(rrr),threshold,scale,shape,lambda,betas(bbb),365.25);
        end
    end
    
    %%%
    
    workyears=targyears(1):endyear;
    
    fid = fopen([filename '.tsv'],'w');
    fprintf(fid,'Instantaneous sea-level rise allowances\n');
    if length(sitename)>0; fprintf(fid,[sitename '\n']); end
    fprintf(fid,'\n');
    
    fprintf(fid,'All allowances are instantaneous allowances by year.\n');
    fprintf(fid,'For integrated allowances, take the mean over the years of interest.\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'beta\tN0');
    fprintf(fid,'\t%0.0f',workyears);
    fprintf(fid,'\n');

    for bbb=1:length(betas)
        for rrr=1:length(N0s)
            fprintf(fid,'%0.3f',betas(bbb));
            fprintf(fid,'\t%0.3f',N0s(rrr));
            fprintf(fid,'\t%0.3f',interp1(targyears,A(:,rrr,bbb),workyears));
            fprintf(fid,'\n');
        end
    end
    
    
    fclose(fid);
    
