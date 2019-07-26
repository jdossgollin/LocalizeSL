%% Example script for use of LocalizeSL
% with Bamber et al. (2019) projections.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-03-01 10:58:34 -0500

selectedSite = 0; % use PSMSL ID here to select site

% set up path

rootdir='~/Dropbox/Code/LocalizeSL';
corefiles=load(fullfile(rootdir,'IFILES/SLRProjections190725core_SEJ_full.mat'));
addpath(fullfile(rootdir,'MFILES'));

% specify scenario labels and scenarios to use
% important since we differ from defaults here

for ccc=1:2
    if ccc==1
        corefile=corefiles.corefileH;
        selscens=1;
        scenlabs={'RCP85+H'};
    elseif ccc==2
        corefile=corefiles.corefileL;
        selscens=1;
        scenlabs={'2.0C+L'};
    end

    % generate local samples

    [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,selscens);

    % plot curves

    figure;
    hp1=PlotSLRProjection(sampslocrise,targyears,[],scenlabs,selscens,2300);
    title(sitenames{1});
    pdfwrite(['GSLproj_' scenlabs{1}])

    % plot variance decomposition

    figure;
    hp2=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[2050 2300],1,1);
    subplot(2,2,1); title([ sitenames{1} ' - ' scenlabs{1}]);
    pdfwrite(['GSLprojvar_' scenlabs{1}])

    % output quantiles of projections

    quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
    WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,['LSLproj_' scenlabs{1}  '_' ]);

    % output timing of height exceedances
    WriteTableSLRHeightExceedanceTiming(sampslocrise,[],siteids,sitenames,targyears,scens,1,['LSLheights_' scenlabs{1}  '_' ]);

    % output Monte Carlo samples
    WriteTableMC(sampslocrise,[],siteids,sitenames,targyears,scens,['LSLproj_MC_' scenlabs{1}  '_' ]);

    % output decomposition
    WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,['LSLproj_decomp_' scenlabs{1}  '_' ]);

end