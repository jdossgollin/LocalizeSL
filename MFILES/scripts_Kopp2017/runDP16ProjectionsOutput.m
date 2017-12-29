% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-12-29 15:15:20 -0500

% now output DP16 projections, for various combinations
% of Pliocene filters and bias correction 

for jjj=1:size(sampsDP,1);
    for kkk=1:size(sampsDP,2)
        clear substituteDP; substituteDP.samps = sampsDP{jjj,kkk};
        [sampsGSLrise2s{jjj,kkk},sampsGSLcomponents2s{jjj,kkk}]=LocalizeStoredProjections(0,corefile,RDscenmap,substituteDP);
    end
end

for preferredEnsemble=1:length(ensembleLab)
    for preferredBC=1:length(bcsets)
        preferredEnsembleSetLabel=[ensembleLab{preferredEnsemble} bcsets{preferredBC}];

        % output revised GSL

        WriteTableDecomposition(sampsGSLcomponents2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,colsGSL,scensGSL,['GSLproj_decomp_' preferredEnsembleSetLabel '_']);
        WriteTableSLRProjection(sampsGSLrise2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,['GSLproj_' preferredEnsembleSetLabel '_']);

        %%% now generate local projections

        fn=['LSLproj_' preferredEnsembleSetLabel ];
        if exist([fn '.tsv'],'file');
            delete([fn '.tsv']);
        end

        WriteTableSLRProjectionAppend(sampsGSLrise2s{preferredBC,preferredEnsemble},quantlevs,siteidsGSL,sitenamesGSL,targyears,scensGSL,fn);
        
        for ppp=1:length(selectedSites)

            selectedSite=selectedSites(ppp);
            
            % generate local samples
            substituteDP.samps = sampsDP{preferredBC,preferredEnsemble};
            [sampslocrise2,sampsloccomponents2,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile,[],substituteDP);

            % output quantiles of projections
            quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
            WriteTableSLRProjectionAppend(sampslocrise2,quantlevs,siteids,sitenames,targyears,scens,fn);

        end
    end
end