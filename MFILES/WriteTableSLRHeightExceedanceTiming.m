function WriteTableSLRHeightExceedanceTiming(sampslocrise,heights,siteids,sitenames,targyears,scens,oneway,fileprefix) 


% WriteTableSLRHeightExceedanceTiming(sampslocrise,heights,siteids,sitenames,targyears,scens,oneway,[fileprefix])
%
% Write table showing the probability of when heights (in cm) are exceeded.
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections.
%
% If oneway is set to 1 (default), then any samples with sea-level falls
% will be treated as though they stayed at their maximum value. If it is
% set to 0, then they will not be.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-06-29 08:54:03 -0400

defval('fileprefix','LSLheights_');
defval('heights',10:10:500);
defval('oneway',1);


for jjj=1:length(siteids)
    fid=fopen([fileprefix num2str(siteids(jjj)) '.tsv'],'w');
    for kkk=1:length(scens)
        Nsamps=size(sampslocrise{jjj,kkk},1);
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (cm above 2000 CE level)\n']);
        fprintf(fid,'\t %0.1f',heights);
        fprintf(fid,'\n');
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            if oneway
                workslr = max(sampslocrise{jjj,kkk}(:,1:ttt)/10,[],2);
            else
                workslr = sampslocrise{jjj,kkk}(:,ttt)/10;
            end
            u=bsxfun(@gt,workslr,heights(:)');
            u=sum(u,1)/Nsamps;
            fprintf(fid,'\t%0.3f',u);
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');        
    end
    fclose(fid);
end
