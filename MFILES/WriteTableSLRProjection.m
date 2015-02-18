function WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,fileprefix)

% WriteTableSLRProjection(sampslocrise,[quantlevs],siteids,sitenames,targyears,scens,[fileprefix])
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections
% quantlevs are desired quantiles. (Default = [.01 .05 .167 .5 .833 .95 .99 .995 .999])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Feb 18 12:41:22 EST 2015

defval('fileprefix','LSLproj_');
defval('quantlevs',[.01 .05 .167 .5 .833 .95 .99 .995 .999]);


for jjj=1:length(siteids)
    fid=fopen([fileprefix num2str(siteids) '.tsv'],'w');
    for kkk=1:length(scens)
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (cm above 2000 CE level)\n']);
        fprintf(fid,'\t %0.6g',quantlevs);
        fprintf(fid,'\n');
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',quantile(sampslocrise{jjj,kkk}(:,ttt)/10,quantlevs));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');        
    end
    fclose(fid);
end
