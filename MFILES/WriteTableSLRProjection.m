function WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,fileprefix,unitstr,fstr,scalefactor)

% WriteTableSLRProjection(sampslocrise,[quantlevs],siteids,sitenames,targyears,scens,[fileprefix],[unitstr],[fstr],[scalefactor])
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections
% quantlevs are desired quantiles. (Default = [.01 .05 .167 .5 .833 .95 .99 .995 .999])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 23 13:36:34 EDT 2016

defval('fileprefix','LSLproj_');
defval('quantlevs',[.01 .05 .167 .5 .833 .95 .99 .995 .999]);
defval('unitstr','cm above 2000 CE level');
defval('fstr','%0.0f');
defval('scalefactor',.1);

for jjj=1:length(siteids)
    fid=fopen([fileprefix num2str(siteids(jjj)) '.tsv'],'w');
    for kkk=1:length(scens)
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (' unitstr ')\n']);
        fprintf(fid,'\t %0.6g',quantlevs);
        fprintf(fid,'\n');
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,['\t' fstr],quantile(scalefactor*sampslocrise{jjj,kkk}(:,ttt),quantlevs));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');        
    end
    fclose(fid);
end
