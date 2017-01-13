function WriteTableSLRProjectionAppend(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,filenm,unitstr,fstr,scalefactor)

% WriteTableSLRProjectionAppend(sampslocrise,[quantlevs],siteids,sitenames,targyears,scens,[filenm],[unitstr],[fstr],[scalefactor])
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections
% quantlevs are desired quantiles. (Default = [.01 .05 .167 .5 .833 .95 .99 .995 .999])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jan 13 13:53:32 EST 2017

defval('filenm','LSLproj');
defval('quantlevs',[.01 .05 .167 .5 .833 .95 .99 .995 .999]);
defval('unitstr','cm above 2000 CE level');
defval('fstr','%0.0f');
defval('scalefactor',.1);

if ~exist([filenm '.tsv'],'file')
    fid=fopen([filenm '.tsv'],'w');
    fprintf(fid,[unitstr '\n\n']);
    fprintf(fid,'Site Name\tSite ID\tScenario\tYear');
    fprintf(fid,'\t %0.6g',quantlevs);
    fprintf(fid,'\n');
    fclose(fid);
end
    
for jjj=1:length(siteids)
    fid=fopen([filenm '.tsv'],'a');
    for kkk=1:length(scens)
        for ttt=1:length(targyears)
            if ~isnan(quantile(scalefactor*sampslocrise{jjj,kkk}(:,ttt),quantlevs(1)))
            fprintf(fid,sitenames{jjj});
            fprintf(fid,'\t%0.0f',siteids(jjj));
            fprintf(fid,['\t' scens{kkk}]);
            fprintf(fid,'\t%0.0f',targyears(ttt));
            fprintf(fid,['\t' fstr],quantile(scalefactor*sampslocrise{jjj,kkk}(:,ttt),quantlevs));
            fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
end
