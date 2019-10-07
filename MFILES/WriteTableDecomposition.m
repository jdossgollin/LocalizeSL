function WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,fileprefix,subcomp,labls)

% WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens,[fileprefix],[subcomp],[labls])
%
% Writes table with decomposition into major components.
%
% sampsloccomponents, siteids, sitenames, targyears, scens, cols are outputted by LocalizeStoredProjections
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-09-12 18:52:17 -0400

defval('fileprefix','LSLproj_decomp_');
defval('quantlevs',[.01 .05 .167 .5 .833 .95 .99 .995 .999]);

if cols.colOD==cols.colTE
    cols.colTE2=[];
else
    cols.colTE2=cols.colTE;
end

defval('subcomp',{cols.colAIS,cols.colGIS,[cols.colTE2 cols.colOD],[cols.colGIC],[cols.colLS],[cols.colGIA],[1:size(sampsloccomponents{1},2)]});
defval('labls',{'AIS','GIS','Ocean','GIC','LWS','Bkgd','Total'}); 

for jjj=1:length(siteids)
    for kkk=1:length(scens)

        fid=fopen([fileprefix num2str(siteids(jjj)) '_' scens{kkk} '.tsv'],'w');
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (cm above 2000 CE level)\n']);
        fprintf(fid,'\t %0.6g',quantlevs);
        fprintf(fid,'\n');
 
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f\n',targyears(ttt));
            for sss=1:length(labls)
                if length(subcomp{sss})>0
                    u=squeeze(sum(sampsloccomponents{jjj,kkk}(:,subcomp{sss},ttt),2));
                    fprintf(fid,labls{sss});
                    fprintf(fid,'\t%0.1f',quantile(u/10,quantlevs));
                    fprintf(fid,'\n');
                end               
            end
            fprintf(fid,'\n');
        end
    fclose(fid);
    end
end
