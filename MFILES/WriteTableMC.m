function WriteTableMC(sampsloccomponents,selcomponents,siteids,sitenames,targyears,scens,fileprefix)

% WriteTableMC(sampsloccomponents,[selcomponents],siteids,sitenames,targyears,scens,[fileprefix])
%
% sampsloccomponents, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections
% Can take sampslocrise in lieu of sampsloccomponents.
%
% If sampsloccomponents are used, selcomponents specifies the columns of the contributing processes
% that are summed up. By default all (1:24) will be summed. A single column can be specified,
% or a set will be added up. For example, 1:23 will sum up columns 1-23 (all except the geological
% background).
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Feb 18 12:45:02 EST 2015

defval('selcomponents',[]);
defval('fileprefix','LSLproj_MC_');


if length(selcomponents)==0
    if ndims(sampsloccomponents{1})==3
        selcomponents=1:size(sampsloccomponents{1},2);
    end 
end

for jjj=1:length(siteids)
    for kkk=1:length(scens)
        if length(selcomponents)==0
            wsamps=sampsloccomponents{jjj,kkk};
        else
            wsamps=squeeze(sum(sampsloccomponents{jjj,kkk}(:,selcomponents,:),2));
        end
        
    fid=fopen([fileprefix num2str(siteids(jjj)) '_' scens{kkk} '.tsv'],'w');
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (cm above 2000 CE level)\n']);
        fprintf(fid,'\n');
        for ttt=1:length(targyears)
            fprintf(fid,'%0.0f',targyears(ttt));
            fprintf(fid,'\t%0.0f',wsamps(:,ttt)/10);
            fprintf(fid,'\n');
        end
    fclose(fid);
    end
end
