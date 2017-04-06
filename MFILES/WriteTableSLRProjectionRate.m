function WriteTableSLRProjectionRate(sampslocrise,quantlevs,siteids,sitenames,targyears,scens,difftimestep,fileprefix,unitstr,fstr)

% WriteTableSLRProjectionRate(sampslocrise,[quantlevs],siteids,sitenames,targyears,scens,difftimestep,[fileprefix],[unitstr],[fstr])
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections
% difftimestep is the time step of the rate calculation
% quantlevs are desired quantiles. (Default = [.01 .05 .167 .5 .833 .95 .99 .995 .999])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Apr 12 16:37:37 EDT 2016

defval('fileprefix','LSLproj_rate_');
defval('quantlevs',[.01 .05 .167 .5 .833 .95 .99 .995 .999]);
defval('unitstr','mm/yr above 2000 CE level');
defval('fstr','%0.2f');

defval('timestep',targyears(2)-targyears(1));
initialyear=targyears(1)-timestep;
targyears=[initialyear targyears];

Mdiff=(bsxfun(@minus,targyears',targyears)==0)-(bsxfun(@minus,targyears',targyears)==difftimestep);
Mdiff=Mdiff(sum(Mdiff,2)==0,:);
rateyr=abs(Mdiff)*targyears'/2;
Mdiff=Mdiff/difftimestep;

for jjj=1:length(siteids)
    fid=fopen([fileprefix '-' num2str(difftimestep) 'y-' num2str(siteids(jjj)) '.tsv'],'w');
    for kkk=1:length(scens)
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (' unitstr ')\n']);
        fprintf(fid,'%0.0f-year average rate\n',difftimestep);
        fprintf(fid,'\t %0.6g',quantlevs);
        fprintf(fid,'\n');
        wsamps=sampslocrise{jjj,kkk};
        wsamps=[zeros(size(wsamps,1),1) wsamps];
        samps=Mdiff*wsamps'; samps=samps';
        for ttt=1:length(rateyr)
            fprintf(fid,'%0.0f',rateyr(ttt));
            fprintf(fid,['\t' fstr],quantile(samps(:,ttt),quantlevs));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');        
    end
    fclose(fid);
end
