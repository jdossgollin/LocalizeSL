function WriteTableSLRRateExceedanceTiming(sampslocrise,rates,siteids,sitenames,targyears,scens,difftimestep,oneway,fileprefix) 

% WriteTableSLRRateExceedanceTiming(sampslocrise,rates,siteids,sitenames,targyears,scens,difftimestep,oneway,[fileprefix])
%
% Write table showing the probability of when rates (in mm/yr) are exceeded.
%
% sampslocrise, siteids, sitenames, targyears, scens are outputted by LocalizeStoredProjections.
%
% If oneway is set to 1, then any samples with sea-level falls
% will be treated as though they stayed at their maximum value.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Apr 06 11:41:11 EDT 2017

defval('fileprefix','LSLrates_exceedance_');
defval('rates',10:10:500);
defval('oneway',0);

defval('timestep',targyears(2)-targyears(1));
initialyear=targyears(1)-timestep;
targyears=[initialyear targyears];

Mdiff=(bsxfun(@minus,targyears',targyears)==0)-(bsxfun(@minus,targyears',targyears)==difftimestep);
Mdiff=Mdiff(sum(Mdiff,2)==0,:);
rateyr=abs(Mdiff)*targyears'/2;
Mdiff=Mdiff/difftimestep;


for jjj=1:length(siteids)
    fid=fopen([fileprefix num2str(siteids(jjj)) '.tsv'],'w');
    for kkk=1:length(scens)
        Nsamps=size(sampslocrise{jjj,kkk},1);
        fprintf(fid,[sitenames{jjj} ' [' num2str(siteids(jjj)) '] -- ' scens{kkk} ' (mm/yr)\n']);
        fprintf(fid,'\t %0.1f',rates);
        fprintf(fid,'\n');
        
        wsamps=sampslocrise{jjj,kkk};
        wsamps=[zeros(size(wsamps,1),1) wsamps];
        samps=Mdiff*wsamps'; samps=samps';

        for ttt=1:length(rateyr)
            fprintf(fid,'%0.0f',rateyr(ttt));
            if oneway
                workslr = max(samps(:,1:ttt),[],2);
            else
                workslr = samps(:,ttt);
            end
            u=bsxfun(@gt,workslr,rates(:)');
            u=sum(u,1)/Nsamps;
            fprintf(fid,'\t%0.3f',u);
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');        
    end
    fclose(fid);
end
