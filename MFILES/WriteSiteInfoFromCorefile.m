function outtable=WriteSiteInfoFromCorefile(corefile,fn)

% WriteSiteInfoFromCorefile(corefile,fn)
%
% Output site information from corefile.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2020-02-07 14:38:12 -0500

defval('fn',[]);
outtable=sprintf('Site\tSiteID\tLat\tLong\n');

for zzz=1:length(corefile.targregions)
    outtable=[outtable sprintf(corefile.targregionnames{zzz})];
    outtable=[outtable sprintf('\t%0.0f',corefile.targregions(zzz))];
    if isfield(corefile,'targsitecoords')
        outtable=[outtable sprintf('\t%0.3f',corefile.targsitecoords(zzz,:))];
    end
    outtable=[outtable sprintf('\n')];
end

if length(fn)>0
    fid=fopen(fn,'w');
    fprintf(fid,outtable);
    fclose(fid);
end