function [hp,vars,fvars]=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,limyrs,sitesel,scensel,dopooledscens,subcomp,labls,colrs,dopanels)

% PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[limyrs],[sitesel],[scensel],[dopooledscens],[subcomp],[labls],[colrs],[dopanels])
%
% Plot variance decomposition of sea-level rise projections.
%
% sampsloccomponents, targyears, cols are outputted by LocalizeStoredProjections.
% limyrs are the limts of years plotted (default: [2010 2100])
% sitesel is the row id of the site of interest in sampslocrise (default: 1)
% scensel is the sequential id of scenario of interest (default: 1)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Jan 14 17:16:36 EST 2017

% variance plots

defval('colrs','rcbgmyrk');
linew=[1 1 1 1 1 1 1];

defval('dopooledscens',0);
defval('sitesel',1); % scenario selection
defval('scensel',1); % site selection
defval('limyrs',[2010 2100]);
defval('dopanels',[1 2]);

vars=[];
fvars=[];

subyears=find((targyears>=limyrs(1)).*(targyears<=limyrs(2)));

if cols.colOD==cols.colTE
    cols.colTE2=[];
else
    cols.colTE2=cols.colTE;
end

cols.colIS = [cols.colAIS cols.colGIS];
defval('subcomp',{cols.colAIS,cols.colIS,[cols.colIS cols.colTE2 cols.colOD],[cols.colGIC cols.colIS cols.colTE2 cols.colOD],[cols.colGIC cols.colIS cols.colLS cols.colTE2 cols.colOD],[cols.colGIC cols.colIS cols.colLS cols.colTE2 cols.colOD cols.colGIA]});
defval('labls',{'AIS','GIS','Ocean','GIC','LWS','Bkgd'});

% add variance if a field
if isfield(cols,'colvar')
    labls={'Var',labls{:}};
    subcomp={cols.colvar,subcomp{:}};
    for qq=2:length(subcomp)
        subcomp{qq}=[subcomp{qq} cols.colvar];
    end
end


clear hl hp;
if ismember(1,dopanels)
    if length(dopanels)>1
        hp(1)=subplot(2,2,1);
    else
        hp=gca;
    end

    yrs=[targyears(subyears)];
    vlast=zeros(1,length(yrs));

    for i=1:length(labls)		
        u=squeeze(sum(sampsloccomponents{sitesel,scensel}(:,subcomp{i},subyears),2));
        vcur=[var(u)]/1e6;
        vcur(vcur<(1e-9*max(vcur)))=0;
        hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
        vlast=vcur;
        vars(i,:)=vcur;
    end

    if dopooledscens
        u0A = [];
        for scensel2=1:size(sampsloccomponents,2)
            u0A=[u0A ; squeeze(sum(sampsloccomponents{sitesel,scensel2}(:,subcomp{end},subyears),2))];
        end

        for jj=1:size(u0A,2)
            sub=find(~isnan(u0A(:,jj)));
            varu0A(jj)=var(u0A(sub,jj));
        end
        vcur=[varu0A]/1e6;
        hl(i+1)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i+1)); hold on;
        vars(end+1,:)=vcur;
    end

    ylabel('m^2');
    if dopooledscens
        legend(hl(end:-1:1),'Scen',labls{end:-1:1},'Location','Northwest');
    else
        legend(hl(end:-1:1),labls{end:-1:1},'Location','Northwest');
    end
    longticks(gca);
end

if ismember(2,dopanels)
    if length(dopanels)>1
        hp(2)=subplot(2,2,2);
    else
        hp=gca;
    end
    yrs=[targyears(subyears)];
    vlast=zeros(1,length(yrs));
    clear hl;

    if ~dopooledscens
        varu0A=var(squeeze(sum(sampsloccomponents{sitesel,scensel}(:,subcomp{end},subyears),2)));
    end

    for i=1:length(labls)		
        u=squeeze(sum(sampsloccomponents{sitesel,scensel}(:,subcomp{i},subyears),2));
        vcur=[var(u)./varu0A];
        vcur(vcur<(1e-9*max(vcur)))=0;
        hl(i)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i)); hold on;
        vlast=vcur;
        fvars(i,:)=vcur;
    end
    if dopooledscens
        vcur=[varu0A./varu0A];
        hl(i+1)=patch([yrs yrs(end:-1:1)],[vcur vlast(end:-1:1)],colrs(i+1)); hold on;
        fvars(end+1,:)=vcur;

    end
    ylabel('Fraction of variance');
    ylim([0 1]);
    if length(dopanels)==1
        if dopooledscens
            legend(hl(end:-1:1),'Scen',labls{end:-1:1},'Location','Northwest');
        else
            legend(hl(end:-1:1),labls{end:-1:1},'Location','Northwest');
        end
    end
    longticks(gca);
    set(hp,'xlim',limyrs);

end