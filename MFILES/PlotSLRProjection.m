function hp=PlotSLRProjection(sampslocrise,targyears,sitesel,scens,subscens)

% PlotSLRProjection(sampslocrise,targyears,[sitesel],[scens],[subscens])
%
% Plot sea-level rise projections.
%
% sampslocrise, targyears, and scens are outputted by LocalizeStoredProjections.
% sitesel is the row id of the site of interest in sampslocrise (default: 1)
% subscens can be used to select scenarios by sequential id (default: [1 3 4])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Feb 18 12:31:09 EST 2015

%%%

defval('sitesel',1);
defval('scens',{'RCP85','RCP60','RCP45','RCP26'});
defval('subscens',[1 3 4]);
colrs='rcbg';

scenlab={};
for kk=subscens
	scenlab={scenlab{:},[upper(scens{kk}(1:3)) ' ' scens{kk}(4) '.' scens{kk}(5)]};
end

for kk=subscens
    
    hp(kk)=plot([2000 targyears],[0  quantile(sampslocrise{sitesel,kk}/10,.5)],[colrs(kk)],'linew',2); hold on;
    plot([2000 targyears],[0  quantile(sampslocrise{sitesel,kk}/10,.05)],[colrs(kk),'--']);
    plot([2000 targyears],[0  quantile(sampslocrise{sitesel,kk}/10,.95)],[colrs(kk),'--']);
    plot([2000 targyears],[0  quantile(sampslocrise{sitesel,kk}/10,.01)],[colrs(kk),':']);
    plot([2000 targyears],[0 quantile(sampslocrise{sitesel,kk}/10,.99)],[colrs(kk),':']);
end
ylabel('cm'); ylim([0 300]); xlim([2000 targyears(end)]);
legend(hp(subscens),scenlab{:},'location','NorthWest');