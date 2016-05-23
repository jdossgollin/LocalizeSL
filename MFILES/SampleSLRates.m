function [sampsdlocrise,targyearrates]=SampleSLRates(sampslocrise,targyears,difftimestep,zerotime)

% [sampsdlocrise,targyearrates]=SampleSLRates(sampslocrise,targyears,difftimestep)
%
% Generate samples of DIFFTIMESTEP-average rate of change from
% SAMPSLOCRISE. 
% 
% OUTPUTS
% -------
% SAMPSDLOCRISE: Cell array of samples of rate of change, same dimensions as sampslocrise
% TARGYEARRATES: Central years of timesteps
% DIFFTIMESTEP: Time step used for rate average (default: 20)
%
% EXAMPLE
% -------
% rootdir='~/Dropbox/Code/LocalizeSL';
% corefile=fullfile(rootdir,'IFILES/SLRProjections140523core.mat');
% selectedSite=180;
% difftimestep=20;
% quantlevs=[.5 .167 .833 .05 .95 .005 .995 .999];
% [sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);
% [sampsdlocrise,targyearrates]=SampleSLRates(sampslocrise,targyears,difftimestep);
% WriteTableSLRProjection(sampsdlocrise,quantlevs,siteids,sitenames,targyearrates,scens,['LSLproj_rates_' num2str(difftimestep) '_'],'mm/yr','%0.1f',1);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 23 13:31:20 EDT 2016

defval('zerotime',2000);
defval('difftimestep',20);

targyears=[zerotime targyears];
Mdiff=(bsxfun(@minus,targyears',targyears)==0)-(bsxfun(@minus,targyears',targyears)==difftimestep);
Mdiff=Mdiff(sum(Mdiff,2)==0,:);
targyearrates=abs(Mdiff)*targyears'/2;
Mdiff=Mdiff/difftimestep;
clear sampsdlocrise;
for www=1:size(sampslocrise,1)
    for jjj=1:size(sampslocrise,2)
        wsamps=sampslocrise{www,jjj};
        wsamps=[zeros(size(wsamps,1),1) wsamps];
        sampsdlocrise{www,jjj}=[Mdiff*wsamps']';
    end   
end