function [fpsitenew]=SubstituteFingerprint(fpsiteold,fp,lo,la,fptosub,fptouse,sitecoords)

% [fpsitenew]=SubstituteFingerprint(fpsiteold,fp,lo,la,fptosub,fptouse,sitecoords)
%
% INPUT:
% fpsiteold - from core file
% fp - gridded fingerprints from which to draw for substitutions
% lo - longitude for gridded fingerprints
% la - latitude for gridded fingerprints
% fptosub - columns of fpsiteold to replace
% fptouse - fingerprints to draw upon for replacement
% sitecoords - coordinates of site
%
% EXAMPLE:
% core=load(corefile)
% [fp,fpname,lo,la] = readFingerprint(subdir)
% substitutep.fpsite=SubstituteFingerprint(core.fpsiteold,fp,lo,la,fptosub,fptouse,core.targsitecoords)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Mar 09 22:44:23 EST 2017

fpsitenew=fpsiteold;

for ii=1:length(fptosub)
    fpsitenew(:,ii) = 1000*interp2(lo,la,fp(:,:,fptouse(ii)),mod(sitecoords(:,2),360),sitecoords(:,1));
    sub=find(sitecoords(:,1)>1e3); fpsite(sub,ii) = 1;
end
