function ConditionalDistributionsPlotGSLComponents(p,condtargs,proj,projhi,projlo,projCONT,projCONThi,projCONTlo,colsCONT,colsCONTlab)

% ConditionalDistributionsPlotGSLComponents(p,condtargs,proj,projhi,projlo,projCONT,projCONThi,projCONTlo,colsCONT,colsCONTlab)
%
% Generate output plots and tables for conditional scenarios.
%
% INPUT
% -----
% p: core sea-level structure 
% condtargs: target heights (in mm) upon which conditioned
% proj: median GSL projection for each scenario
% projhi: high GSL projection for each scenario
% projlo: low GSL projection for each scenario
% targyearrates: years for rates
% projCONT: median contribution projections for each scenario
% projCONThi: high contribution projections for each scenario
% projCONTlo: low contribution projections for each scenario
% colsCONT: columns of core files used for contribution breakdown
% colsCONTlab: labels for contribution breakdown
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Oct 01 01:07:49 EDT 2016
%

defval('filesuffix','');

disp('Global Scenarios component table');
fid=fopen(['GSLScenariosComponents' filesuffix '.tsv'],'w');

fprintf(fid,'Global scenarios (cm)\n');
today=date;
fprintf(fid,['Produced by Robert Kopp on ' today '\n\n']);
fprintf(fid,'Site\tID\tLatitude\tLongitude\tScenario\tComponent\t2000');
fprintf(fid,'\t%0.0f',p.targyears);
fprintf(fid,'\n');

for qqq=1:size(projCONT,1)
    fprintf(fid,'GSL');
    fprintf(fid,'\t%0.0f',NaN);
    fprintf(fid,'\t%0.2f',[NaN NaN]);
    fprintf(fid,'\t%0.1f - MED',condtargs(1,qqq)/1000);
    fprintf(fid,'\tTotal');
    fprintf(fid,'\t0');
    fprintf(fid,'\t%0.1f',proj(qqq,:)/10);
    fprintf(fid,'\n');
    
    for ttt=1:length(colsCONT)
        fprintf(fid,'GSL');
        fprintf(fid,'\t%0.0f',NaN);
        fprintf(fid,'\t%0.2f',[NaN NaN]);
        fprintf(fid,'\t%0.1f - MED',condtargs(1,qqq)/1000);
        fprintf(fid,['\t' colsCONTlab{ttt}]);
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projCONT(qqq,:,ttt)/10);
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'GSL');
    fprintf(fid,'\t%0.0f',NaN);
    fprintf(fid,'\t%0.2f',[NaN NaN]);
    fprintf(fid,'\t%0.1f - LOW',condtargs(1,qqq)/1000);
    fprintf(fid,'\tTotal');
    fprintf(fid,'\t0');
    fprintf(fid,'\t%0.1f',projlo(qqq,:)/10);
    fprintf(fid,'\n');
    
    for ttt=1:length(colsCONT)
        fprintf(fid,'GSL');
        fprintf(fid,'\t%0.0f',NaN);
        fprintf(fid,'\t%0.2f',[NaN NaN]);
        fprintf(fid,'\t%0.1f - LOW',condtargs(1,qqq)/1000);
        fprintf(fid,['\t' colsCONTlab{ttt}]);
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projCONTlo(qqq,:,ttt)/10);
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'GSL');
    fprintf(fid,'\t%0.0f',NaN);
    fprintf(fid,'\t%0.2f',[NaN NaN]);
    fprintf(fid,'\t%0.1f - HIGH',condtargs(1,qqq)/1000);
    fprintf(fid,'\tTotal');
    fprintf(fid,'\t0');
    fprintf(fid,'\t%0.1f',projhi(qqq,:)/10);
    fprintf(fid,'\n');
    
    for ttt=1:length(colsCONT)
        fprintf(fid,'GSL');
        fprintf(fid,'\t%0.0f',NaN);
        fprintf(fid,'\t%0.2f',[NaN NaN]);
        fprintf(fid,'\t%0.1f - HI',condtargs(1,qqq)/1000);
        fprintf(fid,['\t' colsCONTlab{ttt}]);
        fprintf(fid,'\t0');
        fprintf(fid,'\t%0.1f',projCONThi(qqq,:,ttt)/10);
        fprintf(fid,'\n');
    end

    

end
fclose(fid);
