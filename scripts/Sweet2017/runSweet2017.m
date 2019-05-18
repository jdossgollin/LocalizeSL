% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jan 26 16:57:11 EST 2017

addpath(pwd);
runSeaLevelConditionalDistributions
cd ..
runSeaLevelConditionalDistributionUncorrelatedLow
cd ..

addpath(pwd)
workdir='workdir-170111';
if ~exist(workdir,'dir')
    mkdir(workdir)
end
cd(workdir);
runCombineNonclimaticBackground