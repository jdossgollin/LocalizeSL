% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-08-03 21:50:32 -0400

selectedSite = 12;

corefiles = {'SLRProjections170113GRIDDEDcore.mat','SLRProjections170113GRIDDEDcore.mat','SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat','SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC.mat','SLRProjections190726core_SEJ_full.mat','SLRProjections190726core_SEJ_full.mat'}; % specify corefiles to use
corefilelabs = {'K14_8.5', 'K14_2.6','DP16_8.5','DP16_2.6','B19H','B19L'}; % specify corefile labels
subcore = {'', '', '','','corefileH','corefileL'}; % specify if corefile file contains multiple cores
selscens = [1 4 1 4 1 1];
linespecs = {'-', '--', '-.', ':'};
colorspecs = [217 95 2] / 255;
truncatesamplesat = 3.0; % m at which to truncate samples
truncatesamplesyear = 2100; % year to use for truncation

rootdir = '~/Dropbox/Code/LocalizeSL';
addpath(fullfile(rootdir, 'MFILES'));

extremesIFILES = fullfile(rootdir, 'IFILES/extremes');
datDir = fullfile(extremesIFILES, 'declustered');

% import table of calibration parameters
parmdat = importdata(fullfile(extremesIFILES, 'GPDfits_withunc.tsv'), '\t', 1);
extracols = size(parmdat.textdata, 2) - size(parmdat.data, 2);
lambdas = parmdat.data(:, find(strcmpi('lambda', parmdat.textdata(1, :))) - extracols);
thresholds = parmdat.data(:, find(strcmpi('u', parmdat.textdata(1, :))) - extracols);
scales = parmdat.data(:, find(strcmpi('scale', parmdat.textdata(1, :))) - extracols);;
shapes = parmdat.data(:, find(strcmpi('shape', parmdat.textdata(1, :))) - extracols);
AEP10pts = parmdat.data(:, find(strcmpi('AEP0.1', parmdat.textdata(1, :))) - extracols);
Vscale = parmdat.data(:, find(strcmpi('Vscale', parmdat.textdata(1, :))) - extracols);
Vshape = parmdat.data(:, find(strcmpi('Vshape', parmdat.textdata(1, :))) - extracols);
Vscaleshape = parmdat.data(:, find(strcmpi('Vscaleshape', parmdat.textdata(1, :))) - extracols);
psmslids = parmdat.data(:, 1);
tgids = parmdat.data(:, 2);
NOAAnames = parmdat.textdata(2:end, 1);

qqq = find(psmslids == selectedSite);
wfile = fullfile(datDir, ['maxtofit.dclist.' num2str(tgids(qqq)) '_xdat.dc.tsv']);
historicaldata = importdata(wfile);
longname = NOAAnames{qqq};
shortname = longname(setdiff(1:length(longname), strfind(longname, ' ')));
shortname(strfind(shortname, ',')) = '-';

clear targyears effcurve;

for ccc = 1:length(corefiles)

    corefile = load(fullfile(rootdir, ['IFILES/' corefiles{ccc}]));

    if length(subcore{ccc}) > 0
        corefile = corefile.(subcore{ccc});
    end

    ccclab = corefilelabs{ccc};

    % find instantaneous allowances with specified GPD fits

    [sampslocrise, ~, siteids, sitenames, targyears{ccc}, scens, cols] = LocalizeStoredProjections(selectedSite, corefile, selscens(ccc));
    legstr{ccc} = [ccclab '-' scens{1}];

    samps = [zeros(size(sampslocrise{1}, 1), 1) sampslocrise{1, 1}] / 1000; % add base year and convert to meters
    targyears{ccc} = [2000 targyears{ccc}];

    % truncate samples
    [s, si] = sort(samps(:, find(targyears{ccc} == truncatesamplesyear)));
    [mi] = find(s > truncatesamplesat);

    if length(mi) > 0
        subsi = (length(si) - mi(1)):mi(1);
        samps = samps(si(subsi), :);
    end

    acov = [Vscale(qqq) Vscaleshape(qqq); Vscaleshape(qqq) Vshape(qqq)];
    parmsamps = lhsnorm([scales(qqq) shapes(qqq)], acov, 1000);
    parmsamps(:, 1) = max(eps, parmsamps(:, 1));

    clf;
    clear pm;
    pm.showuncertainty = 1; pm.historicaldata = historicaldata;
    [effcurve{ccc}, testz, histcurve, histcurvesamps, effcurveESLR, effcurve999, integratecurve] = ...
        SLRFloodNexpVsLevelCurves(samps, targyears{ccc}, thresholds(qqq), ...
        parmsamps(:, 1), parmsamps(:, 2), lambdas(qqq), longname, pm);

    pdfwrite([shortname '_returncurves_' ccclab]);

    clf;
    clear pm; pm.betas = [0 .5 .9 1];
    [Ainst, ALDC, ADLfromstart, ADLfp, ADLLDCfromstart, ADLendyears, z0, hp] = SLRAllowancePlot(samps, targyears{ccc}, effcurve{ccc}, testz, histcurve, effcurve999, integratecurve, longname, pm);
    pdfwrite([shortname '_Allowances_' ccclab]);

    [A, ADL, z0] = SLRAllowanceWriteTable([shortname '_Allowances'], targyears{ccc}, effcurve{ccc}, testz, histcurve, effcurve999, integratecurve, [.01 .1 .002], [1 .9 .5 .0], [2100], longname)
end

historicalcolor = 'c';

clear hl;
clf;
subplot(2, 1, 1);
ct = sum(bsxfun(@gt, historicaldata, testz)) / (length(historicaldata) / 365.25);
subct = find(diff(ct) < 0);
subct = intersect(subct, find(testz > thresholds(qqq)));
plot(testz(subct), ct(subct), 's', 'Color', historicalcolor); hold on;
hl(1) = plot(testz, histcurve, 'k');
plot(testz, quantile(histcurvesamps, [.17 .83], 1), 'color', [.6 .6 .6]);
hold on;
set(gca, 'yscale', 'log'); ylim([.002 10]);

for ccc = 1:length(corefiles)
    subyr = find(targyears{ccc} == 2050);
    hl(end + 1) = plot(testz, effcurve{ccc}(subyr, :), 'linestyle', linespecs{ccc}, 'color', colorspecs(1, :)); hold on;
end

ylabel({'Expected number of', 'annual exceedances'});
xlabel('Extreme sea level (m)');
title([longname ' - 2050']);
u = legstr';
u = {'Historical', u{1:end}};
hld = legend(hl, u, 'location', 'northeast')
set(hld, 'fontsize', 7)
xlim([0 4.75]);

subplot(2, 1, 2);
ct = sum(bsxfun(@gt, historicaldata, testz)) / (length(historicaldata) / 365.25);
subct = find(diff(ct) < 0);
subct = intersect(subct, find(testz > thresholds(qqq)));
plot(testz(subct), ct(subct), 's', 'Color', historicalcolor); hold on;
hl(1) = plot(testz, histcurve, 'k');
plot(testz, quantile(histcurvesamps, [.17 .83], 1), 'color', [.6 .6 .6]);
hold on;
set(gca, 'yscale', 'log'); ylim([.002 10]);

for ccc = 1:length(corefiles)
    subyr = find(targyears{ccc} == 2100);
    hl(end + 1) = plot(testz, effcurve{ccc}(subyr, :), 'linestyle', linespecs{ccc}, 'color', colorspecs(1, :)); hold on;
end

ylabel({'Expected number of', 'annual exceedances'});
xlabel('Extreme sea level (m)');
title([longname ' - 2100']);
u = legstr';
u = {'Historical', u{1:end}};
hld = legend(hl, u, 'location', 'northeast')
set(hld, 'fontsize', 7)
xlim([0 4.75]);
pdfwrite([shortname '_ESL']);

fid = fopen([shortname '_ESL.tsv'], 'w');
fprintf(fid, '\t%0.3f', testz);
fprintf(fid, '\n');
fprintf(fid, 'historical');
fprintf(fid, '\t%0.3g', histcurve);
fprintf(fid, '\n');
fprintf(fid, 'historical - 17th');
fprintf(fid, '\t%0.3g', quantile(histcurvesamps, .17));
fprintf(fid, '\n');
fprintf(fid, 'historical - 83rd');
fprintf(fid, '\t%0.3g', quantile(histcurvesamps, .83));
fprintf(fid, '\n');

for ccc = 1:length(corefiles)

    doyr = 2050;
    subyr = find(targyears{ccc} == doyr);
    fprintf(fid, [legstr{ccc} ' - %0.0f'], doyr);
    fprintf(fid, '\t%0.3g', effcurve{ccc}(subyr, :));
    fprintf(fid, '\n');

    doyr = 2100;
    subyr = find(targyears{ccc} == doyr);
    fprintf(fid, [legstr{ccc} ' - %0.0f'], doyr);
    fprintf(fid, '\t%0.3g', effcurve{ccc}(subyr, :));
    fprintf(fid, '\n');

end

fclose(fid);
