# LocalizeSL: Offline sea-level localization code for Kopp et al. (2014)

README file last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Nov 30 22:32:51 EST 2016

## Citation

This code is intended to accompany the results of

	R. E. Kopp, R. M. Horton, C. M. Little, J. X. Mitrovica, M. Oppenheimer,
	D. J. Rasmussen, B. H. Strauss, and C. Tebaldi (2014). Probabilistic 21st
	and 22nd century sea-level projections at a global network of tide	gauge
	sites. Earth's Future 2: 287–306, doi:10.1002/2014EF000239. 

Please cite that paper when using any results generated with this code.

## Overview

This code requires MATLAB with the Statistics Toolbox.

This MATLAB code is intended to help end-users who wish to work with the sea-level rise projections of Kopp et al. (2014) in greater detail than provided by the supplementary tables accompanying that table but without re-running the full global analysis using the supplementary code accompanying the paper. Key functionality these routines provide include:

1. Local sea-level rise projections at decadal time points and arbitrary quantiles
2. Localized Monte Carlo samples, disaggregatable by contributory process
3. Localized variance decomposition plots 

These routines do not provide the extreme flood level analysis in Kopp et al. (2014), but the Monte Carlo time series samples they produce can be combined with other analyses to look at probabilistic changes in flood frequency over time.

The IFILES directory contains the ~200 MB file SLRProjections140523core.mat, which stores 10,000 Monte Carlo samples for each of the processes contributing to global sea-level change, along with metadata. The code loads these samples without regenerating them and then localizes them. (Note that this files is stored in the Github archive using Git Large File Storage; if you do not have git-lfs set up, cloning the archive will only get you a pointer to this file, which you can download using the Github web interface.)

Functions are stored in the MFILES directory.

The most important function is **LocalizeStoredProjections**:

 	[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] =
	LocalizeStoredProjections(focussites,storefile)

LocalizeStoredProjections takes as input two parameters. STOREFILE is the path of the SLRProjections140523core.mat file. FOCUSSITES is the PSMSL ID or IDs of the site(s) of interest. (Please see psmsl.org or the supplementary tables to Kopp et al. (2014) to identify the IDs corresponding to your site of interest. Specify 0 if you want GSL samples returned in the same format.)

LocalizeStoredProjections outputs two M x N cell arrays of localized Monte Carlo samples, SAMPSLOCRISE and SAMPSLOCCOMPONENTS. In each cell array, the m rows correspond to the sites specified in FOCUSSITES and the N columns to different RCPs (specifically, RCP 8.5, RCP 6.0, RCP 4.5, and RCP 2.6). 

The individual cells of SAMPSLOCRISE are P x Q arrays, with the P rows being 10,000 Monte Carlo samples and the Q columns corresponding to decadal time points. The individual cells of SAMPSLOCRISE are P x Q arrays, with the P rows being 10,000 Monte Carlo samples and the Q columns corresponding to decadal time points. The individual cells of SAMPSLOCCOMPONENTS are P x R x Q arrays. The 1st and 3rd dimensions correspond to the rows and columns of SAMPSLOCRISE; the R columns represent 24 different factors contributing to sea-level rise. Specifically, these factors are:

	1 - GIC: Alaska
	2 - GIC: Western Canada/US
	3 - GIC: Arctic Canada North
	4 - GIC: Arctic Canada South
	5 - GIC: Greenland peripheral glaciers
	6 - GIC: Iceland
	7 - GIC: Svalbard
	8 - GIC: Scandinavia
	9 - GIC: Russian Arctic
	10 - GIC: North Asia
	11 - GIC: Central Europe
	12 - GIC: Caucasus
	13 - GIC: Central Asia
	14 - GIC: South Asia West
	15 - GIC: South Asia East
	16 - GIC: Low Latitude
	17 - GIC: Southern Andes
	18 - GIC: New Zealand
	19 - Greenland Ice Sheet
	20 - West Antarctic Ice Sheet
	21 - East Antarctic Ice Sheet
	22 - Land water storage
	23 - Oceanographic processes (thermal expansion and ocean dynamics)
	24 - GIA, tectonics, and other background processes
	
The other outputs of LocalizeStoredProjections are identifying information that can be passed out to the output commands. SITEIDS returns the PSMSL site IDs of selected sites; SITENAMES the names of those sites; TARGYEARS the years of the output; SCENS the RCPs; and COLS are column labels.

Several other provided functions produce output, with detailed parameter specification described in the headers.

**PlotSLRProjection** generates a time series plot analogous to Figure 3 of Kopp et al. (2014).

**PlotSLRProjectionVariance** generates a variance decomposition plot analogous to Figure 4 of Kopp et al. (2014).
**WriteTableMC** outputs Monte Carlo samples.
**WriteTableSLRProjection** outputs desired quantiles of the projections.

## Example usage

	selectedSite = 12; % use PSMSL ID here to select site

	% set up path
	
	rootdir='~/Dropbox/Code/LocalizeSL';
	corefile=fullfile(rootdir,'IFILES/SLRProjections140523core.mat');
	addpath(fullfile(rootdir,'MFILES'));
	
	% generate local samples
	
	[sampslocrise,sampsloccomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(selectedSite,corefile);
	
	% plot curves
	
	figure;
	hp1=PlotSLRProjection(sampslocrise,targyears,[],scens);
	xlim([2000 2100]); ylim([0 200]);
	title(sitenames{1});
	
	% plot variance decomposition
	
	figure;
	hp2=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1);
	subplot(2,2,1); title([ sitenames{1} ' - RCP 8.5']);
	
	figure;
	hp3=PlotSLRProjectionVariance(sampsloccomponents,targyears,cols,[],1,4);
	subplot(2,2,1); title([sitenames{1} ' - RCP 2.6']);
	
	% output quantiles of projections
	
	quantlevs=[.01 .05 .167 .5 .833 .95 .99 .995 .999];
	WriteTableSLRProjection(sampslocrise,quantlevs,siteids,sitenames,targyears,scens);
	
	% output Monte Carlo samples
	WriteTableMC(sampslocrise,[],siteids,sitenames,targyears,scens);
	
	% output Monte Carlo samples without background trend,
	% to allow incorporation of alternative estimates of background trend
	
	WriteTableMC(sampsloccomponents,1:23,siteids,sitenames,targyears,scens,'LSLProj_nobkgd_');
	
	% output decomposition
	WriteTableDecomposition(sampsloccomponents,quantlevs,siteids,sitenames,targyears,cols,scens);

	% pull GSL samples
	[sampsGSLrise,sampsGSLcomponents,siteids,sitenames,targyears,scens,cols] = LocalizeStoredProjections(0,corefile);
	WriteTableDecomposition(sampsGSLcomponents,quantlevs,siteids,sitenames,targyears,cols,scens);

	

----

    Copyright (C) 2015 by Robert E. Kopp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
