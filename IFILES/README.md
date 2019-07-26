# Core files for LocalizeSL

![by-nc-sa](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)

Unless otherwise noted, core files and other input files for LocalizeSL are licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0).

* **SLRProjections190726core_SEJ_full.mat**: Core files for RCP8.5+H and 2.0degree+L scenarios in the structured expert judgement of Bamber et al. (2019), 24 July 2019. This version uses a cubic spline to interpolate between the elicited time periods, and updates the previous version (SLRProjections190301core_SEJ.mat), which provided only the elicited time periods. This file contains two corefile variables, one associated with the RCP8.5+H scenario and one with the 2.0degree+L scenario. Using these corefiles with processing routines requires some careful attention to the scenario specification passed to various functions, as each corefile contains only one scenario and the scenario labels differ from other corefiles. For post-2100 RCP8.5+H projections, note that some caution is required, as these combine RCP 8.5 for non-ice sheet components with a 5°C stabilization for ice sheets; the RCP8.5+H scenario may underestimate the ice-sheet contribution under RCP 8.5 (or, conversely, overestimate other components for 5°C stabilization.) This is less of a problem with the 2.0°C stabilization scenario (see Figure 1 of Rasmussen et al., 2018, for the associated temperature trajectories for the non-ice sheet components).

* **SLRProjections180124GRIDDEDcore_Tscens.mat**: Core file for 1.5°C, 2.0°C and 2.5°C scenarios derived from K14, as developed in Rasmussen et al. (2018), 24 January 2018. Using this corefile with processing routines requires some careful attention to the scenario specification passed to various functions, as the number of scenarios and scenario labels differ from other corefiles.

* **SLRProjections170113GRIDDEDcore.mat**: Core file for gridded K14 projections used in Kopp et al. (2017), 13 January 2017. This files supplants SLRProjections140523core.mat (23 May 2014), used in Kopp et al. (2014).

* **SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC**: Pre-composed core file for gridded DP16 (with bias correction) projections used in Kopp et al. (2017), 13 January 2017. This can also be regenerated from SLRProjections170113GRIDDEDcore.mat and DecontoPollard-AIS-160411.mat using DecontoPollardEnsembleGSLCompose.

* **SLRProjections161027GRIDDEDcore.mat**: Core file for gridded North American/US territory projections used in Sweet et al. (2017), 27 October 2016