# IFILES for LocalizeSL

* **SLRProjections190301core_SEJ.mat**: Core files for RCP8.5+H and 2.0degree+L scenarios in the structured expert judgement of Bamber et al. (2019), 1 March 2019. This file contains two corefile variables, one associated with the RCP8.5+H scenario and one with the 2.0degree+L scenario. Using these corefiles with processing routines requires some careful attention to the scenario specification passed to various functions, as each corefile contains only one scenario and the scenario labels differ from other corefiles.

* **SLRProjections180124GRIDDEDcore_Tscens.mat**: Core file for 1.5°C, 2.0°C and 2.5°C scenarios derived from K14, as developed in Rasmussen et al. (2018), 24 January 2018. Using this corefile with processing routines requires some careful attention to the scenario specification passed to various functions, as the number of scenarios and scenario labels differ from other corefiles.

* **SLRProjections170113GRIDDEDcore.mat**: Core file for gridded K14 projections used in Kopp et al. (2017), 13 January 2017. This files supplants SLRProjections140523core.mat (23 May 2014), used in Kopp et al. (2014).

* **SLRProjections170113GRIDDEDcore-DP16-Pl5_15-BC**: Pre-composed core file for gridded DP16 (with bias correction) projections used in Kopp et al. (2017), 13 January 2017. This can also be regenerated from SLRProjections170113GRIDDEDcore.mat and DecontoPollard-AIS-160411.mat using DecontoPollardEnsembleGSLCompose.

* **SLRProjections161027GRIDDEDcore.mat**: Core file for gridded North American/US territory projections used in Sweet et al. (2017), 27 October 2016

* **DecontoPollard-AIS-160411.mat**: DeConto and Pollard (2016) Antarctic ice sheet projections, as used in Kopp et al. (2016). (Note: not a core file for LocalizeSL)