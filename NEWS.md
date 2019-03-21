# island 0.2.3

* Tests sample data has been hardcoded to avoid problems of RNG versions.
* Citation added

# island 0.2.2

* Duplicated hyperlinks to the vignettes with duplicate titles removed.

# island 0.2.1

* Vignettes island, IBDmodels and detectability added.
* Added functions sss_cedp, mss_cedp and upgma_model_selection that deal with imperfect detectability when data for replicates is present. 
* Added function ibd_models that produce population data for three different models.
* Added functions NLL_rss, NLL_isd, NLL_imd and NLL_env for ready calculation of Negative Log-Likelihoods for specified parameters in each of the sampling schemes treated here.
* Changed and unified output of regular_sampling_scheme and irregular_single_dataset, irregular_multiple_datasets.
* Fixed bugs in PA_simulation, regular_sampling_scheme, irregular_single_dataset, irregular_multiple_datasets and weight_of_evidence.
* Added unitary tests for the package.

# island 0.1.3

## Minor changes:
* Added a Doi for the original source of dataset idaho.
* Changed functions cetotrans, data_generation and PA_simulation to further integrate the environmental covariates framework in the package.
* Corrected some typos.
* Renamed dataset alonso to alonso15 (to reflect that it comes from Alonso's article of 2015).
