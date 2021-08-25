# Bayesian state space model for time series forecast with application in eczema severity prediction

This project is intended to build on [**Hurault et al. (2020), "Personalised prediction of daily eczema severity scores  using a mechanistic machine learning model"**](https://doi.org/10.1111/cea.13717), published in Clinical & Experimental Allergy. The broad structure of the project has made some references to the file structure of Hurault et al. The model files are written in [Stan](https://mc-stan.org/), a probabilistic programming language. Other than that, all files are written in R. R also provides the interface for Stan. To run a Stan model in R, please install the [rstan](https://cran.r-project.org/web/packages/rstan/index.html) package which is like the interface for Stan in R. Dr Reiko Tanaka and Mr Guillem Hurault are my supervisors for this project.

## File structure

[`model_fit.R`](model_fit.R): contains the codes for fitting the six Bayesian state-space models in my project. Indicate the model name and hyperparameters at the beginning of the file and run the remaining codes to estimate the model and obtain the diagnostics. The bulkier models could take around 9 hours to run. Results are stored in `Result` folder.

[`model_check.R`](model_check.R): contains the codes for running prior predictive check, fake data check and posterior predictive check. The model Stan files to be checked should be put in `Model_check` folder. Results are stored in `Result` folder.

[`model_eval.R`](model_eval.R): contains the codes for evaluating models via forward chaining to predict Bother scores. Indicate the model name and hyperparameters at the beginning of the file and run the remaining codes to carry out forward chaining. Results are stored in `Result_eval` folder.

[`model_eval-sassad.R`](model_eval-sassad.R): contains the code for evaluating models that integrated SASSAD via forward chaining to predict SASSAD scores. Results are stored in `Result_eval-sassad` folder.

[`graph_plot.R`](graph_plot.R): contains the codes for plotting the graphs in the report. Please ensure that you have obtained the result files in the respective folders before running this file.

[`function_data.R`](function_data.R): contains the functions for data pre-processing and formatting data into a Stan list to run Stan models.

[`function_eval.R`](function_eval.R): contains the functions for forward chaining and evaluating the results of forward chaining like calculating the metrics and doing posterior checks. Some functions needed for `model_check.R` are also in this file. 

The `Model` folder contains the Stan codes of the models. These files do not have the required generated quantities for checks and predictions. They are for fitting only to speed up fitting slightly. The models are exactly the same in all model files apart from generated quantities specifications.

The `Model_check` folder contains the Stan codes of the models for prior predictive check, fake data check and posterior predictive check. They have the required generated quantities written. 

The `Model_eval` folder contains the Stan codes of the models for forward chaining to predict Bother, with the appropriate generated quantities.

The `Model_eval-sassad` folder contains the Stan codes of the models for forward chaining to predict SASSAD.

The `Model_misc` folder contains the Stan codes of the models that had been attempted but did not produce desired results. It is for documentation purpose.

The `Result` folder stores the results from the checks and fitting the model entirely.

The `Result_eval` folder stores the results from forward chaining to predict Bother.

The `Result_eval-sassad` folder stores the results from forward chaining to predict SASSAD.

The `Study_preliminary` folder contains the codes for the preliminary study stage of the project, it contains 4 files:
 - [`fit_bother.R`](Study_preliminary/fit_bother.R) contains the codes for exploring Bother scores.
 - [`fit_sassad.R`](Study_preliminary/fit_sassad.R) contains the codes for exploring SASSAD scores.
 - [`functions.R`](Study_preliminary/functions.R) stores the functions used in data preprocessing of preliminary analysis
 - [`misc_codes.R`](Study_preliminary/misc_codes.R) contains the other codes written during preliminary study but the results are not directly relevant to stage two of the project.
 
The SWET dataset is loaded from a proprietary package called `TanakaData` which holds the data. Due to the private and sensitive nature of the data, it is not available to the public according to the data sharing agreement.

As this project builds on the foundation of Hurault et al., I have reused or adapted some functions from Hurault et al. The reused functions include:
 - `plot_patient_coef`: uses ggplot2 to plot patients' individual parameter estimates from a fitted model object.
 - `generate_treatment`: this function is for creating fake data of the binary treatment variable. 
 - `compute_pmf`: this functions computes the probability mass of each of the eleven categories of Bother with a posterior sample.

The adapted functions include: 
 - `get_index`: creates a series of "patient, day" pairs.
 - `extract parameters`: extracts the parameters from a fitted Stan object and assign "patient, day" pairs to the parameters for indexing.
 - `prepare_ppc`: prepares a data frame for plotting posterior predictive checks using a fitted Stan object and its parameters.
 - `plot_ppc`: uses ggplot2 to plot the posterior check graph.


