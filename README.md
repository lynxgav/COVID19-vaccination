# COVID19-vaccination

## DOI for this GitHub repository
doi: http://doi.org/10.5281/zenodo.4636679

## Overview

We investigate relaxation scenarios during the SARS-CoV-2 vaccination rollout using an age-structured transmission model that has been fitted to age-specific seroprevalence data, hospital admissions, and projected vaccination coverage for Portugal.

All details are described in the preprint that is currently under review by Nature Communications

> João Viana, Christiaan van Dorp, Ana Nunes, Manuel C. Gomes, Michiel van Boven, Mirjam E. Kretzschmar, Marc Veldhoen, and Ganna Rozhnova (2021). Controlling the pandemic during the SARS-CoV-2 vaccination rollout: a modeling study, PREPRINT (Version 1) available at Research Square https://doi.org/10.21203/rs.3.rs-358417/v1.

### Data 

The data was added to the ``` data ``` folder for convenience.

### Seroprevalence data

The seroprevalence data has been published at:

> Kislaya I, Gonçalves P, Barreto M, Sousa R, Garcia AC, Matos R, Guiomar R, Rodrigues AP; ISNCOVID-19 Group. Seroprevalence of SARS-CoV-2 Infection in Portugal in May-July 2020: Results of the First National Serological Survey (ISNCOVID-19). Acta Med Port. 2021 Feb 1;34(2):87-94. doi: 10.20344/amp.15122. Epub 2021 Feb 1. PMID: 33641702.

### Contact matrices

The baseline (pre-pandemic) contact matrix is from the publication:

>Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. Inferring high-resolution human mixing patterns for disease modeling. Nat Commun 12, 323 (2021). https://doi.org/10.1038/s41467-020-20544-y

The contact matrix after the first lockdown in April 2020 was inferred using the Dutch matrix from the publication:

> Backer JA, Mollema L, Vos ER, Klinkenberg D, van der Klis FR, deMelker HE, et al. Impact of physical distancing measures against COVID-19 on contacts and mixing patterns: repeated cross-sectional surveys, the Netherlands, 201617, April 2020 and June 2020. Eurosurveillance. 2021;26(8). doi:https://doi.org/10.2807/1560-7917.ES.2021.26.8.2000994.

### Demographic data

We used publicly available [data](https://www.pordata.pt/Portugal/Popula%C3%A7%C3%A3o+residente++m%C3%A9dia+anual+total+e+por+grupo+et%C3%A1rio-10) from the Contemporary Portugal Database (PORDATA): https://www.pordata.pt/.

### Hospitalization data

The hospitalization data are by the Central Administration of the Health System and the Shared Services of the Ministry of Health, covering all public hospitals in Portugal receiving COVID-19 patients. This set is in the data folder.

### Vaccination coverage data 

The vaccination coverage data for Portugal are by ECDC https://www.ecdc.europa.eu/en/publications-data/data-covid-19-vaccination-eu-eea.

Similar data can be downloaded from [Our World in Data](https://ourworldindata.org/covid-vaccinations) at https://github.com/owid/covid-19-data/tree/master/public/data/vaccinations.

The same data can be downloaded from https://github.com/dssg-pt/covid19pt-data.

## Model

### Inference

Parameter estimation was done with *R version 3.6.0* using *R Studio Version 1.3.1056* (Interface to *R*) and Stan using *rstan* *R* package version 2.21.1 (*R* interface to *Stan*) and *cmdstanr* *R* package Version 0.1.3 on *Windows 10 Home Version 2004*.

The scripts can be found in the ``` scripts ``` directory. The R and Stan scripts are based on scripts used for the publication:

> Rozhnova G, van Dorp CH, Bruijning-Verhagen P, Bootsma MCJ, van de Wijgert JHHM, Bonten MJM, Kretzschmar ME. Model-based evaluation of school- and non-school-related measures to control the COVID-19 pandemic. Nature Communications. 2021;12(1):1614. https://doi.org/10.1038/s41467-021-21899-6.

### Analysis

Analysis of the model was performed on a *Mac OS X El Capitan Version 10.11.5* and *Windows 10 Home Version 2004* using *Mathematica Version Number 10.0.2.0*. 

The notebooks *.nb, where the analyses for different scenarios were performed, can be found in the ```notebooks``` directory.

### Figures

Figures for the manuscript were produced in the notebooks *.nb. Figures can be found in the ```figures``` directory.

### Outputs

Output files produced in *R* or *Mathematica* can be found in the ```outputs``` directory.

## Other

### OS requirements

This package is supported for *macOS*, *Windows 10* and *Linux*. The package has been tested on:

- Windows 10 Home Version 2004
- Mac OS X El Capitan Version 10.11.5

### Hardware requirements

Our study requires only a standard computer with enough RAM to support the in-memory operations.

## Installation guide

Dependencies:

- R Version 3.6.0 https://www.r-project.org/
- R Studio Version 1.3.1056 (Interface to R) https://rstudio.com/
- rstan R package Version 2.21.1 (R interface to Stan) https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html
- cmdstanr R package Version 0.1.3 on Windows 10 Home Version 2004 https://mc-stan.org/cmdstanr/
- Mathematica 10.0.2.0 https://www.wolfram.com/mathematica/

If you use *Windows* you probably need *RTools40* and/or *Git Bash* to be able to compile the model.

- [*Windows 10*] RTools40 https://cran.r-project.org/bin/windows/Rtools/
- [*Windows 10*] Git Bash https://git-scm.com/downloads

Some *R package* that also needed to be installed for the code to run are specified in the beginning of the *XXXXXX.R* file.

## Instructions

You should proceed as follows: 1) use the *R* and its packages to fit the model to the data; 2) export the parameter estimates and use the *Mathematica* notebooks to perform the analyses, run scenarios and create the figures.

The necessary files are:

- Age stratified demography
- Age distribution of morbidities
- Baseline (pre-pandemic) contact matrix
- Contact matrix after the first lockdown
- Age stratified hospitalization data
- Age stratified seroprevalence 
- Vaccination rollout data
- Vaccination plan

### R Studio

After all the dependencies are installed in *R* it is necessary to change the directories in the ```XXXXXXXX.R``` file so that the *R Studio* finds all the necessary files. 

Once that is done run the code sequentially from the start until the end.

If *cmdstanr* was installed correctly the line below should create an executable in the working directory. The compilation should take a couple of minutes.

``` sm <- cmdstan_model(stan_model_file) ```

The following line is where the fitting is done. Depending on the precision, number of iterations, etc, it can take from minutes to days to complete the fitting.

```
fit <- sm$sample(
  (...)
)
```

To save the fit done run:

```
save(fit.rstan, file = "outputs/my_fit.rda")
```

To export the parameters found in the fit into a *.csv* file run:

```
write.csv(output, file = "outputs/my_parameters.csv", row.names = FALSE)
```
After this line of code there are plotting utilities in the *R* script that allow for a quick analysis of the fitting. If desired, they can be skipped and the analysis can be done solely in *Mathematica*.

### Mathematica

To perform the analyses change the directories in the notebooks ```*.nb``` such that *Mathematica* finds all necessary files and run the code sequentially from the start until the end.

Depending on the number of parameter samples, complexity of the scenarios, etc. the time of computation for the notebooks ```*.nb``` is not clearly determined. Most computations run within a matter of seconds but Re calculation can take up to few hours depending on the number samples and operating system used.

***Warning***  : The ```*.nb``` notebook is very RAM hungry on Windows and can cause *Mathematica* to crash. It is advised to store the results, clear the definitions and import the results to prevent this from happening. This is already implemented in the code.
