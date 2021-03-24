# COVID19-vaccination

## DOI for this GitHub repository
doi: XXX

## Overview

We investigate relaxation scenarios using an age-structured transmission model that has been fitted to age-specific seroprevalence data, hospital admissions, and projected vaccination coverage for Portugal.

All details are described in the preprint

> XXX

### Data 

The data was added to the ``` data ``` folder for convenience.

### Seroprevalence data

The seroprevalence data has been published at:

> Kislaya I, Gonçalves P, Barreto M, Sousa R, Garcia AC, Matos R, Guiomar R, Rodrigues AP; ISNCOVID-19 Group. Seroprevalence of SARS-CoV-2 Infection in Portugal in May-July 2020: Results of the First National Serological Survey (ISNCOVID-19). Acta Med Port. 2021 Feb 1;34(2):87-94. doi: 10.20344/amp.15122. Epub 2021 Feb 1. PMID: 33641702.

### Contact matrices

The baseline (pre-pandemic) contact matrix is from the publication:

>Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. Inferring high-resolution human mixing patterns for disease modeling. Nat Commun 12, 323 (2021). https://doi.org/10.1038/s41467-020-20544-y

The contact matrix after the first lockdown in April 2020 was inferred using the Dutch matrix from the publication:

> Jantien A. Backer, Liesbeth Mollema, R.A. Eric Vos, Don Klinkenberg, Fiona R.M. van der Klis, Hester E. de Melker, Susan van den Hof, Jacco Wallinga The impact of physical distancing measures against COVID-19 transmission on contacts and mixing patterns in the Netherlands: repeated cross-sectional surveys in 2016/2017, April 2020 and June 2020 medRxiv 2020.05.18.20101501; doi: https://doi.org/10.1101/2020.05.18.20101501

### Demographic data

We used publicly available [data](https://www.pordata.pt/Portugal/Popula%C3%A7%C3%A3o+residente++m%C3%A9dia+anual+total+e+por+grupo+et%C3%A1rio-10) from the Contemporary Portugal Database (PORDATA) : https://www.pordata.pt/.

### Hospitalization data

The hospitalization data are by the Central Administration of the Health System and the Shared Services of the Ministry of Health, covering all public hospitals in Portugal receiving COVID-19 patients.

### Vaccination coverage data 

The vaccination coverage data are by ECDC https://www.ecdc.europa.eu/en/publications-data/data-covid-19-vaccination-eu-eea

## Model

### Inference

Parameter estimation was done with *R version 3.6.0* using *R Studio Version 1.3.1056* (Interface to *R*) and Stan using *rstan* *R* package version 2.21.1 (*R* interface to *Stan*) and *cmdstanr* *R* package Version 0.1.3 on *Windows 10 Home Version 2004*.

The scripts can be found in the ``` scripts ``` directory. The R and Stan scripts are based on scripts used for the publication :

> Rozhnova G, van Dorp CH, Bruijning-Verhagen P, Bootsma MCJ, van de Wijgert JHHM, Bonten MJM, Kretzschmar ME. Model-based evaluation of school- and non-school-related measures to control the COVID-19 pandemic. Nature Communications. 2021;12(1):1614. https://doi.org/10.1038/s41467-021-21899-6

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
- [*Windows 10*] RTools40 https://cran.r-project.org/bin/windows/Rtools/
- [*Windows 10*] Git Bash https://git-scm.com/downloads

If you use *Windows* you probably need *RTools40* and/or *Git Bash* to be able to compile the model.

Some *R package* that also needed to be installed for the code to run are specified in the beginning of the *XXXXXX.R* file.

## Instructions

You should proceed as follows: 1) use the *R* and its packages to fit the model to the data; 2) export the parameter estimates and use the *Mathematica* notebooks to perform the analyses, run scenarios and create the figures.

The necessary files are:

- Age stratified hospitalization data
- Age stratified demography
- Age stratified seroprevalence 
- Baseline (pre-pandemic) contact matrix
- Contact matrix after the lockdown

### R Studio

After all the dependencies are installed in *R* it is necessary to change the directories in the ```XXXXXXXX.R``` file so that the *R Studio* finds all the necessary files. 

Once that is done run the code sequentially from the start until the end.

If *cmdstanr* was installed correctly the line below should create an executable in the working directory. The compilation should take a couple of minutes.

``` sm <- cmdstan_model(stan_model_file) ```

The following line is where the fitting is done. Depending on the precision, number of iterations, etc, it can take from minutes to weeks to complete the fitting.

```
fit <- sm$sample(
  (...)
)
```

To save the fit done run:

```
save(fit.rstan, file = "output/my_fit.rda")
```

To export the parameters found in the fit into a *.csv* file run:

```
write.csv(output, file = "output/my_parameters.csv", row.names = FALSE)
```
After this line of code there are plotting utilities in the *R* script that allow for a quick analysis of the fitting. If desired, they can be skipped and the analysis can be done solely on *Mathematica*.

### Mathematica

To perform the analyses change the directories in the notebook ```XXXXXXX.nb``` such that *Mathematica* finds all necessary files and run the code sequentially from the start until the end.

Depending on the number of parameter samples, complexity of the scenarios, etc the time of computation for the ```XXXXXXX.nb``` notebook is not clearly determined. Most computations run very fast, Re calculation can take up to few hours depending on the number samples used.

***Warning***  : The ```XXXXXXX.nb``` notebook is very RAM hungry and can cause *Mathematica* to crash. It is advised to store the results, clear the definitions and import the results to prevent this from happening. This is already implemented in the code.
