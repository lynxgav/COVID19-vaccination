# COVID19-vaccination

## DOI for this GitHub repository
doi: XXX

## Overview

We used an age-structured model fitted to age-specific seroprevalence and hospital admission data to investigate the impact of vaccination on controlling the COVID-19 pandemic in Portugal.

All details are described in the preprint

> XXX

### Data 

The data was added to the ``` data ``` folder for convenience.

### Seroprevalence data

The seroprevalence data has been published at:

> Kislaya I, Gonçalves P, Barreto M, Sousa R, Garcia AC, Matos R, Guiomar R, Rodrigues AP; ISNCOVID-19 Group. Seroprevalence of SARS-CoV-2 Infection in Portugal in May-July 2020: Results of the First National Serological Survey (ISNCOVID-19). Acta Med Port. 2021 Feb 1;34(2):87-94. doi: 10.20344/amp.15122. Epub 2021 Feb 1. PMID: 33641702.

### Contact matrices

The unperturbed contact matrix and the school contact matrix are from the publication:

>Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. Inferring high-resolution human mixing patterns for disease modeling. Nat Commun 12, 323 (2021). https://doi.org/10.1038/s41467-020-20544-y

The lockdown contact matrix was extrapolated from the publication:

> Jantien A. Backer, Liesbeth Mollema, R.A. Eric Vos, Don Klinkenberg, Fiona R.M. van der Klis, Hester E. de Melker, Susan van den Hof, Jacco Wallinga The impact of physical distancing measures against COVID-19 transmission on contacts and mixing patterns in the Netherlands: repeated cross-sectional surveys in 2016/2017, April 2020 and June 2020 medRxiv 2020.05.18.20101501; doi: https://doi.org/10.1101/2020.05.18.20101501

### Demographic data

We used publicly available [data](https://www.pordata.pt/Portugal/Popula%C3%A7%C3%A3o+residente++m%C3%A9dia+anual+total+e+por+grupo+et%C3%A1rio-10) from the Contemporary Portugal Database (PORDATA) : https://www.pordata.pt/.

### Hospitalization data

The hospitalization data was shared by a DGS collaborator and is not avaiable for public distribution. ***???***

## Model

### Inference

Model inference was done with *R version 3.6.0* using *R Studio Version 1.3.1056* (Interface to *R*) and Stan using *rstan* *R* package version 2.21.1 (*R* interface to *Stan*) and *cmdstanr* *R* package Version 0.1.3 on *Windows 10 Home Version 2004*.

The scripts can be found in the ``` scripts ``` directory. The R and Stan scripts are based on scripts used for the publication :

> van Boven M, Teirlinck AC, Meijer A, Hooiveld M, van Dorp CH, Reeves RM, Campbell H, van der Hoek W; RESCEU Investigators. Estimating Transmission Parameters for Respiratory Syncytial Virus and Predicting the Impact of Maternal and Pediatric Vaccination. J Infect Dis. 2020 Oct 7;222(Supplement_7):S688-S694. doi: https://doi.org/10.1093/infdis/jiaa424 

### Analysis

Analysis of the model was performed on a *Mac OS X El Capitan Version 10.11.5* and *Windows 10 Home Version 2004* using *Mathematica Version Number 12.1.0.0*. 

The notebook XXXXXXXXX.nb, where the analysis was performed, can be found in the ```notebooks``` directory.

### Figures

Figures for the manuscript were produced in the notebook XXXXXXXXX.nb. Figures can be found in the ```figures``` directory.

### Output

Output files produced in *R* or *Mathematica* can be found in the ```output``` directory.

## Other

### OS requirements

This package is supported for *macOS*, *Windows 10* and *Linux* (altought I didnt test it i think it is possible to run it) . The package has been tested on:

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
- Mathematica 12.1.0.0 https://www.wolfram.com/mathematica/
- [*Windows 10*] RTools40
- [*Windows 10*] Git Bash

If you use *Windows* you probably need *RTools40* and/or *Git Bash* to be able to compile the model.

After all the dependencies are installed it is necessary to change the directories in the *XXXXXXXX.R* file so that the *R Studio* finds all the necessary files. 

Once that is done you only need to sequentially run the code from the start.

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
