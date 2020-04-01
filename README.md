
# Mortality-dynamics

by Katrien Antonio, Bavo Campo and Sander Devriendt.

Course materials for the *Advanced Life Insurance Mathematics* course
taught in academic year 2019-2020 at KU Leuven.

📆 March - June, 2020 <br> 🕜 approx. 2-3h per computer lab <br> 📌
Advanced Life Insurance Mathematics class at KU Leuven

## Goals of the computer labs

You’ll work through the essential steps of the implementation in `R` of
the concepts discussed in the Advanced Life Insurance Mathematics
course.

## Schedule and Course Material

The schedule is subject to small changes.

| Description                                                 | Lecture material                                                                                                        | R script                                                                                    | R solutions |
| ----------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | ----------- |
| Computer lab 1 - Getting mortality data                     | [sheets data](https://katrienantonio.github.io/mortality-dynamics/sheets/ALIM_computer_lab_1.html#data)                 | [script](https://katrienantonio.github.io/mortality-dynamics/scripts/ALIM_computer_lab_1.R) |             |
| Computer lab 1 - Parametric mortality laws                  | [sheets mortality laws](https://katrienantonio.github.io/mortality-dynamics/sheets/ALIM_computer_lab_1.html#parametric) | [script](https://katrienantonio.github.io/mortality-dynamics/scripts/ALIM_computer_lab_1.R) |             |
| Computer lab 2 - Fitting Lee Carter with iterative LS       | [sheets iterative LS](https://katrienantonio.github.io/mortality-dynamics/sheets/ALIM_computer_lab_2.html#LS)           | [script](https://katrienantonio.github.io/mortality-dynamics/scripts/ALIM_computer_lab_2.R) |             |
| Computer lab 2 - Fitting Lee Carter with Poisson likelihood | [sheets Poisson](https://katrienantonio.github.io/mortality-dynamics/sheets/ALIM_computer_lab_2.html#POI)               | [script](https://katrienantonio.github.io/mortality-dynamics/scripts/ALIM_computer_lab_2.R) |             |
| Computer lab 2 - Forecasting with Lee Carter                | [sheets forecasting](https://katrienantonio.github.io/mortality-dynamics/sheets/ALIM_computer_lab_2.html#forecasting)   | [script](https://katrienantonio.github.io/mortality-dynamics/scripts/ALIM_computer_lab_2.R) |             |

## Software requirements

Please bring a laptop with a recent version of R and RStudio installed.
Make sure you can connect your laptop to the internet (or download the
course material one day before the start of the workshop). You will
need:

  - R (at least 3.5.2 <https://cloud.r-project.org/bin/windows/base/> )
  - RStudio (
    <https://www.rstudio.com/products/rstudio/download/#download> )

You should install and load the packages that will be used throughout
the workshop. You can use the following instructions to install (if
necessary) and load the packages. These instructions are also available
in `prework_installation_packages.R` from the `scripts` folder.

``` r
packages <- c("tidyverse", "demography")
suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))
```

## RStudio Cloud

We provide a back-up plan in case your local installation of R (and the
necessary packages) is not working properly. Please join our project on
RStudio Cloud via the following link

<https://rstudio.cloud/spaces/55584/join?access_code=hnHmPmSNsbu8BdEBJkeMXNPtHBelpnnJBAmF9XTH>

After creating an account for RStudio you will be able to work with the
scripts and data sets in the cloud.

## Instructor

<img src="img/Katrien.jpg" width="110"/>

<p align="justify">

[Katrien Antonio](https://katrienantonio.github.io/) is professor in
insurance data science at KU Leuven and associate professor at
University of Amsterdam. She teaches courses on data science for
insurance, life and non-life insurance mathematics and loss models.
Research-wise Katrien puts focus on pricing, reserving and fraud
analytics, as well as mortality dynamics.

## Let’s go\!

You are now ready to load the data and build predictive models for life
insurance applications.
