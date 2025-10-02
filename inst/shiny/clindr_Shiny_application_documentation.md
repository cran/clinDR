Purpose of this document
========================

The purpose of this document is to provide supporting documentation for
the {clinDR} package Shiny application.

{clinDR} Shiny application
==========================

Code for the {clinDR} Shiny app was originally written by
Mike K Smith, Pfizer (mike.k.smith@pfizer.com), and currently
maintained by Neal Thomas

{clinDR} package
----------------

The {clinDR} package used in the Shiny app was written by Neal Thomas,
Pfizer (snthomas99@gmail.com)

Description
-----------

The Shiny app provides a wrapper around trial simulation functions
within the {clinDR} package.

The package description for {clinDR} is as follows:

    packageDescription("clinDR")

Package: clinDR
Version: 2.4.01
Date: 2022-05-31
Title: Simulation and Analysis Tools for Clinical Dose Response Modeling
Authors@R: c(person("Neal","Thomas",role=c("aut","cre"),
          email="snthomas99@gmail.com",comment=c(ORCID = '0000-0002-1915-8487')),
          person("Jing","Wu",role="aut", email="Jing_Wu@uri.edu"), person("Mike
          K.","Smith",role="aut", email="mike.k.smith@pfizer.com"))
Maintainer: Neal Thomas <snthomas99@gmail.com>
Description: Bayesian and ML Emax model fitting, graphics and simulation for clinical dose
          response.  The summary data from the dose response meta-analyses in Thomas, Sweeney,
          and Somayaji (2014) <doi:10.1080/19466315.2014.924876> and Thomas and Roy (2016)
          <doi:10.1080/19466315.2016.1256229> Wu, Banerjee, Jin, Menon, Martin, and
          Heatherington(2017) <doi:10.1177/0962280216684528> are included in the package.  The
          prior distributions for the Bayesian analyses default to the posterior predictive
          distributions derived from these references.
Depends: R (>= 3.5.0), rstan (>= 2.17.3), shiny
Imports: foreach,graphics,DoseFinding,stats,mvtnorm,utils, parallel,doParallel, ggplot2, tidyr,
          purrr, tibble, dplyr, glue, waiter
License: GPL (>= 2)
NeedsCompilation: no
Packaged: 2021-03-23 20:23:03 UTC; ThomasN
Author: Neal Thomas [aut, cre] (<https://orcid.org/0000-0002-1915-8487>), Jing Wu [aut], Mike K.
          Smith [aut]

Intended scope of use
---------------------

The {clinDR} Shiny app is intended for use in planning of dose-response
studies

### {clinDR} package overview (as it used in the {clinDR} Shiny app)

The application uses the `emaxsim` and `emaxsimB` functions from the
{clinDR} package to assess trial operating characteristics - data is
simulated using a 3- or 4- parameter Emax model with parameters requested by the
user. This simulated data is then analysed using a 3 (or 4
parameter) Emax model. For the maximum likelihood approach implemented in 
emaxsim, if the Emax model fit does not converge,
then the estimation “steps down” to try MLE for 3 parameter Emax model
(if the user specifies trying a 4-parameter Emax), Exponential, Log,
log-linear and linear models. More detail is found in the explanation of
the `emaxsim` and `emaxalt` functions within the {clinDR} package. The
`emaxsimB` function fits a Bayesian Emax model using Stan with 
priors specified by the user through the Shiny app.

Binary data (response rates) can also be simulated and modelled in a
similar fashion.

Further information can be found in {clinDR} documentation.

In particular, users should read the documentation for `emaxsim`,
`emaxsimB` and `emaxPrior.control` carefully to understand how to specify
parameters and priors for the Emax model parameters. These are provided
in tabs within the shiny app for convenience.

### {clinDR} Shiny app

The application provides the user with basic settings for {clinDR}
`emaxsim` and `emaxsimB` simulations. Additional code is used to
graphically summarise the simulation output and fitted results.

The application also provides code used to produce results within the
app so that the user can recreate results in a separate R session by
copy-pasting this code into an R script. The user can then specify
additional arguments to {clinDR} functions to have finer control over the
simulation and analysis processes.

The Shiny app is intended to providing a starting point for colleagues
wishing to evaluate different dose-response or dose-finding designs. It
could provide a means of quick feedback and discussion within MIDD teams
involving Statistics, Clin Pharm, Pharmacometrics and Clinicians. If the
settings within the app are sufficient for a final design decision to be
made, the app provides sufficient reporting and reproducibility with the
user being provided with downloadable code.

Pre-requisites
--------------

This code requires the following packages in addition to {Shiny}:

-   {clinDR}
    -   Need to run the clinDR function compileStanModels() before
        execution.  
-   {glue}
-   {tidyr}  
-   {waiter}
-   {dplyr}
-   {purrr}
-   {tibble}


Version information
-------------------

The version of the {clinDR} Shiny app is v1.4.

### Version history

<table>
<thead>
<tr class="header">
<th style="text-align: left;">Version</th>
<th>Date</th>
<th style="text-align: left;">Comments</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">v1.0</td>
<td>2018/12/03</td>
<td style="text-align: left;">Initial version.</td>
</tr>
<tr class="even">
<td style="text-align: left;">v1.1</td>
<td>2019/01/02</td>
<td style="text-align: left;">Added markdown reporting download functionality</td>
</tr>
<tr class="odd">
<td style="text-align: left;">v1.2</td>
<td>2020/02/20</td>
<td style="text-align: left;">Updated to work with clinDR v2.2.1</td>
</tr>
<tr class="even">
<td style="text-align: left;">v1.3</td>
<td>2020/04/22</td>
<td style="text-align: left;">Refactored</td>
</tr>
<tr class="odd">
<td style="text-align: left;">v1.4</td>
<td>2020/05/19</td>
<td style="text-align: left;">Updates following review by Neal Thomas</td>
</tr>
<tr class="even">
<td style="text-align: left;">v1.5</td>
<td>2020/09/02</td>
<td style="text-align: left;">Incorporate input from reviewers</td>
</tr>
<tr class="odd">
<td style="text-align: left;">v1.6</td>
<td>2020/11/17</td>
<td style="text-align: left;">Handle user feedback if prior settings not reviewed</td>
<tr class="even">
<td style="text-align: left;">v1.7</td>
<td>2021/03/22</td>
<td style="text-align: left;">Minor changes to handling variable selection in
     fit.quantiles and associated glue code for reproducibility</td>
<tr class="odd">
<td style="text-align: left;">v1.8</td>
<td>2021/03/31</td>
<td style="text-align: left;">Revising prior to CRAN submission</td>
</tr>
<tr class="even">
<td style="text-align: left;">v2.1</td>
<td>2023/07/30</td>
<td style="text-align: left;">Major changes to make output updating match inputs</td>
</tr>
</tbody>
</table>


Source code
-----------

Source code for the application is held in a github repository as part of 
the clinDR package.



Testing
-------

### Scope

Since the Shiny app is a wrapper around the {clinDR} package, it is
assumed that the {clinDR} package testing is handled by the package
author, and that this is adequate for its intended use.

### Results

#### Algorithms

The Shiny app has been run with inputs specified in the {clinDR} package
documentation and has been found to match results when running these
examples external to the app.

#### Presented code

When the code generated within the Shiny app is run external to the
application, the results match those presented within the application.
The inputs visibile in the shiny app are checked against the inputs in the
generated code to ensure the inputs visible in the shiny interface
are used in the calculations.


