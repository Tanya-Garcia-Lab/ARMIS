ARMIS
================

## Overview

The R package `ARMIS` analyze the repeated measured data from the same
subject in the small sample size of study. This package uses three mixed
models to capture the time effect, while accounting for the correlations
of measurements over time from the same subject. First, we allow the
flexible variance-covariance structures on the model. Second, we use
baseline measurements as a covariate in the model. Third, we use
percent-change from baseline as a data normalization method. These
methods are described in the following manuscript:

**Lee, U., Garcia, T.P., Carroll, R.J., Gilbreth, K.R., Wu, G. (2019).
“Analysis of repeated measures data in nutrition research”, Frontiers
In Bioscience, Landmark, 24, 1378-1390**

## Installation

To install `ARMIS` from GitHub,

``` r
devtools::install_github("unkyunglee/ARMIS")
```

## Example

``` r
We provide a pseudo dataset to run an example to show reproducibility of our methods in the manuscript. 
We consider the dataset `pseudo_steer_data` available from R package `ARMIS`. 

data(pseudo_steer)
head(pseudo_steer)
```

The data consist of 6 subjects’ information with 4 variables. The amino
acids are repeatedly measured at 6 different time points for each
subject.

``` r
head(pseudo_steer)
?pseudo_steer # this gives you more information on the dataset
```

We fit our methods to the `pseudo_steer_data` data. First, we specify
the parameters and run the function `anova.test()`.

``` r
# Specify the parameters
data=pseudo_steer;
num.aa=1; 
n=6; 
time.points=6;
subid="Steer";
group="Group"; 
time="Time"; 
resp.var="citrulline"; 
amino.names="citrulline";
name.tt=c("time0", "time1", "time2", "time3", "time4", "time5");
name.steer=c("steer1","steer2","steer3","steer4","steer5","steer6");
interv.length=7; 
num.method=4; 
corStruct="gen.ar1"; 
hetero=TRUE;
file="cit_data.csv"

# produce one table and three figures of our manuscript
result<-anova.test(data=speudo_steer, num.aa=1, n=6, time.points=6, subid="Steer",
                   group="Group", time="Time", resp.var="citrulline", amino.names="citrulline",
                   interv.length=7, num.method=4, corStruct="gen.ar1", hetero=TRUE,
                  file="rpaa_cit_data.csv")
```
