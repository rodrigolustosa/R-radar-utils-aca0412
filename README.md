# R-radar-utils-aca0412

This repository is a partial translation of the IDL scripts used in the Radar Meteorology (ACA0412) course in 2022 at the Institute of Astronomy, Geophysics and Atmospheric Sciences of the University of São Paulo (IAG/USP) taught by professor Carlos Morales. 

## Introduction

There is not a known R package with all tools needed to manipulate weather radar data, from common mathematical and geographical operations to make plots, as opposed to Python language (e.g.: [Py-ART](https://arm-doe.github.io/pyart/)). This repository aims to translate some IDL scripts so that R language can be used in the Radar Meteorology (ACA0412) course at University of São Paulo, and maybe serve as a base to make a complete weather radar R package in the future. It was made as a final project from the same course.


## Packages
Some R packages were used in the scripts. You can install them executing the following code in your R console:
```
install.packages("tidyverse")
install.packages("fields")
install.packages("lubridate")
install.packages("geobr")
install.packages("shiny")
```
If there is any problem with `geobr`, you can install it with `devtools::install_github("ipeaGIT/geobr", subdir = "r-package")` in place of `install.packages("geobr")`. And if you are using Linux, some Linux packages should also be installed ([as explainned here](https://blog.zenggyu.com/en/post/2018-01-29/installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/)). Running the following code in your terminal should install all the ones you need. 
```
sudo apt install r-base-dev
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
```
