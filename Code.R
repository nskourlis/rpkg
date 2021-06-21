
# creates description and namespace files
usethis::use_description()
usethis::use_namespace()

# Create R directory
base::dir.create("R")

# creates Package-level documentation so you can run ?nameofpackage
usethis::use_package_doc()

# created README.Rmd for Github landing page
# an .Rbuildignore file gets created
usethis::use_readme_rmd()

# creates license file
usethis::use_mit_license("Sahir Bhatnagar")

# creates news file
usethis::use_news_md()

# setup continuous integration via travis-ci
usethis::use_travis()

# sets up testing infrastructure
usethis::use_testthat()



## code to prepare `DATASET` dataset goes here
usethis::use_data(DATASET, overwrite = TRUE)

usethis::use_data_raw()

# load required packages ----
if (!require("pacman")) install.packages("pacman") 
library("pacman")
pacman::p_load(magrittr, dplyr, usethis, data.table, here)
pacman::p_load(here)
### An MSM exapmle ####
library("msm")
library("mstate")
library("tidyverse")
head(cav)
# clean data ----
epil <- cav
DT <- epil %>% as.data.table
DT.base <- DT %>% distinct(PTNUM, .keep_all = TRUE)
DT.base[,`:=`(years=0)]
DT.epil <- rbind(DT, DT.base)
df_epil <-DT.epil

# write data in correct format to data folder ----
usethis::use_data(df_epil, overwrite = TRUE)

#' Seizure Counts for Epileptics
#'
#' @description Thall and Vail (1990) give a data set on two-week 
#' seizure counts for 59 epileptics. The number of seizures was 
#' recorded for a baseline period of 8 weeks, and then patients 
#' were randomly assigned to a treatment group or a control group. 
#' Counts were then recorded for four successive two-week periods. 
#' The subject's age is the only covariate.
#'
#' @format his data frame has 295 rows and the following 5 columns:
#' \describe{
#'   \item{y}{the count for the 2-week period.}
#'   \item{trt}{treatment, "placebo" or "progabide"}
#'   \item{post}{post treatment. 0 for no, 1 for yes}
#'   \item{subject}{subject id}
#'   \item{tj}{time}
#' }
#' @source \url{https://cran.r-project.org/package=MASS}
"df_epil"


pacman::p_load(sinew)
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/msboxes_R.R")
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/freq_total.R")
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/freq_total_msm.R")
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/msmjson2.R")
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/flexjson2.R")
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/mstatejson2.R")
pacman::p_load(sinew)
sinew::makeOxyFile("C:/Users/niksko/Desktop/rpkg/R/runMSMplus.R")

usethis::use_package("survival", type = "Imports")
usethis::use_package("mstate", type = "Imports")
usethis::use_package("msm", type = "Imports")
usethis::use_package("stringi", type = "Imports")
usethis::use_package("RJSONIO", type = "Imports")


usethis::use_package("visNetwork", type = "Imports")
usethis::use_package("shinyjs", type = "Imports")
usethis::use_package("shiny", type = "Imports")
usethis::use_package("mstate", type = "Imports")
usethis::use_package("tidyverse", type = "Imports")
usethis::use_package("tidyr", type = "Imports")
usethis::use_package("ggplot2", type = "Imports")
usethis::use_package("DiagrammeR", type = "Imports")
usethis::use_package("stringr", type = "Imports")
usethis::use_package("dplyr", type = "Imports")
usethis::use_package("gapminder", type = "Imports")
usethis::use_package("plyr", type = "Imports")
usethis::use_package("viridis", type = "Imports")
usethis::use_package("cowplot", type = "Imports")
usethis::use_package("magick", type = "Imports")
usethis::use_package("StatMeasures", type = "Imports")
usethis::use_package("processx", type = "Imports")
usethis::use_package("webshot", type = "Imports")
usethis::use_package("htmlwidgets", type = "Imports")
usethis::use_package("raster", type = "Imports")
usethis::use_package("jsonlite", type = "Imports")
usethis::use_package("devtools", type = "Imports")
usethis::use_package("usethis", type = "Imports")
usethis::use_package("githubinstall", type = "Imports")
usethis::use_package("shinyMatrix", type = "Imports")
usethis::use_package("dlm", type = "Imports")
usethis::use_package("rsvg", type = "Imports")
usethis::use_package("miniUI", type = "Imports")
usethis::use_package("htmltools", type = "Imports")
usethis::use_package("webshot", type = "Imports")



devtools::check()

usethis::use_vignette(name = "MSMplus_application_input")


rpkg::runMSMplus()

