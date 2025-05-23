% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conduct_CPT.R
\name{conduct_MCCPT}
\alias{conduct_MCCPT}
\title{Run a Monte-Carlo changepoints analysis on supplied paleoclimate record data}
\usage{
conduct_MCCPT(
  sites_data,
  age_lowerbound = 6000,
  age_upperbound = 22000,
  output_folder = paste0(getwd(), "/CPT_outputs/"),
  minseg_len = NULL,
  n_cpts = NULL,
  cpt_calc = "meanvar",
  cpt_method = "BinSeg",
  penalty_value = "2*log(n)",
  penalty_type = "Manual",
  save = TRUE,
  rtn = TRUE,
  rev_y = FALSE,
  verbose = TRUE,
  uncertainty_res = 20
)
}
\arguments{
\item{sites_data}{formatted paleoclimate data, output from import_files.}

\item{age_lowerbound}{minimum age for the analysis in years BP.}

\item{age_upperbound}{maximum age for the analysis in years BP.}

\item{output_folder}{full file path to a folder where the outputs will be save.}

\item{minseg_len}{NULL or specify a minimum segment length. User will be prompted for each record if none is specified.}

\item{n_cpts}{number of changepoints to attempt. User will be prompted on a per-record basis if none is specified.}

\item{cpt_calc}{which type of changepoint calculation? One of "meanvar" or "mean".}

\item{cpt_method}{passed to method in cpt.mean or cpt.meanvar. See \link[changepoint]{cpt.mean}}

\item{penalty_value}{value used for pen.value in cpt.mean or cpt.meanvar.}

\item{penalty_type}{character; value passed to penalty in cpt.mean or cpt.meanvar.}

\item{save}{TRUE/FALSE to save outputs to the output_folder.}

\item{rtn}{TRUE/FALSE to return changepoint data.}

\item{rev_y}{TRUE/FALSE to reverse y axis on plots}

\item{uncertainty_res}{Resolution in years for uncertainty resolution within the PrC calculation. Adjust this for shorter series.}
}
\description{
Apply MCCPT to paleoclimate data
}
