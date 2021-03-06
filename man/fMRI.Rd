% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fMRI}
\alias{fMRI}
\title{fMRI Data}
\format{A list of two items:
\describe{
  \item{\code{Y}}{Three-way tensor of fMRI activations for 36 distinct trials. For each trial, a time series of 55 fMRI activations were measured at 4,498 spatial locations.}
  \item{\code{map}}{Matrix of (x, y, z) locations for each level of the third mode of \code{Y}}
  ...
}}
\source{
fMRI data for subject 04847, downlaoded from \url{http://www.cs.cmu.edu/afs/cs.cmu.edu/project/theo-81/www/},
containing only time series with 55 time points corresponding to sentence/picture trials as in \url{https://arxiv.org/pdf/1202.2476.pdf}, processed as follows: \cr
\cr
\code{library(R.matlab)} \cr
\code{data <- readMat("data-starplus-04847-v7.mat")} \cr
\code{map <- data$meta[[1]]} \cr
\code{time.55 <- which(unlist(lapply(data$data, function(x) {nrow(x[[1]])})) == 55)} \cr
\code{conds <- numeric(length(data$data))} \cr
\code{for (i in 1:length(conds)) conds[i] <- data$info[, , i]$cond \%in\% c(2, 3)} \cr
\code{cond.23 <- which(conds == 1)} \cr
\code{use <- cond.23[cond.23 \%in\% time.55]} \cr
\code{Y <- array(dim = c(length(use), 55, 4698))} \cr
\code{for (u in use) Y[which(u == use), , ] <- data$data[[u]][[1]]} \cr
\code{fMRI <- list("Y" = Y, "map" = map)}
}
\usage{
\code{data(fMRI)} \cr
\code{Y <- fMRI[["Y"]]} \cr
\code{map <- fMRI[["map"]]}
}
\description{
A dataset containing fMRI measurements for a single subject during 36 distinct trials, during which the subject viewed a picture and a sentence that may or may not correctly describe the picture.
}
\keyword{datasets}
