#' Brain Tumor Data
#'
#' A dataset containing gene expression data for 43 tumors
#'
#' @usage
#' \code{data(braintumors)} \cr
#' \code{Y <- braintumors[["Y"]]} \cr
#' \code{class <- braintumors[["class"]]}
#'
#' @format A list of two items:
#' \describe{
#'   \item{\code{Y}}{Matrix of 356 continuous gene expression measurements for 43 brain tumors.}
#'   \item{\code{class}}{Vector giving tumor type for each tumor represented in \code{Y}}
#'   ...
#' }
#' @source Data constructed as follows from the \code{denoiseR} package (\url{https://github.com/julierennes/denoiseR}): \cr
#' \code{library(devtools)} \cr
#' \code{install_github("julierennes/denoiseR")} \cr
#' \code{library(denoiseR)} \cr
#' \code{data(tumors)} \cr
#' \code{Y <- as.matrix(tumors[order(tumors[, ncol(tumors)]), -ncol(tumors)])} \cr
#' \code{class <- tumors[order(tumors[, ncol(tumors)]), ncol(tumors)]} \cr
#' \code{braintumors <- list("Y" = Y, "class" = class)}
"braintumors"
#' fMRI Data
#'
#' A dataset containing fMRI measurements for a single subject during 36 distinct trials, during which the subject viewed a picture and a sentence that may or may not correctly describe the picture.
#'
#' @usage
#' \code{data(fMRI)} \cr
#' \code{Y <- fMRI[["Y"]]} \cr
#' \code{map <- fMRI[["map"]]}
#'
#' @format A list of two items:
#' \describe{
#'   \item{\code{Y}}{Three-way tensor of fMRI activations for 36 distinct trials. For each trial, a time series of 55 fMRI activations were measured at 4,498 spatial locations.}
#'   \item{\code{map}}{Matrix of (x, y, z) locations for each level of the third mode of \code{Y}}
#'   ...
#' }
#' @source fMRI data for subject 04847, downlaoded from \url{http://www.cs.cmu.edu/afs/cs.cmu.edu/project/theo-81/www/},
#'containing only time series with 55 time points corresponding to sentence/picture trials as in \url{https://arxiv.org/pdf/1202.2476.pdf}, processed as follows: \cr
#' \cr
#' \code{library(R.matlab)} \cr
#' \code{data <- readMat("data-starplus-04847-v7.mat")} \cr
#' \code{map <- data$meta[[1]]} \cr
#' \code{time.55 <- which(unlist(lapply(data$data, function(x) {nrow(x[[1]])})) == 55)} \cr
#' \code{conds <- numeric(length(data$data))} \cr
#' \code{for (i in 1:length(conds)) conds[i] <- data$info[, , i]$cond \%in\% c(2, 3)} \cr
#' \code{cond.23 <- which(conds == 1)} \cr
#' \code{use <- cond.23[cond.23 \%in\% time.55]} \cr
#' \code{Y <- array(dim = c(length(use), 55, 4698))} \cr
#' \code{for (u in use) Y[which(u == use), , ] <- data$data[[u]][[1]]} \cr
#' \code{fMRI <- list("Y" = Y, "map" = map)}
"fMRI"
#' Fusarium Data
#'
#' A dataset containing severity of disease incidence ratings for 20 varieties of wheat infected with 7 strains of Fusarium head blight over 4 years, 1990-1993.
#'
#' @usage
#' \code{data(fusarium)} \cr
#' \code{Y <- fusarium[["Y"]]}
#'
#' @format A list of one item:
#' \describe{
#'   \item{\code{Y}}{Three-way tensor of disease ratings for 20 varieties of wheat infected with 7 strains of Fusarium head blight over 4 years, 1990-1993.}
#'   ...
#' }
#' @source Fusarium data used and printed in \url{http://www.jstor.org/stable/2533660?seq=1#page_scan_tab_contents}, processed as follows: \cr
#' \cr
#' \code{fus <- read.table("fusarium.txt", header = FALSE) # Text file of data from paper} \cr
#' \code{Y <- array(NA, dim = c(4, 7, 20))} \cr
#' \code{for (i in 1:4) Y[i, , ] <- as.matrix(fus[(i - 1)*7 + 1:7, ])} \cr
#' \code{Y <- log((Y/100)/(1 - (Y/100)))} \cr
#' \code{fusarium <- list("Y" = Y)} \cr
"fusarium"
