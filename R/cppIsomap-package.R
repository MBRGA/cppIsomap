#' cppIsomap: RCpp port of RDRToolbox's Isomap embedding implementation
#' 
#' Provides cppIsomap function that computes the Isomap embedding as introduced in 2000 by Tenenbaum, de Silva and Langford. 
#' Direct port of RDRToolbox version. Does not implement residual plot.
#' 
#' @details The only function exported by this package is cppIsomap, which has the same parameters as RDRToolbox::Isomap.
#' @author Marc Burgess
#' @keywords internal
#' @docType package
#' @name cppIsomap-package
#' @useDynLib cppIsomap
#' @importFrom Rcpp evalCpp
#' @exportPattern "ˆ[[:alpha:]]+"
#' @title cppIsomap
NULL