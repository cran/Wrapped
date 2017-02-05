\name{wrappedg}
\alias{dwrappedg}
\alias{pwrappedg}
\alias{qwrappedg}
\alias{rwrappedg}
\alias{mwrappedg}
\title{Wrapped G Distribution}
\description{Computes the pdf, cdf, quantile, random numbers  for wrapped G distributions. Also computes the Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status when the distribution is fitted  to some data}
\usage{
dwrappedg(x, spec, K = K, ...)
pwrappedg(x, spec, K = K, ...)
qwrappedg(p, spec, K = K, ...)
rwrappedg(n, spec, ...)
mwrappedg(g, data, starts, K = K, method = "BFGS")
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{K}{the limit used to approximate the wrapped G distribution, for example, the wrapped G pdf is approximated by the sum of g(x+2*pi*k, ...) from k=-K to k=K, where g denotes the pdf specified by spec}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{...}{other parameters}
  \item{g}{same as spec but must be one of normal ("norm", package base), Gumbel ("gumbel, package evd), logistic ("logis", package base), Student t ("t", package base), Cauchy ("cauchy", package base), skew normal ("sn", package sn), skew t ("st", package st), asymmetric Laplace ("ald", package ald), normal Laplace ("nl", NormalLaplace), skew Laplace ("skewlap", GeneralizedHyperbolic), generalized logistic ("glogis", package glogis), four parameter logistic ("4pl", package irtProb), skew logistic ("sl", package sld), exponential power ("normalp", package normalp), skew exponential power of type 1 ("SEP1", package gamlss.dist), skew exponential power of type 2 ("SEP2", package gamlss.dist), skew exponential power of type 3 ("SEP3", package gamlss.dist), skew exponential power of type 4 ("SEP4", package gamlss.dist), normal exponential t ("NET", package gamlss.dist),
skew normal type 2 ("SN2", package gamlss.dist), skew generalized t ("sgt", package sgt), skew hyperbolic ("skewhyp", package  SkewHyperbolic), generalized hyperbolic Student t ("ght", package fBasics), power hyperola logistic ("logishp", package FatTailsR), Kiener ("kiener1", package  FatTailsR), Laplace mixture ("laplacem", package LaplacesDemon), skew Laplace ("skewlaplace", package  rmutil), asymmetric Laplace ("asl", package cubfits), asymmetric Laplace ("asla", package cubfits), asymmetric Laplace ("al", package lqmm), polynomial tail Laplace ("ptl", package LCA), generalized asymmetric t ("gat", pacakge GEVStableGarch), variance gamma ("vg", package  VarianceGamma), normal inverse gamma ("nig", package  GeneralizedHyperbolic), skew Cauchy ("sc", package sn), slash ("slash", package VGAM), ex Gaussian ("exGAUS", package  gamlss.dist), stable ("stable", package stable) and log gamma ("lgamma", package  ordinal).  For details  including the density function and parameter specifications see Nadarajah and Zhang (2017)}
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{starts}{initial values for the parameters of the distributions specified by \code{g}}
  \item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated or an object giving the values of Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status.}
\references{S. Nadarajah and Y. Zhang, Wrapped: A new R package for wrapped distributions, submitted}
\author{Saralees Nadarajah, Yuanyuan Zhang}
\examples{x=runif(10,min=-1,max=1)
dwrappedg(x,"norm",mean=1,sd=1,K=1000)
pwrappedg(x,"norm",mean=1,sd=1,K=1000)
qwrappedg(runif(10,min=0,max=1),"norm",mean=1,sd=1,K=1000)
rwrappedg(10,"norm",mean=1,sd=1)
mwrappedg("norm",runif(100,min=-1,max=1),starts=c(1,1),K=10,method="BFGS")}
\keyword{Wrapped G distribution}
