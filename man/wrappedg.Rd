\name{wrappedg}
\alias{dwrappedg}
\alias{pwrappedg}
\alias{qwrappedg}
\alias{rwrappedg}
\alias{mwrappedg}
\alias{plotfour}
\title{Wrapped G Distribution}
\description{Computes the pdf, cdf, quantile, random numbers  for any wrapped G distributions. Computes maximum likelihood estimates of the parameters, standard errors, 95 percent confidence intervals, value of Cramer-von Misses statistic, value of Anderson Darling statistic, value of Kolmogorov Smirnov test statistic and its $p$-value, value of Akaike Information Criterion, value of Consistent Akaike Information Criterion, value of Bayesian Information Criterion, value of Hannan-Quinn information criterion, minimum value of the negative log-likelihood function and convergence status when the wrapped distribution is fitted  to some data}
\usage{
dwrappedg(x, spec, K = K, ...)
pwrappedg(x, spec, K = K, ...)
qwrappedg(p, spec, K = K, ...)
rwrappedg(n, spec, ...)
mwrappedg(g, data, starts, K = K, method = "BFGS")
plotfour(g, K = K, para, plotit)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{number of random numbers to be generated}
  \item{K}{the limit used to approximate the wrapped G distribution, for example, the wrapped G pdf is approximated by the sum of g(x+2*pi*k, ...) from k=-K to k=K, where g denotes the pdf specified by spec}
  \item{spec}{a character string specifying the distribution of G and g (for example, "norm" if G and g correspond to the standard normal).}
  \item{...}{other parameters}
  \item{g}{same as spec but must be one of normal ("norm", package base), Gumbel ("gumbel, package evd), logistic ("logis", package base), Student t ("t.scaled", package metRology), Cauchy ("cauchy", package base), skew normal ("sn", package sn), skew t ("st", package sn), asymmetric Laplace ("ald", package ald), normal Laplace ("nl", NormalLaplace), skew Laplace ("skewlap", GeneralizedHyperbolic), generalized logistic ("glogis", package glogis), skew logistic ("sl", package sld), exponential power ("normalp", package normalp), skew exponential power of type 1 ("SEP1", package gamlss.dist), skew exponential power of type 2 ("SEP2", package gamlss.dist), skew exponential power of type 3 ("SEP3", package gamlss.dist), skew exponential power of type 4 ("SEP4", package gamlss.dist), normal exponential t ("NET", package gamlss.dist), skew normal type 2 ("SN2", package gamlss.dist), skew generalized t ("sgt", package sgt), skew hyperbolic ("skewhyp", package  SkewHyperbolic), asymmetric Laplace ("asl", package cubfits), asymmetric Laplace ("asla", package cubfits), asymmetric Laplace ("al", package lqmm), polynomial tail Laplace ("ptl", package LCA), generalized asymmetric t ("gat", pacakge GEVStableGarch), variance gamma ("vg", package  VarianceGamma), normal inverse gamma ("nig", package  GeneralizedHyperbolic), skew Cauchy ("sc", package sn), slash ("slash", package VGAM), generalized hyperbolic Student t ("ght", package fBasics), ex Gaussian ("exGAUS", package  gamlss.dist) and log gamma ("lgamma", package  ordinal). For details  including the density function and parameter specifications see Nadarajah and Zhang (2017)}
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{starts}{initial values for the parameters of the distributions specified by \code{g}}
  \item{method}{the method for optimizing the log likelihood function.  It can be one of \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.  The default is \code{"BFGS"}.  The details of these methods can be found in the manual pages for  \code{optim}}
  \item{para}{a list with four components, each component is a vector specifying
 the parameter values of the specified wrapped distribution}
  \item{plotit}{a character string saying what is to be plotted.  It should take one of "pdf", "cdf", "quantile" or "random".}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random numbers generated or an object giving the values of Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value, maximum likelihood estimates, Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion, standard errors of the maximum likelihood estimates, minimum value of the negative log-likelihood function and convergence status.  \code{plotfour} draws four plots of the PDF, CDF, quantile function
or the histogram of the radii of 100 random numbers of a specified wrapped distribution.}
\references{S. Nadarajah and Y. Zhang, Wrapped: An R package for wrapped distributions, submitted}
\author{Saralees Nadarajah, Yuanyuan Zhang}
\examples{x=runif(10,min=-1,max=1)
dwrappedg(x,"norm",mean=1,sd=1,K=1000)
pwrappedg(x,"norm",mean=1,sd=1,K=1000)
qwrappedg(runif(10,min=0,max=1),"norm",mean=1,sd=1,K=1000)
rwrappedg(10,"norm",mean=1,sd=1)
mwrappedg("norm",runif(100,min=-1,max=1),starts=c(1,1),K=10,method="BFGS")
plotfour("norm",K=100,para=list(c(0,0.1),c(0,2),c(0,5),c(0,20)),
plotit="pdf")}
\keyword{Wrapped G distribution}
