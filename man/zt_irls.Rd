% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zt_irls.R
\name{zt_irls}
\alias{zt_irls}
\title{The IRLS Function}
\usage{
zt_irls(formula, data, epsilon = 0.001, maxit = 10, ...)
}
\arguments{
\item{formula}{A formula of the equation that should be used for regression.}

\item{data}{The data that should be used. Formula should use variable names as seen in the data.}

\item{epsilon}{This is the criterion for convergence: set manually, for getting a feel for how it works}

\item{maxit}{The maximum number of iterations you are willing to sit through before you blow your brains out...}

\item{...}{Not currently supported.}
}
\description{
This function is designed to perform an IRLS algorithm to calculate parameters for regression on a given dataset (this will be designed to handle count data).
}
\details{
This is a practice function to see if I'm understanding the concepts (and if it actually runs). Documentation is provided purely for practice.
}
\examples{
# using the "fish" dataset, downloaded with this package
# ztfish <- read.csv("http://www.ats.ucla.edu/stat/data/fish.csv") #if you want to grab it manually
# ztfish <- ztfish[ztfish$count>0,]

#an example with a two way interaction
zt_irls(count~persons*camper, ztfish)

#an example with a three way interaction
zt_irls(count~persons*camper*nofish, ztfish)

}
\seealso{
\code{\link[pscl]{zeroinfl}}
}
\author{
Michael Floren
}
