% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scip_dat.R
\name{scip_dat}
\alias{scip_dat}
\title{Generation of Small Count Inflated Poisson Data}
\usage{
scip_dat(n, desmat, sam, kop, gam, c)
}
\arguments{
\item{n}{The number of data points to be generated.}

\item{desmat}{The design matrix to use for data generation.}

\item{sam}{A vector of parameters to use for the small count distribution.}

\item{kop}{A vector of parameters to use for the large count distribution.}

\item{gam}{A vector of parameters to use for the binomial distribution (1 is small count group).}

\item{c}{The cutoff of inflation}
}
\value{
This function currently returns a list of two vectors
  \item{y}{A vector of the generated count variables}
  \item{g}{A vector of the group each generated count belongs to}
}
\description{
This function generates small count inflated Poisson (SCIP) data, giving a list of values and group membership.
}
\details{
This function takes a set of arguments and returns a vector of two elements. The first element is the generated count from the SCIP distribution, the second element is the group of membership (for dissertation purposes). The returning of the group membership element may be removed at a later date.
}
\examples{
#show a histogram of small count inflated data
hist(scip_dat(n=3000, desmat=matrix(rep(1, 3000), ncol=1), sampi=log(.1), koppa=log(4), gamma=log(1), c=3)$y)

}
\author{
Michael Floren
}
