\name{spT.MCMC.plot}
\alias{spT.MCMC.plot}
\title{MCMC plots for the posterior samples}
\description{
This function is used to obtain MCMC trace plots, density plots, auto-correlation plots, partial auto-correlation plots of the posterior samples. 
}
\usage{
spT.MCMC.plot(x, nBurn=0, ACF=FALSE, PARTIAL.acf=FALSE)
}
\arguments{
\item{x}{The posterior output obtained from \code{\link{spT.Gibbs}} function.}
\item{nBurn}{Number of burn-in required.}
\item{ACF}{Logical value for the auto-correlation function (ACF) plot. Takes TRUE or FALSE.}
\item{PARTIAL.acf}{Logical value for the partial auto-correlation function (PACF) plot. Takes TRUE or FALSE.}
%\item{package}{Input: coda or spTimer.}
%\item{thin}{Thinning interval, see coda package.}
}
\note{
This function can use the package coda to produce plots. 
}
\seealso{
\code{\link{spT.Gibbs}, \link{spT.MCMC.stat}}.
}
\examples{
\dontrun{
##

# for GP models
spT.MCMC.plot(post.gp,nBurn=1000)

# for AR models
spT.MCMC.plot(post.ar,nBurn=1000)

# for GPP based models
spT.MCMC.plot(post.gpp,nBurn=1000)

##
}
}
\keyword{utility}
