\name{spT.MCMC.stat}
\alias{spT.MCMC.stat}
\title{MCMC summary for the posterior samples}
\description{
This function is used to obtain the MCMC summary of the posterior samples.  
}
\usage{
spT.MCMC.stat(x, nBurn)
}
\arguments{
\item{x}{The posterior MCMC samples obtained from \code{\link{spT.Gibbs}}.}
\item{nBurn}{The number of burn-in.}
}
\seealso{
\code{\link{spT.Gibbs}, \link{spT.MCMC.plot}}.
}
\examples{
\dontrun{
##

# for GP models
spT.MCMC.stat(post.gp,nBurn=1000)

# for AR models
spT.MCMC.stat(post.ar,nBurn=1000)

# for GPP based models
spT.MCMC.stat(post.gpp,nBurn=1000)

##
}
}
\keyword{utility}