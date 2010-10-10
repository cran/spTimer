\name{spT.Summary.Stat}
\alias{spT.Summary.Stat}
\title{Summary Statistics}
\description{
This function is used to find the summary statistics for the predictions. Here, we get `Mean', `Median, `SD', and Lower and Upper limit of the 95\% interval.
}
\usage{
spT.Summary.Stat(y)
}
\arguments{
\item{y}{The data matrix.}
}
\examples{
\dontrun{
##

# Create `y', which is the prediction values of the observations.

y <- matrix(rnorm(25*2000,5,2),25,2000)
spT.Summary.Stat(y)

## 
}
}
\keyword{utility}