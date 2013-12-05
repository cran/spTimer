\name{spT.time}
\alias{spT.time}
\title{Timer series information.}
\description{
This function defines the time series in the spatio-temporal data.
}
\usage{
spT.time(t.series, segments=1)
}
\arguments{
\item{t.series}{Number of times within each segment in each series. Can take both regular and irregular time-series. }
\item{segments}{Number of segments in each time series. This should  be a constant. }
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
## 

# regular time-series in each year
time.data<-spT.time(t.series=30,segments=2)

# irregular time-series in each year
time.data<-spT.time(t.series=c(35,25),segments=2)

##
}
}
\keyword{spT}