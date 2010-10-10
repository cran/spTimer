\name{spT.forecast}
\alias{spT.forecast}
\title{Forecast for the spatio-temporal models.}
\description{
This function is used to get the K-step forecast values using MCMC samples.
}
\usage{
spT.forecast(nBurn, K=1, fore.data, fore.coords, 
     posteriors, pred.samples.ar=NULL, Summary=TRUE)
}
\arguments{
\item{nBurn}{Number of burn-in. Initial MCMC samples to discard before making inference.}
\item{K}{Number of K-step (time points) ahead forecast.}
\item{fore.data}{The data set for the covariate values for forecasts. This data should have the same space-time structure as the original data frame.}
\item{fore.coords}{The coordinates for the forecast sites. The locations are in similar format to coords.}
\item{posteriors}{Posterior MCMC samples obtained from function \code{\link{spT.Gibbs}}.}
\item{pred.samples.ar}{The prediction output of the function \code{\link{spT.prediction}}, if forecasts are in the prediction locations. Only used for the model: "AR".}
\item{Summary}{Logical, if TRUE then provide MCMC summary statistics for forecast.}
}
\value{
\item{forecast.samples}{Forecast MCMC samples.}
\item{n.fore.sites}{Total number of sites for forecasts.}
\item{Mean}{Average of the MCMC predictions, if Summary=TRUE.}
\item{Median}{Median of the MCMC predictions, if Summary=TRUE.}
\item{SD}{Standard deviation of the MCMC predictions, if Summary=TRUE.}
\item{Low}{Lower limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{Up}{Upper limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{computation.time}{The computation time.}
\item{model}{The model method used for prediction.}
}
\seealso{
\code{\link{spT.Gibbs}}, \code{\link{spT.prediction}}.
}
\examples{
\dontrun{
##

###########################
## The GP models:
###########################

##
##  Foreecasts: GP models
##

# define forecast coordinates and variables
# forecast in the validation sites (hold-out sites)

fore.coords<-as.matrix(val.site[,2:3])

# one-step ahead forecast, i.e., in day 62
fore.gp <- spT.forecast(nBurn=100, K=1, fore.data=val.fore, 
           fore.coords=fore.coords, posteriors=post.gp, Summary=TRUE)

names(fore.gp)

# forecast validations for the validation sites
spT.validation(val.fore$o8hrmax,c(fore.gp$Mean)) # 
spT.pCOVER(val.fore$o8hrmax,c(fore.gp$Up),c(fore.gp$Low)) # 
spT.segment.plot(val.fore$o8hrmax,c(fore.gp$Mean),
      c(fore.gp$Up),c(fore.gp$Low))
abline(a=0,b=1)


# forecast in the fitted sites 
fore.coords<-as.matrix(fit.site[,2:3])

fore.gp <- spT.forecast(nBurn=100, K=1, fore.data=fit.fore, 
           fore.coords=fore.coords, posteriors=post.gp, Summary=TRUE)
names(fore.gp)

# forecast validations for the fitted sites
spT.validation(fit.fore$o8hrmax,c(fore.gp$Mean)) # 
spT.pCOVER(fit.fore$o8hrmax,c(fore.gp$Up),c(fore.gp$Low)) # 
spT.segment.plot(fit.fore$o8hrmax,c(fore.gp$Mean),
      c(fore.gp$Up),c(fore.gp$Low))
abline(a=0,b=1)


###########################
## The AR models:
###########################


##
## Forecasts: AR models
##

# define forecast coordinates and variables
# forecast in the validation sites (hold-out sites)

fore.coords<-as.matrix(val.site[,2:3])

# one-step ahead forecast, i.e., in day 62
fore.ar <- spT.forecast(nBurn=100, K=1, fore.data=val.fore, 
           fore.coords, pred.samples.ar=pred.ar,
           posteriors=post.ar, Summary=TRUE)

names(fore.ar)

# forecast validations for the validation sites
spT.validation(val.fore$o8hrmax,c(fore.ar$Mean)) # 
spT.pCOVER(val.fore$o8hrmax,c(fore.ar$Up),c(fore.ar$Low)) # 
spT.segment.plot(val.fore$o8hrmax,c(fore.ar$Mean),
      c(fore.ar$Up),c(fore.ar$Low))
abline(a=0,b=1)


# forecast in the fitted sites 
fore.coords<-as.matrix(fit.site[,2:3])

fore.ar <- spT.forecast(nBurn=100, K=1, fore.data=fit.fore, 
           fore.coords, pred.samples.ar=NULL,
           posteriors=post.ar, Summary=TRUE)
names(fore.ar)


# forecast validations for the fitted sites
spT.validation(fit.fore$o8hrmax,c(fore.ar$Mean)) # 
spT.pCOVER(fit.fore$o8hrmax,c(fore.ar$Up),c(fore.ar$Low)) # 
spT.segment.plot(fit.fore$o8hrmax,c(fore.ar$Mean),
      c(fore.ar$Up),c(fore.ar$Low))
abline(a=0,b=1)


###########################
## Models with GPP approximations:
###########################


##
## Forecasts: GPP based models
##

# define forecast coordinates and variables
# forecast in the validation sites (hold-out sites)
fore.coords<-as.matrix(val.site[,2:3])

# one-step ahead forecast, i.e., in day 62
fore.gpp <- spT.forecast(nBurn=100, K=1, fore.data=val.fore, 
           fore.coords=fore.coords, posteriors=post.gpp, 
           Summary=TRUE)

names(fore.gpp)

# forecast validations for the validation sites
spT.validation(val.fore$o8hrmax,c(fore.gpp$Mean)) #
spT.pCOVER(val.fore$o8hrmax,c(fore.gpp$Up),c(fore.gpp$Low)) # 
spT.segment.plot(val.fore$o8hrmax,c(fore.gpp$Mean),
      c(fore.gpp$Up),c(fore.gpp$Low))
abline(a=0,b=1)


# forecast in the fitted sites 
forecast.coords<-as.matrix(fit.site[,2:3])

fore.gpp <- spT.forecast(nBurn=100, K=1, fore.data=fit.fore, 
           fore.coords=forecast.coords, posteriors=post.gpp, 
           Summary=TRUE)
names(fore.gpp)

# forecast validations for the fitted sites
spT.validation(fit.fore$o8hrmax,c(fore.gpp$Mean)) # 
spT.pCOVER(fit.fore$o8hrmax,c(fore.gpp$Up),c(fore.gpp$Low)) # 
spT.segment.plot(fit.fore$o8hrmax,c(fore.gpp$Mean),
      c(fore.gpp$Up),c(fore.gpp$Low))
abline(a=0,b=1)


##
}
}
\keyword{spT}
