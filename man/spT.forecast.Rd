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
\item{fore.samples}{Forecast MCMC samples.}
\item{fore.coords}{Forecast coordinates.}
\item{Mean}{Average of the MCMC predictions, if Summary=TRUE.}
\item{Median}{Median of the MCMC predictions, if Summary=TRUE.}
\item{SD}{Standard deviation of the MCMC predictions, if Summary=TRUE.}
\item{Low}{Lower limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{Up}{Upper limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{computation.time}{The computation time.}
\item{model}{The model method used for prediction.}
}
\seealso{
\code{\link{predict.spT}}, \code{\link{spT.Gibbs}}, \code{\link{spT.prediction}}.
}
\examples{
\dontrun{
##

###########################
## The GP models:
###########################

##
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
data(DataValFore);

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.gp <- spT.forecast(nBurn=100, K=2, fore.data=DataValFore, 
           fore.coords=fore.coords, posteriors=post.gp, Summary=TRUE)
print(fore.gp)

# forecast validations for the validation sites
spT.validation(DataValFore$o8hrmax,c(fore.gp$Mean)) # 
spT.pCOVER(DataValFore$o8hrmax,c(fore.gp$Up),c(fore.gp$Low)) # 
spT.segment.plot(DataValFore$o8hrmax,c(fore.gp$Mean),
      c(fore.gp$Up),c(fore.gp$Low))
abline(a=0,b=1)

##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
data(DataFitFore)

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.gp <- spT.forecast(nBurn=100, K=2, fore.data=fit.fore, 
           fore.coords=fore.coords, posteriors=post.gp, Summary=TRUE)
print(fore.gp)

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
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
data(DataValFore);

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.ar <- spT.forecast(nBurn=100, K=2, fore.data=DataValFore, 
           fore.coords, pred.samples.ar=pred.ar,
           posteriors=post.ar, Summary=TRUE)
print(fore.ar)

# forecast validations for the validation sites
spT.validation(DataValFore$o8hrmax,c(fore.ar$Mean)) # 
spT.pCOVER(DataValFore$o8hrmax,c(fore.ar$Up),c(fore.ar$Low)) # 
spT.segment.plot(DataValFore$o8hrmax,c(fore.ar$Mean),
      c(fore.ar$Up),c(fore.ar$Low))
abline(a=0,b=1)

##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
data(DataFitFore)

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.ar <- spT.forecast(nBurn=100, K=2, fore.data=fit.fore, 
           fore.coords, pred.samples.ar=NULL,
           posteriors=post.ar, Summary=TRUE)
print(fore.ar)

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
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
data(DataValFore);

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.gpp <- spT.forecast(nBurn=100, K=2, fore.data=DataValFore, 
           fore.coords=fore.coords, posteriors=post.gpp, 
           Summary=TRUE)

names(fore.gpp)

# forecast validations for the validation sites
spT.validation(DataValFore$o8hrmax,c(fore.gpp$Mean)) #
spT.pCOVER(DataValFore$o8hrmax,c(fore.gpp$Up),c(fore.gpp$Low)) # 
spT.segment.plot(DataValFore$o8hrmax,c(fore.gpp$Mean),
      c(fore.gpp$Up),c(fore.gpp$Low))
abline(a=0,b=1)

##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
data(DataFitFore)

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.gpp <- spT.forecast(nBurn=100, K=2, fore.data=fit.fore, 
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
