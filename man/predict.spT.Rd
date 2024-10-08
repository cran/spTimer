\name{predict.spT}
\alias{predict.spT}
\title{Spatial and temporal predictions for the spatio-temporal models.}
\description{
This function is used to obtain spatial predictions in the unknown locations and also to get the temporal forecasts using MCMC samples.
}
\usage{
\method{predict}{spT}(object, newdata, newcoords, foreStep=NULL, type="spatial", 
        nBurn, tol.dist, predAR=NULL, Summary=TRUE, ...)

}

\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{newdata}{The data set providing the covariate values for spatial prediction or temporal forecasts. This data should have the same space-time structure as the original data frame.}
\item{newcoords}{The coordinates for the prediction or forecast sites. The locations are in similar format to \code{coords}, see \code{\link{spT.Gibbs}}.}
\item{foreStep}{Number of K-step (time points) ahead forecast, K=1,2, ...; Only applicable if type="temporal".}
\item{type}{If the value is "spatial" then only spatial prediction will be performed at the \code{newcoords} which must be different from the fitted sites provided by \code{coords}.  When the "temporal" option is specified then forecasting will be performed  and in this case the \code{newcoords} may also contain elements of the fitted sites in which case only temporal forecasting beyond the last fitted time point will be performed. } 
\item{nBurn}{Number of burn-in. Initial MCMC samples to discard before making inference.}
\item{tol.dist}{Minimum tolerance distance limit between fitted and predicted locations.}
\item{predAR}{The prediction output, if forecasts are in the prediction locations. Only applicable if type="forecast" and data fitted with the "AR" model.}
\item{Summary}{To obtain summary statistics for the posterior predicted MCMC samples. Default is TRUE.}
\item{...}{Other arguments.}
}
\value{
\item{pred.samples or fore.samples}{Prediction or forecast MCMC samples.}
\item{pred.coords or fore.coords}{prediction or forecast coordinates.}
\item{Mean}{Average of the MCMC predictions}
\item{Median}{Median of the MCMC predictions}
\item{SD}{Standard deviation of the MCMC predictions}
\item{Low}{Lower limit for the 95 percent CI of the MCMC predictions}
\item{Up}{Upper limit for the 95 percent CI of the MCMC predictions}
\item{computation.time}{The computation time.}
\item{model}{The model method used for prediction.}
\item{type}{"spatial" or "temporal".}
\item{...}{Other values "obsData", "fittedData" and "residuals" are provided only for temporal prediction. This is to analyse the \code{spTimer} forecast output using package \code{forecast} through function \code{\link{as.forecast.object}}.}
}
\references{

 Bakar, K. S. and Sahu, S. K. (2014) spTimer: Spatio-Temporal Bayesian Modelling Using R. Technical Report, University of Southampton, UK. To appear in the Journal of Statistical Software. \cr
 Sahu, S. K. and Bakar, K. S. (2012) A comparison of Bayesian Models for Daily Ozone Concentration Levels Statistical Methodology , 9, 144-157. \cr
 Sahu, S. K. and Bakar, K. S. (2012) Hierarchical Bayesian auto-regressive models for large space time data with applications to ozone concentration modelling. Applied Stochastic Models in Business and Industry, 28, 395-415. \cr
}
\seealso{
\code{\link{spT.Gibbs}, \link{as.forecast.object}}.
}
\examples{
\dontshow{
# Read data
data(NYdata)
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))
#
nItr <- 100
nBurn <- 0
# GP: MCMC via Gibbs using default choices
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        nItr=nItr,  nBurn=nBurn,
        data=DataFit, model="GP", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
pred.gp <- predict(post.gp, newdata=DataValPred, newcoords=~Longitude+Latitude)
fore.gp <- predict(post.gp, newdata=DataValFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2)
# AR: MCMC via Gibbs using default choices
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        nItr=nItr,  nBurn=nBurn,
        data=DataFit, model="AR", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
pred.ar <- predict(post.ar, newdata=DataValPred, newcoords=~Longitude+Latitude)
fore.ar <- predict(post.ar, newdata=DataValFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2)
# GPP: MCMC via Gibbs using default choices
# Define knots
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
set.seed(11)
post.gpp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
         nItr=nItr,   nBurn=nBurn, 
         data=DataFit, model="GPP", coords=~Longitude+Latitude, 
         knots.coords=knots, scale.transform="SQRT")
pred.gpp <- predict(post.gpp, newdata=DataValPred, newcoords=~Longitude+Latitude)
fore.gpp <- predict(post.gpp, newdata=DataValFore, newcoords=~Longitude+Latitude, 
           type="temporal", foreStep=2)
# GP (truncated): MCMC via Gibbs 
data(NYdata)
# Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))
out <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH,data=DataFit,
  model="truncatedGP",coords=~Longitude+Latitude,
  distance.method="geodetic:km",nItr=nItr, nBurn=nBurn,
  truncation.para = list(at = 30,lambda = 4))
pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
# GPP (truncated):MCMC via Gibbs 
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="truncatedGPP",coords=~Longitude+Latitude,
         knots.coords=knots, distance.method="geodetic:km",
         nItr=nItr, nBurn=nBurn,
         truncation.para = list(at = 30,lambda = 4))
pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
##	   
}
\donttest{
##

###########################
## The GP models:
###########################

##
## Spatial prediction/interpolation
##

# Read data
data(NYdata)
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))

# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        data=DataFit, model="GP", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
print(post.gp)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# Spatial prediction using spT.Gibbs output
set.seed(11)
pred.gp <- predict(post.gp, newdata=DataValPred, newcoords=pred.coords)
print(pred.gp)
names(pred.gp)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(pred.gp$Mean))  

##
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.gp <- predict(post.gp, newdata=DataValFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gp)
names(fore.gp)

# Forecast validations 
spT.validation(DataValFore$o8hrmax,c(fore.gp$Mean)) 

# Use of "forecast" class
#library(forecast)
#tmp<-as.forecast.object(fore.gp, site=1) # default for site 1
#plot(tmp)
#summary(tmp)

##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.gp <- predict(post.gp, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gp)
names(fore.gp)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.gp$Mean)) # 

# Use of "forecast" class
#library(forecast)
#tmp<-as.forecast.object(fore.gp, site=5) # for site 5
#plot(tmp)

###########################
## The AR models:
###########################

##
## Spatial prediction/interpolation
##

# Read data
data(NYdata)
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))

# MCMC via Gibbs using default choices
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="AR", coords=~Longitude+Latitude, 
         scale.transform="SQRT")
print(post.ar)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# Spatial prediction using spT.Gibbs output
set.seed(11)
pred.ar <- predict(post.ar, newdata=DataValPred, newcoords=pred.coords)
print(pred.ar)
names(pred.ar)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(pred.ar$Mean))  

##
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.ar <- predict(post.ar, newdata=DataValFore, newcoords=fore.coords, 
           type="temporal", foreStep=2, predAR=pred.ar)
print(fore.ar)
names(fore.ar)

# Forecast validations 
spT.validation(DataValFore$o8hrmax,c(fore.ar$Mean)) 

# Use of "forecast" class
#tmp<-as.forecast.object(fore.ar, site=1) # default for site 1
#plot(tmp)


##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.ar <- predict(post.ar, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.ar)
names(fore.ar)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.ar$Mean)) # 

# Use of "forecast" class
#tmp<-as.forecast.object(fore.ar, site=1) # default for site 1
#plot(tmp)

#################################
## The GPP approximation models:
#################################

##
## Spatial prediction/interpolation
##

# Read data
data(NYdata)
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))

# Define knots
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

# MCMC via Gibbs using default choices
set.seed(11)
post.gpp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="GPP", coords=~Longitude+Latitude, 
         knots.coords=knots, scale.transform="SQRT")
print(post.gpp)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# Spatial prediction using spT.Gibbs output
set.seed(11)
pred.gpp <- predict(post.gpp, newdata=DataValPred, newcoords=pred.coords)
print(pred.gpp)
names(pred.gpp)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(pred.gpp$Mean))  

##
## Temporal  prediction/forecast 
## 1. In the unobserved locations
##

# Read data
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))

# define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataValFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62 
# in the unobserved locations using output from spT.Gibbs
set.seed(11)
fore.gpp <- predict(post.gpp, newdata=DataValFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gpp)
names(fore.gpp)

# Forecast validations 
spT.validation(DataValFore$o8hrmax,c(fore.gpp$Mean)) 

# Use of "forecast" class
#tmp<-as.forecast.object(fore.gpp, site=1) # default for site 1
#plot(tmp)

##
## Temporal  prediction/forecast 
## 2. In the observed/fitted locations
##

# Read data
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

# Define forecast coordinates
fore.coords<-as.matrix(unique(cbind(DataFitFore[,2:3])))

# Two-step ahead forecast, i.e., in day 61 and 62, 
# in the fitted locations using output from spT.Gibbs
set.seed(11)
fore.gpp <- predict(post.gpp, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gpp)
names(fore.gpp)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.gpp$Mean)) # 

# Use of "forecast" class
#tmp<-as.forecast.object(fore.gpp, site=1) # default for site 1
#plot(tmp)

##

######################################################
## The Truncated/Censored GP models:
######################################################

##
## Model fitting
##

data(NYdata)

# Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

#
nItr <- 5000 # number of MCMC samples for each model
nBurn <- 1000 # number of burn-in from the MCMC samples
# Truncation at 30 
# fit truncated GP model 
out <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH,data=DataFit,
  model="truncatedGP",coords=~Longitude+Latitude,
  distance.method="geodetic:km",nItr=nItr,nBurn=nBurn,report=5,
  truncation.para = list(at = 30,lambda = 4),
  fitted.values="ORIGINAL")
#  
summary(out)
head(fitted(out))
plot(out,density=FALSE)
#
head(cbind(DataFit$o8hrmax,fitted(out)[,1]))
plot(DataFit$o8hrmax,fitted(out)[,1])
spT.validation(DataFit$o8hrmax,fitted(out)[,1])

##
## prediction (spatial)
##

pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
names(pred)
plot(DataValPred$o8hrmax,c(pred$Mean)) 
spT.validation(DataValPred$o8hrmax,c(pred$Mean)) 
#pred$prob.below.threshold

##
## forecast (temporal)
##

# unobserved locations
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataValFore$o8hrmax,c(fore$Mean)) 
plot(DataValFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold

# observed locations 
fore <- predict(out,newdata=DataFitFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataFitFore$o8hrmax,c(fore$Mean)) 
plot(DataFitFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold


######################################################
## The Truncated/Censored GPP models:
######################################################

##
## Model fitting
##

data(NYdata)

# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
# Define knots
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

# Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

#
nItr <- 5000 # number of MCMC samples for each model
nBurn <- 1000 # number of burn-in from the MCMC samples
# Truncation at 30 
# fit truncated GPP model 
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="truncatedGPP",coords=~Longitude+Latitude,
         knots.coords=knots, distance.method="geodetic:km",
         nItr=nItr,nBurn=nBurn,report=5,fitted="ORIGINAL",
         truncation.para = list(at = 30,lambda = 4))
#  
summary(out)
head(fitted(out))
plot(out,density=FALSE)
#
head(cbind(DataFit$o8hrmax,fitted(out)[,1]))
plot(DataFit$o8hrmax,fitted(out)[,1])
spT.validation(DataFit$o8hrmax,fitted(out)[,1])

##
## prediction (spatial)
##

pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
names(pred)
plot(DataValPred$o8hrmax,c(pred$Mean)) 
spT.validation(DataValPred$o8hrmax,c(pred$Mean)) 
#pred$prob.below.threshold

##
## forecast (temporal)
##

# unobserved locations
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataValFore$o8hrmax,c(fore$Mean)) 
plot(DataValFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold

# observed locations 
fore <- predict(out,newdata=DataFitFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataFitFore$o8hrmax,c(fore$Mean)) 
plot(DataFitFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold


######################################################
######################################################
##
}
}
\keyword{spT}