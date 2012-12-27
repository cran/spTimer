\name{predict.spT}
\alias{predict.spT}
\title{Spatial and temporal predictions for the spatio-temporal models.}
\description{
This function is used to obtain spatial predictions in the unknown locations and also to get the temporal forecasts using MCMC samples.
}
\usage{
## S3 method for class 'spT'
\method{predict}{spT}(object, newdata, newcoords, foreStep=NULL, type="spatial", 
        nBurn, tol.dist, predAR=NULL, ...)

## see also \code{\link{spT.prediction}} for spatial prediction and \code{\link{spT.forecast}} for temporal prediction or forecast
}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{newdata}{The data set for the covariate values for spatial prediction or temporal forecasts. This data should have the same space-time structure as the original data frame.}
\item{newcoords}{The coordinates for the prediction or forecast sites. The locations are in similar format to coords.}
\item{foreStep}{Number of K-step (time points) ahead forecast, K=1,2, ...; Only applicable if type="temporal".}
\item{type}{If "spatial" the do spatial prediction and if "temporal" the do temporal prediction or forecast. Default value is "spatial".}
\item{nBurn}{Number of burn-in. Initial MCMC samples to discard before making inference.}
\item{tol.dist}{Minimum tolerance distance limit between fitted and predicted locations.}
\item{predAR}{The prediction output, if forecasts are in the prediction locations. Only applicable if type="forecast" and data fitted with the "AR" model.}
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
}
\seealso{
\code{\link{spT.Gibbs}}, \code{\link{spT.prediction}}, \code{\link{spT.forecast}}.
}
\examples{
\dontrun{
##

###########################
## The GP models:
###########################

##
## Spatial prediction/interpolation
##

# Read data
data(DataValPred)

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
data(DataValFore);

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
fore.gp <- predict(post.gp, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gp)
names(fore.gp)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.gp$Mean)) # 

##
## Fit and spatially prediction simultaneously
##

# Read data 
data(DataFit); 
data(DataValPred)

# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.gp.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="GP", coords=coords, 
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.gp.fitpred)
summary(post.gp.fitpred)
coef(post.gp.fitpred)
plot(post.gp.fitpred,residuals=TRUE)
names(post.gp.fitpred)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gp.fitpred$prediction[,1]))  

###########################
## The AR models:
###########################

##
## Spatial prediction/interpolation
##

# Read data
data(DataValPred)

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
data(DataValFore);

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
fore.ar <- predict(post.ar, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.ar)
names(fore.ar)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.ar$Mean)) # 

##
## Fit and spatially prediction simultaneously
##

# Read data 
data(DataFit); 
data(DataValPred)

# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.ar.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="AR", coords=coords, 
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.ar.fitpred)
summary(post.ar.fitpred)
coef(post.ar.fitpred)
names(post.ar.fitpred)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.ar.fitpred$prediction[,1]))  

#################################
## The GPP approximation models:
#################################

##
## Spatial prediction/interpolation
##

# Read data
data(DataValPred)

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
data(DataValFore);

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
fore.gpp <- predict(post.gpp, newdata=DataFitFore, newcoords=fore.coords, 
           type="temporal", foreStep=2)
print(fore.gpp)
names(fore.gpp)

# Forecast validations 
spT.validation(DataFitFore$o8hrmax,c(fore.gpp$Mean)) # 

##
## Fit and spatially prediction simultaneously
##

# Read data 
data(DataFit); 
data(DataValPred)

# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
knots.coords<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.gpp.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="GP", coords=coords, knots.coords=knots.coords,
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.gpp.fitpred)
summary(post.gpp.fitpred)
coef(post.gpp.fitpred)
plot(post.gpp.fitpred, residuals=TRUE)
names(post.gpp.fitpred)

# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gpp.fitpred$prediction[,1]))  

##
}
}
\keyword{spT}