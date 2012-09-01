\name{spT.prediction}
\alias{spT.prediction}
\title{Predictions for the spatio-temporal models.}
\description{
This function is used to obtain predicted values using the MCMC samples.
}
\usage{
spT.prediction(nBurn, pred.data, pred.coords, posteriors,
              tol.dist=2, Summary=TRUE)
}
\arguments{
\item{nBurn}{Initial MCMC samples to discard before making inference.}
\item{pred.data}{The data set for the covariate values for prediction. This data should have the same space-time structure as the original data frame.}
\item{pred.coords}{The coordinates for the prediction sites. The locations are in similar format to coords.}
\item{posteriors}{Posterior MCMC samples obtained from function \code{\link{spT.Gibbs}}.}
\item{tol.dist}{Minimum tolerance distance limit between fitted and predicted locations.}
\item{Summary}{Logical, if \code{TRUE} then provide MCMC summary statistics for predictions.}
%\item{fitted.X}{Optional: To provide covariate values for the fitted sites.}
}
\value{
\item{pred.samples}{Predicted MCMC samples.}
\item{pred.coords}{Prediction coordinates.}
\item{distance.method}{Name of the distance calculation method.}.
\item{Distance.matrix.pred}{The distance matrix for the prediction locations.}
\item{cov.fnc}{Name of the covariance function used in models.}
\item{predN}{Total number of data points in the prediction.}
\item{Mean}{Average of the MCMC predictions, if Summary=TRUE.}
\item{Median}{Median of the MCMC predictions, if Summary=TRUE.}
\item{SD}{Standard deviation of the MCMC predictions, if Summary=TRUE.}
\item{Low}{Lower limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{Up}{Upper limit for the 95 percent CI of the MCMC predictions, if Summary=TRUE.}
\item{computation.time}{The computation time.}
\item{model}{The model method used for prediction.}
}
\seealso{
\code{\link{predict.spT}}, \code{\link{spT.Gibbs}}, \code{\link{spT.forecast}}.
}
\examples{
\dontrun{
##

###########################
## The GP models:
###########################

##
## Prediction for GP models
##

# Read data
data(DataValPred)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# GP prediction
set.seed(11)
pred.gp <- spT.prediction(nBurn=2000, pred.data=DataValPred, 
         pred.coords, posteriors=post.gp, tol.dist=2, 
         Summary=TRUE)
print(pred.gp)
names(pred.gp)

# validation criteria and plots
spT.validation(DataValPred$o8hrmax,c(pred.gp$Mean)) # 
spT.pCOVER(DataValPred$o8hrmax,c(pred.gp$Up),c(pred.gp$Low)) #
spT.segment.plot(DataValPred$o8hrmax,c(pred.gp$Mean),
      c(pred.gp$Up),c(pred.gp$Low))
abline(a=0,b=1)

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
## Prediction for AR models
##

# Read data
data(DataValPred)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# AR prediction
set.seed(11)
pred.ar <- spT.prediction(nBurn=100, pred.data=DataValPred, 
         pred.coords=pred.coords, 
         posteriors=post.ar, tol.dist=2, Summary=TRUE)
print(pred.ar)

# validation criteria and plots
spT.validation(DataValPred$o8hrmax,c(pred.ar$Mean)) # 
spT.pCOVER(DataValPred$o8hrmax,c(pred.ar$Up),c(pred.ar$Low)) # 
spT.segment.plot(DataValPred$o8hrmax,c(pred.ar$Mean),
      c(pred.ar$Up),c(pred.ar$Low))
abline(a=0,b=1)

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

###########################
## The GPP approximations models:
###########################

##
## prediction for the GPP
##

# Read data
data(DataValPred)

# Define prediction coordinates
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))

# GPP prediction
set.seed(11)
pred.gpp <- spT.prediction(nBurn=100, pred.data=DataValPred, 
       pred.coords=pred.coords, posteriors=post.gpp, tol.dist=2, 
       Summary=TRUE)
 
names(pred.gpp)

# validation criteria and plots
spT.validation(DataValPred$o8hrmax,c(pred.gpp$Mean)) #
spT.pCOVER(DataValPred$o8hrmax,c(pred.gpp$Up),c(pred.gpp$Low)) #
spT.segment.plot(DataValPred$o8hrmax,c(pred.gpp$Mean),
      c(pred.gpp$Up),c(pred.gpp$Low))#,limit=c(0,100))
abline(a=0,b=1)

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
