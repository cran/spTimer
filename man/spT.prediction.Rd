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
\item{predicted.samples}{Predicted MCMC samples.}
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
\code{\link{spT.Gibbs}}, \code{\link{spT.forecast}}.
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

# define prediction coordinates and variables
pred.coords<-as.matrix(val.site[,2:3])

# GP prediction
pred.gp <- spT.prediction(nBurn=100, pred.data=val.dat, 
         pred.coords, posteriors=post.gp, tol.dist=2, 
         Summary=TRUE)

names(pred.gp)

# validation criteria and plots
spT.validation(val.dat$o8hrmax,c(pred.gp$Mean)) # 
spT.pCOVER(val.dat$o8hrmax,c(pred.gp$Up),c(pred.gp$Low)) #
spT.segment.plot(val.dat$o8hrmax,c(pred.gp$Mean),
      c(pred.gp$Up),c(pred.gp$Low))
abline(a=0,b=1)


##
## Fit and Prediction simultaneously
##


pred.coords<-as.matrix(val.site[,2:3])

coords<-as.matrix(fit.site[,2:3])
time.data<-spT.time(t.series=61,segment=1)

priors<-spT.priors(model="GP",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))
initials<-spT.initials(model="GP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.01)
spatial.decay<-spT.decay(type="MH", tuning=0.08)
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.01,0.02),segments=5)

nItr<-5000 

# MCMC via Gibbs
post.gp.fit.pred <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=fit.dat, model="GP", time.data=time.data, 
         coords=coords, pred.coords=pred.coords, pred.data=val.dat,
         priors=priors, initials=initials, 
         nItr=nItr, nBurn=1000, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay,
         annual.aggregation="NONE")

names(post.gp.fit.pred)

post.gp.fit.pred$parameter


###########################
## The AR models:
###########################


##
## Prediction for AR models
##

# define prediction coordinates and variables
pred.coords<-as.matrix(val.site[,2:3])

# AR prediction
pred.ar <- spT.prediction(nBurn=100, pred.data=val.dat, 
         pred.coords=pred.coords, 
         posteriors=post.ar, tol.dist=2, Summary=TRUE)

names(pred.ar)

# validation criteria and plots
spT.validation(val.dat$o8hrmax,c(pred.ar$Mean)) # 
spT.pCOVER(val.dat$o8hrmax,c(pred.ar$Up),c(pred.ar$Low)) # 
spT.segment.plot(val.dat$o8hrmax,c(pred.ar$Mean),
      c(pred.ar$Up),c(pred.ar$Low))
abline(a=0,b=1)


##
## fit and predict combinedly for AR models with text output
##

pred.coords<-as.matrix(val.site[,2:3])

coords<-as.matrix(fit.site[,2:3])
time.data<-spT.time(t.series=61,segment=1)

priors<-spT.priors(model="AR",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))

initials<-spT.initials(model="AR", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.01)
spatial.decay<-spT.decay(type="MH", tuning=0.08)
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.01,0.02),segments=5)

nItr<-500
nBurn<-100

# MCMC via Gibbs
post.ar.fit.pred <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=fit.dat, model="AR", time.data=time.data, 
         coords=coords, pred.coords=pred.coords, pred.data=val.dat,
         priors=priors, initials=initials, 
         nItr=nItr, nBurn=nBurn, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay,
         annual.aggregation="an4th")

names(post.ar.fit.pred)

post.ar.fit.pred$para


###########################
## Models with GPP approximations:
###########################


##
## prediction for the GPP
##

# define prediction coordinates and variables
pred.coords<-as.matrix(val.site[,2:3])

# GPP prediction
pred.gpp <- spT.prediction(nBurn=100, pred.data=val.dat, 
       pred.coords=pred.coords, posteriors=post.gpp, tol.dist=2, 
       Summary=TRUE)
 
names(pred.gpp)


# validation criteria and plots
spT.validation(val.dat$o8hrmax,c(pred.gpp$Mean)) #
spT.pCOVER(val.dat$o8hrmax,c(pred.gpp$Up),c(pred.gpp$Low)) #
spT.segment.plot(val.dat$o8hrmax,c(pred.gpp$Mean),
      c(pred.gpp$Up),c(pred.gpp$Low))#,limit=c(0,100))
abline(a=0,b=1)


##
## fit and predict together for the GPP with text output
##

pred.coords<-as.matrix(val.site[,2:3])

coords<-as.matrix(fit.site[,2:3])
knots.coords<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
time.data<-spT.time(t.series=61,segment=1)

priors<-spT.priors(model="GPP",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))
#priors<-NULL
initials<-spT.initials(model="GPP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)
#initials<-NULL

# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.01)
spatial.decay<-spT.decay(type="MH", tuning=0.0008)
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.001,0.002),segments=5)


nItr<-500 

# MCMC via Gibbs
post.gpp.fit.pred <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=fit.dat, model="GPP", time.data=time.data, 
         coords=coords, knots.coords=knots.coords,
         pred.coords=pred.coords, pred.data=val.dat,
         priors=priors, initials=initials, 
         nItr=nItr, nBurn=100, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay,
         annual.aggregation="NONE")

names(post.gpp.fit.pred)

post.gpp.fit.pred$para


##
}
}
\keyword{spT}
