\name{spT.Gibbs}
\alias{spT.Gibbs}
\title{MCMC sampling for the spatio-temporal models.}
\description{
This function is used to draw MCMC samples using the Gibbs sampler. 
}
\usage{
spT.Gibbs(formula, data=parent.frame(), model, time.data, 
coords, knots.coords, pred.coords, pred.data, priors, initials, 
nItr, nBurn, report, tol.dist, distance.method, cov.fnc, 
scale.transform, spatial.decay, annual.aggrn)
%X.out=TRUE, Y.out=TRUE)
}
\arguments{
\item{formula}{The symnbolic description of the model equation of the regression part of the space-time model.}
\item{data}{An optional data frame containing the variables in the model. If omitted, the variables are taken from environment(formula), typically the environment from which spT.Gibbs is called. The data should be ordered first by the time and then by the sites specified by the \code{coords} below.}
\item{model}{The spatio-temporal models to be fitted, current input: "GP", "AR", and "GPP".}
\item{time.data}{Defining the segments of the time-series set up using the function \code{\link{spT.time}}.}
\item{coords}{The n by 2 matrix defining the locations (e.g., longitude/easting, latitude/northing) of the fitting sites, where n is the number of fitting sites.}
\item{knots.coords}{The locations of the knots in similar format to coords above, only required if \code{model="GPP"}.}
\item{pred.coords}{The locations of the prediction sites in similar format to coords above, only required if fit and predictions are done simultaneously.}
\item{pred.data}{The covariate values at the prediction sites specified by \code{pred.coords}. This should have same space-time structure as the original data frame.}
\item{priors}{The prior distributions for the parameters.  Default distributions are specified if these are not specified. See details in \code{\link{spT.priors}}.}
\item{initials}{The preferred initial values for the parameters. If omitted, default values are provided automatically. Further  details are provided in  \code{\link{spT.initials}}.}
\item{nItr}{Number of MCMC iterations.}
\item{nBurn}{Number of burn-in samples. This number of samples will be discarded before making any inference.}
\item{report}{Number of reports to display while running the Gibbs sampler. Defaults to number of iterations. }
\item{distance.method}{The preferred method to calculate the distance between any two locations. The available options are "geodetic:km", "geodetic:mile", "euclidean", "maximum", "manhattan", and "canberra". See details in \code{\link{dist}}.}
\item{tol.dist}{Minimum separation distance between any two locations out of those  specified by coords, knots.coords and pred.coords. The default is 0.005. The programme will exit if the minimum distance is less than the non-zero specified value. This will ensure non-singularity of the covariance matrices.}
\item{cov.fnc}{Covariance function for the spatial effects. The available options  are "exponential", "gaussian", "spherical" and "matern".}
\item{scale.transform}{The transformation method for the response variable. Currently implemented options are: "NONE", "SQRT", and "LOG" with  "NONE" as the deault.}
\item{spatial.decay}{Provides the options for sampling the spatial decay parameter \eqn{\phi}. Currently implemented options  are "DISCRETE", "MH" or "FIXED" and further options for each of these are specified by \code{\link{spT.decay}}. The default is "MH".}
\item{annual.aggrn}{This provides the options for calculating annual summary statistics by aggregating different time segments (e.g., annual mean). Currently implemented options are: "NONE", "ave" and "an4th", where "ave" = annual average, "an4th"= annual 4th highest. Only applicable if \code{\link{spT.time}} inputs more than one segment and when fit and predict are done simultaneously.}
%\item{X.out}{Logical value, if TRUE then the covariates are present in the output.}
%\item{Y.out}{Logical value, if TRUE then the dependent variable is present in the output.}
}
\value{
%\item{X}{The covariate values, if X.out=TRUE.}
%\item{Y}{The observed values, if Y.out=TRUE.}
\item{accept}{The acceptance rate for the \eqn{\phi} parameter if the "MH" method of sampling is chosen.}
\item{phip}{MCMC samples for the parameter \eqn{\phi}.}
\item{sig2ep}{MCMC samples for the parameter \eqn{\sigma^2_\epsilon}.}
\item{sig2etap}{MCMC samples for the parameter \eqn{\sigma^2_\eta}.}
\item{betap}{MCMC samples for the parameter \eqn{\beta}.}
\item{op}{MCMC samples for the true observations.}
\item{fitted}{MCMC summary (mean and sd) for the fitted values.}
\item{sig2ep}{MCMC samples for the parameter \eqn{\sigma^2_\epsilon}.}
\item{tol.dist}{Minimum tolerance distance limit between the locations.}
\item{distance.method}{Name of the distance calculation method.}
\item{cov.fnc}{Name of the covariance function used in model fitting.}
\item{scale.transform}{Name of the scale.transformation method.}
\item{sampling.sp.decay}{The method of sampling for the spatial decay parameter \eqn{\phi}.}
\item{covariate.names}{Name of the covariates used in the model.}
\item{Distance.matrix}{The distance matrix.}
\item{coords}{The coordinate values.}
\item{n}{Total number of sites.}
\item{r}{Total number of segments in time, e.g., years.}
\item{T}{Total points of time, e.g., days within each year.}
\item{p}{Total number of model coefficients, i.e., \eqn{\beta}'s including the intercept.}
\item{initials}{The initial values used in the model.}
\item{priors}{The prior distributions used in the model.}
\item{PMCC}{The predictive model choice criteria obtained by minimising the expected value of a loss function, see Gelfand and Ghosh (1998). Results for both goodness of fit and penalty are given.}
\item{iterations}{The number of samples for the MCMC chain, without burn-in.}
\item{nBurn}{The number of burn-in period for the MCMC chain.}
\item{computation.time}{The computation time required for the fitted model.}
\item{model}{The spatio-temporal model used for analyse the data.}

\item{Text Output}{This option is only applicable when fit and predictions are done simultaneously.\cr

For GP models:\cr
OutGP_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, sig2eps, sig2eta, and phi. \cr
OutGP_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations, where N=nrT. \cr
OutGP_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGP_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutGP_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutGP_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr


For AR models: \cr
OutAR_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, rho, sig2eps, sig2eta, mu_l's, sig2l's and phi. \cr
OutAR_Stats_TrueValue.txt: (N x 2) matrix of true summary values, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutAR_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutAR_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr


For models using GPP approximations: \cr
OutGPP_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, rho, sig2eps, sig2eta, and phi. \cr
OutGPP_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGPP_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGPP_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutGPP_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutGPP_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr

}
}
%\note{
%}
\seealso{
\code{\link{spT.priors}, \link{spT.initials}, \link{spT.geodist}, \link{dist}, \link{spT.MCMC.stat}, \link{spT.MCMC.plot}, \link{spT.prediction}, \link{spT.forecast}}.
}
\examples{
\dontrun{
##

###########################
## Read data:
###########################

library(spTimer)

data(NYdata)
data(NYsite)

head(NYdata)
head(NYsite)

# map
library(maps)
map(database="state",regions="new york")
points(NYsite[,2:3],pch=19)

# descriptive stat

summary(NYdata[,7:10])

# 
n.sites <- length(unique(NYsite[,1])) # Total number of sites
n.valsites <- 4                       # Let 4 sites are set aside  
nall.sites <- 1:n.sites               # Numbered sites 
nrand.sites <- sort(sample(nall.sites, n.valsites)) # Select 4 sites randomly
nval.sites <- nrand.sites              # Numbered validation sites 
nfit.sites <- nall.sites[-nrand.sites] # Numbered fitted sites

# split the NYdata into fitted and validation part

val.dat <- spT.data.selection(data=NYdata,random=F,s=nval.sites)
fit.dat <- spT.data.selection(data=NYdata,random=F,s=nfit.sites)

val.site <- spT.data.selection(data=NYsite,random=F,s=nval.sites)
fit.site <- spT.data.selection(data=NYsite,random=F,s=nfit.sites)

# set aside last day of August for forecast validation

fit.fore<-fit.dat[fit.dat$Month==8 & fit.dat$Day==31,]
val.fore<-val.dat[val.dat$Month==8 & val.dat$Day==31,]

# new fitted and validation data with 61 days

fit.dat<-fit.dat[fit.dat$Month !=8 | fit.dat$Day !=31,]
val.dat<-val.dat[val.dat$Month !=8 | val.dat$Day !=31,]


#
library(maps)
map(database="state",regions="new york")
points(fit.site[,2:3],pch=2,col=4)
points(val.site[,2:3],pch=19,col=2)
#

###########################
## The GP models:
###########################

# define the coordinates
coords<-as.matrix(fit.site[,2:3])

# define the time-series 
time.data<-spT.time(t.series=61,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="GP",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))
#priors<-NULL

# initial values for the model parameters
initials<-spT.initials(model="GP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)
#initials<-NULL

# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.01)
spatial.decay<-spT.decay(type="MH", tuning=0.08)
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.01,0.02),segments=5)


# input for the MCMC algorithms
nItr<-5000

# MCMC via Gibbs
post.gp <- spT.Gibbs(formula=o8hrmax ~ cMAXTMP+WDSP+RH, 
         data=fit.dat, model="GP", time.data=time.data, 
         coords=coords, priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)

#object.size(post.gp)/(1024*1024) #MB
names(post.gp)

# model selection criteria
post.gp$PMCC # 

# MCMC summary and plots
spT.MCMC.stat(post.gp,nBurn=1000)
spT.MCMC.plot(post.gp,nBurn=1000)

# Use of coda pakcage

library(coda)
tmp<-as.mcmc(post.gp)
plot(tmp)
summary(tmp)


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

# define the coordinates
coords<-as.matrix(fit.site[,2:3])

# define the time-series 
time.data<-spT.time(t.series=61,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="AR",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))

# initial values for the model parameters
initials<-spT.initials(model="AR", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.01)
spatial.decay<-spT.decay(type="MH", tuning=0.08)
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.01,0.02),segments=5)

# input for the MCMC algorithms
nItr<-500


# MCMC via Gibbs
post.ar <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=fit.dat, model="AR", time.data=time.data, 
         coords=coords, priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)

object.size(post.ar)/(1024*1024) #MB
names(post.ar)

# model selection criteria
post.ar$PMCC # 

# MCMC summary and plots
spT.MCMC.stat(post.ar,nBurn=100)
spT.MCMC.plot(post.ar,nBurn=100)

# Use of coda pakcage
library(coda)
tmp<-as.mcmc(post.ar)
plot(tmp)
summary(tmp)


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


# define the coordinates and knots
coords<-as.matrix(fit.site[,2:3])
knots.coords<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

library(maps)
map(database="state",regions="new york")
points(coords,pch=2,col=2)
points(knots.coords,pch=19,col=4)

# define the time-series 
time.data<-spT.time(t.series=61,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="GPP",var.prior=Gam(2,1),
        beta.prior=Nor(0,10^4))
#priors<-NULL

# initial values for the model parameters
initials<-spT.initials(model="GPP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)
#initials<-NULL


# input for spatial decay
#spatial.decay<-spT.decay(type="FIXED", value=0.001)
spatial.decay<-spT.decay(type="MH", tuning=0.05) # 
#spatial.decay<-spT.decay(type="DISCRETE",limit=c(0.001,0.009),segments=10)


# input for the MCMC algorithms
nItr<-5000 

# MCMC via Gibbs
post.gpp <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=fit.dat, model="GPP", time.data=time.data, 
         coords=coords, knots.coords=knots.coords,
         priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr/10, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)


object.size(post.gpp)/(1024*1024) #MB
names(post.gpp)


# model selection criteria
post.gpp$PMCC # 

# MCMC summary and plots
spT.MCMC.stat(post.gpp,nBurn=1000)
spT.MCMC.plot(post.gpp,nBurn=1000,ACF=T)

# Use of coda pakcage
library(coda)
tmp<-as.mcmc(post.gpp)
plot(tmp)
summary(tmp)


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
