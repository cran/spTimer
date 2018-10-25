\name{spT.validation2}
\alias{spT.validation2}
\title{Validation Commands}
\description{
The following function is used to validate the predicted observations with the actual values based on some threshold.
}
\usage{
spT.validation2(z,zhat,cutoff,names=FALSE)
}
\arguments{
\item{z}{The original values (matrix or vector).}
\item{zhat}{The predicted values (matrix or vector).}
\item{cutoff}{The threshold value or cut-off point.}
\item{names}{Logical, if TRUE then print the names of the validation statistics.}
}
\value{
\item{TPR}{True Positive Rate, Sensitivity, Hit rate, Recall}
\item{FPR}{False Positive Rate, False alarm}
\item{FNR}{False Negative Rate, Miss rate}
\item{TNR}{True Negative Rate, Specificity}
\item{Prevalence}{Prevalence}
\item{Accuracy}{Accuracy}
\item{Precision}{Precision, Positive Predictive Value}
\item{FOR}{False Ommission Rate}
\item{LRp}{Positive Likelihood Ratio}
\item{LRn}{Negative Likelihood Ratio}
\item{FDR}{False Discovery Rate}
\item{NPV}{Negative Predictive Value}
\item{DOR}{Diagnostic Odds Ratio}
\item{F1score}{F1 score} 
\item{Heidke.Skill}{Heidke Skill}
}
\seealso{
\code{\link{spT.pCOVER},\link{spT.validation}}.
}
\examples{
%\dontrun{
##

# Create `x', which is the true values.
# Create `y', which is the predicted values.

x <- rnorm(100,0,0.1)
y <- rnorm(100,0,1)
spT.validation2(x, y, cutoff=0,names=TRUE)

##
%}
}
\keyword{utility}