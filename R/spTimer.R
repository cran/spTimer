library(coda)

.First.lib <-
    function(libname, pkgname)
{
	## figureout this year automatically
	#this.year <- substr(as.character(Sys.Date()),1,4)

	## output to screen
	#cat("##\n## Spatio-Temporal Bayesian Modelling using R \n")
	#cat("## 2010-", this.year, 
	#	", License: GPL\n##\n", sep="")
	#library.dynam(pkgname, pkgname, lib.loc=libname)
	library.dynam("coda", pkgname, lib.loc=libname)
}

#.onLoad <-
#    function(libname, pkgname)
#{
#	## figureout this year automatically
#	this.year <- substr(as.character(Sys.Date()),1,4)
#
#	## output to screen
#	cat("##\n## Spatio-Temporal Bayesian Modelling using R \n")
#	cat("## 2010-", this.year, 
#		", License: GPL\n##\n", sep="")
#	#cat("##\n## Spatio-Temporal Bayesian Modelling using R \n##\n")
#	library.dynam(pkgname, pkgname, lib.loc=libname)
#}



.onLoad <-
    function(libname, pkgname)
{
	## figureout this year automatically
	this.year <- substr(as.character(Sys.Date()),1,4)

	## output to screen
       packageStartupMessage("##\n## Spatio-Temporal Bayesian Modelling using R")
       packageStartupMessage("## 2010-", this.year, 
		", License: GPL\n##", sep="")
      #packageStartupMessage("## initializing ... spTimer ... ", appendLF = FALSE)
      Sys.sleep(0.2)
      #packageStartupMessage(" done\n##")
      library.dynam(pkgname, pkgname, lib.loc=libname)
}



.onAttach <-
    function(libname, pkgname)
{
	library.dynam(pkgname, pkgname, lib.loc=libname)
}

