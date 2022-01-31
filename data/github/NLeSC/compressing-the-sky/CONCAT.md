\name{approxdata}
\alias{approxdata}
\title{approximate data based on the models as derived with function getmodel}
\usage{
approxdata(models=c(),do.log=TRUE, dependent="f_int",independent="freq_eff",
group="band",id="runcat",modeltime=TRUE)
}
\arguments{
\item{models}{dataframe of models as produced by function getmodel}
\item{do.log}{Turn independent and devependent variables into log}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
}
\description{
This function takes the models dataframe as produced by getmodel
and uses it to replicate the original data based on the model, which will
be an approximation.
}
\value{
A dataframe with the same size as the dataframe from which models was derived
}
\examples{
\dontrun{
rep.data = approxdata(models=models)
}
}
\name{testcompr}
\alias{testcompr}
\title{test data compression}
\usage{
testcompr(data=c(),path=c())
}
\arguments{
\item{data}{vector of data}
\item{path}{directory that is used for temporarily creating files}
}
\description{
This function takes a vector and then tests the compression of the
vector using nine different compression techniques:
 gzip low (level 1), high (level 9) or default (level NA)
 bzip2 low (level 1), high (level 9) or default (level NA)
 xz low (level 1), high (level 9) or default (level NA)
}
\examples{
\dontrun{
result = testcompr(data=D,path=filepath)
}
}
\name{studyshell}
\alias{studyshell}
\title{shell function for studying the data}
\usage{
studyshell(data=c(), resultpath= c(),CX2=0.95,minN = 30,
                      dependent="f_int",independent="freq_eff",group="band",id="runcat",
                      timecol="taustart_ts",mcmc=FALSE,include_ext1=FALSE,do.res=FALSE,
                      do.log=TRUE,imidcol="imageid",xtrsrc="xtrsrc",
                      error.low=c(),error.up=c(),prec=c(),
                      modeltime=TRUE,do.bar=TRUE,
                      deplabel="Intensity (Jy)")
}
\arguments{
\item{data}{dataframe including the variables (columns): source identifiers (runcat, xtrsrc),
independent variable (freq_eff), timestamps (taustart_ts), data type (extract_type),
independent variables (f_int), 1-sigma error in the independent variable (f_int_err),
grouping variable (band), image id (id)
}
\item{resultpath}{directory where code stores the results}
\item{CX2}{Confidence interval for the chi-square test}
\item{minN}{Minimum number of measurements per frequency band}
\item{dependent}{Name of the dependent variable in the model (f_int in our case)}
\item{independent}{Name of the independent variable in the model (freq_eff in our case)}
\item{group}{Grouping variable (band in our case)}
\item{id}{Identifier of the image, this is mainly used as a key to aid tracing back the
origin of the data}
\item{timecol}{Name of variable to indicate time}
\item{mcmc}{If TRUE then use Markov Chain Monte Carlo regression instead
of ordinary least squares}
\item{include_ext1}{If true then keep the datapoints in radio astronomy for which
the extract type equals 1. Extract type equals 1 data points are interpolations and
not real data. A model is likely to be more accurate if these datapoints are ommited}
\item{do.res}{if true then also output the function residuals as a seperate dataframe}
\item{do.log}{whether to take the log from the indpenedent and dependent variables}
\item{imidcol}{Name of column in which image id is stored}
\item{xtrsrc}{Name of column in which xtrsrc is stored}
\item{prec}{Number of decimal places at which the data needs to be roudned. If left empty
then function will estimate this from 10 percent of standard deviation in dependent
variable}
\item{error.low}{lower boundary of error range or uncertainty itself if error.up
is not provided (see error.up)}
\item{error.up}{upper boundary of error range (if not provided then it will assume that
error.low equals the uncertainty itself and calculates from this the upper and lower
boundary}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
\item{do.bar}{whether to plot the error bars}
\item{deplabel}{label name for the dependent variable in the visualisations}
}
\description{
Shell function for investigating the impact of model representations on compression,
reproduction of residuals, approximation of original data. This function is mainly
as a sanity check that the procedure is doing what it is supposed to do
}
\examples{
\dontrun{
studyshell(data=D,resultpath= "~/CTS/results/",CX2=0.95)
}
}
\name{getvis}
\alias{getvis}
\title{generate visualisation from the data as store these in a pdf}
\usage{
getvis(data=c(),resultpath  = c(),dependent="f_int",independent="freq_eff",
group="band",id="runcat",timecol="taustart_ts",do.log=TRUE,
do.bar=TRUE,error.low="err_low",error.up="err_up",deplabel="Intensity (Jy)")
}
\arguments{
\item{data}{dataframe based on a merge between original data and the models}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{resultpath}{directory where where the pdf with visualisations will be stored}
\item{timecol}{Name of variable to indicate time}
\item{do.log}{Turn independent and devependent variables into log}
\item{do.bar}{whether to plot the error bars}
\item{error.low}{lower boundary of error range}
\item{error.up}{upper boundary of error range}
\item{deplabel}{label name for the dependent variable in the visualisations}
}
\description{
This function takes a merge between the original data and the models to
create a visualation of both models and data
}
\value{
No output is stored
}
\examples{
\dontrun{
getvis(data=mydata2)
}
}
\name{getmodel}
\alias{getmodel}
\title{derive models from data}
\usage{
getmodel(data=c(),dependent=c(),independent=c(),group=c(),id=c(),
                    mcmc=FALSE,include_ext1=FALSE,minN=30,do.res=FALSE,CX2=0.95,
                    do.log=TRUE,imidcol="imageid",xtrsrc="xtrsrc",modeltime=TRUE)
}
\arguments{
\item{data}{dataframe}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{mcmc}{If TRUE then use Markov Chain Monte Carlo regression instead
of ordinary least squares}
\item{include_ext1}{If true then keep the datapoints in radio astronomy for which
the extract type equals 1. Extract type equals 1 data points are interpolations and
not real data. A model is likely to be more accurate if these datapoints are ommited
}
\item{minN}{minimum number of measurements required per group (frequency band)}
\item{do.res}{if true then also output the function residuals as a seperate dataframe}
\item{CX2}{Conficende interval for for Chi-square test, default = 0.95}
\item{do.log}{whether to take the log from the indpenedent and dependent variables}
\item{imidcol}{Name of column in which image id is stored}
\item{xtrsrc}{Name of column in which xtrsrc is stored}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
}
\description{
This function takes the data and derives one models per source (the subset of data
identified by id). Here, the function converts the data to log space
and then fits a linear model. This is what is assumed to work in radio astronomogy.
Additionally the function tests whether the dependent values are normally distributed
per frequency band.
}
\value{
Dataframe with models and if do.res is set to TRUE it also includes the model residuals
}
\examples{
\dontrun{
V = getmodel(data=D,dependent="f_int",independent="freq_eff",group="band",id="runcat",
                 mcmc=FALSE,include_ext1=FALSE,minN=minN,doinR=doin.R,CX2=CX2)
}
}
\name{ctsky-package}
\alias{ctsky-package}
\alias{ctsky}
\docType{package}
\title{
\packageTitle{ctsky}
}
\description{
\packageDescription{ctsky}
}
\details{
\packageDESCRIPTION{ctsky}
\packageIndices{ctsky}
The package holds a set of functions I developed as part of the
path-finding project Compressing the sky into a large number of statistical
models.\cr
\cr
The package requires that MonetDB is installed with embedded R.\cr
\cr
The function \link{studyshell} is mainly a sanity check for all the code
as it reproduces all the performance evaluations The function \link{getmodel}
is probably the most important function: It generates the models and tests for
model assumptions.
}
\author{
\packageAuthor{ctsky}
}
\examples{
\dontrun{
library(MonetDB.R)
conn = dbConnect(MonetDB.R(),host="localhost", dbname="rsm", user="monetdb",
                            password="monetdb")
D = dbGetQuery(conn,paste("SELECT a1.runcat,a1.xtrsrc,i2.freq_eff,
                            i2.taustart_ts,x1.extract_type,x1.f_int,
                            x1.f_int_err,i2.band,i2.id
                            FROM assocxtrsource a1,
                            (SELECT t1.runcat
                            FROM
                            (SELECT a.runcat, i.band, count(*)
                            FROM assocxtrsource a, runningcatalog r,
                            extractedsource x, image i
                            WHERE a.runcat = r.id and a.xtrsrc = x.id and x.image = i.id
                            GROUP BY runcat, band having count(*) > 30) t1
                            GROUP BY t1.runcat) t2, runningcatalog r1,
                            extractedsource x1, image i2
                            WHERE
                            a1.runcat = r1.id and a1.runcat = t2.runcat and
                            a1.xtrsrc=x1.id and x1.image= i2.id
                            ORDER BY a1.runcat,a1.xtrsrc",sep=""))
dbDisconnect(conn)

D = D[order(D$runcat,D$band),]
studyshell(data=D,dependent="f_int",independent="freq_eff",group="band",id="runcat",
            resultpath= "~/CTS/results/",do.res=TRUE,CX2=0.95,minN = 30)

}
}
\name{approxdata}
\alias{approxdata}
\title{approximate data based on the models as derived with function getmodel}
\usage{
approxdata(models=c(),do.log=TRUE, dependent="f_int",independent="freq_eff",
group="band",id="runcat",modeltime=TRUE)
}
\arguments{
\item{models}{dataframe of models as produced by function getmodel}
\item{do.log}{Turn independent and devependent variables into log}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
}
\description{
This function takes the models dataframe as produced by getmodel
and uses it to replicate the original data based on the model, which will
be an approximation.
}
\value{
A dataframe with the same size as the dataframe from which models was derived
}
\examples{
\dontrun{
rep.data = approxdata(models=models)
}
}
\name{testcompr}
\alias{testcompr}
\title{test data compression}
\usage{
testcompr(data=c(),path=c())
}
\arguments{
\item{data}{vector of data}
\item{path}{directory that is used for temporarily creating files}
}
\description{
This function takes a vector and then tests the compression of the
vector using nine different compression techniques:
 gzip low (level 1), high (level 9) or default (level NA)
 bzip2 low (level 1), high (level 9) or default (level NA)
 xz low (level 1), high (level 9) or default (level NA)
}
\examples{
\dontrun{
result = testcompr(data=D,path=filepath)
}
}
\name{studyshell}
\alias{studyshell}
\title{shell function for studying the data}
\usage{
studyshell(data=c(), resultpath= c(),CX2=0.95,minN = 30,
                      dependent="f_int",independent="freq_eff",group="band",id="runcat",
                      timecol="taustart_ts",mcmc=FALSE,include_ext1=FALSE,do.res=FALSE,
                      do.log=TRUE,imidcol="imageid",xtrsrc="xtrsrc",
                      error.low=c(),error.up=c(),prec=c(),
                      modeltime=TRUE,do.bar=TRUE,
                      deplabel="Intensity (Jy)")
}
\arguments{
\item{data}{dataframe including the variables (columns): source identifiers (runcat, xtrsrc),
independent variable (freq_eff), timestamps (taustart_ts), data type (extract_type),
independent variables (f_int), 1-sigma error in the independent variable (f_int_err),
grouping variable (band), image id (id)
}
\item{resultpath}{directory where code stores the results}
\item{CX2}{Confidence interval for the chi-square test}
\item{minN}{Minimum number of measurements per frequency band}
\item{dependent}{Name of the dependent variable in the model (f_int in our case)}
\item{independent}{Name of the independent variable in the model (freq_eff in our case)}
\item{group}{Grouping variable (band in our case)}
\item{id}{Identifier of the image, this is mainly used as a key to aid tracing back the
origin of the data}
\item{timecol}{Name of variable to indicate time}
\item{mcmc}{If TRUE then use Markov Chain Monte Carlo regression instead
of ordinary least squares}
\item{include_ext1}{If true then keep the datapoints in radio astronomy for which
the extract type equals 1. Extract type equals 1 data points are interpolations and
not real data. A model is likely to be more accurate if these datapoints are ommited}
\item{do.res}{if true then also output the function residuals as a seperate dataframe}
\item{do.log}{whether to take the log from the indpenedent and dependent variables}
\item{imidcol}{Name of column in which image id is stored}
\item{xtrsrc}{Name of column in which xtrsrc is stored}
\item{prec}{Number of decimal places at which the data needs to be roudned. If left empty
then function will estimate this from 10 percent of standard deviation in dependent
variable}
\item{error.low}{lower boundary of error range or uncertainty itself if error.up
is not provided (see error.up)}
\item{error.up}{upper boundary of error range (if not provided then it will assume that
error.low equals the uncertainty itself and calculates from this the upper and lower
boundary}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
\item{do.bar}{whether to plot the error bars}
\item{deplabel}{label name for the dependent variable in the visualisations}
}
\description{
Shell function for investigating the impact of model representations on compression,
reproduction of residuals, approximation of original data. This function is mainly
as a sanity check that the procedure is doing what it is supposed to do
}
\examples{
\dontrun{
studyshell(data=D,resultpath= "~/CTS/results/",CX2=0.95)
}
}
\name{getvis}
\alias{getvis}
\title{generate visualisation from the data as store these in a pdf}
\usage{
getvis(data=c(),resultpath  = c(),dependent="f_int",independent="freq_eff",
group="band",id="runcat",timecol="taustart_ts",do.log=TRUE,
do.bar=TRUE,error.low="err_low",error.up="err_up",deplabel="Intensity (Jy)")
}
\arguments{
\item{data}{dataframe based on a merge between original data and the models}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{resultpath}{directory where where the pdf with visualisations will be stored}
\item{timecol}{Name of variable to indicate time}
\item{do.log}{Turn independent and devependent variables into log}
\item{do.bar}{whether to plot the error bars}
\item{error.low}{lower boundary of error range}
\item{error.up}{upper boundary of error range}
\item{deplabel}{label name for the dependent variable in the visualisations}
}
\description{
This function takes a merge between the original data and the models to
create a visualation of both models and data
}
\value{
No output is stored
}
\examples{
\dontrun{
getvis(data=mydata2)
}
}
\name{getmodel}
\alias{getmodel}
\title{derive models from data}
\usage{
getmodel(data=c(),dependent=c(),independent=c(),group=c(),id=c(),
                    mcmc=FALSE,include_ext1=FALSE,minN=30,do.res=FALSE,CX2=0.95,
                    do.log=TRUE,imidcol="imageid",xtrsrc="xtrsrc",modeltime=TRUE)
}
\arguments{
\item{data}{dataframe}
\item{dependent}{name of dependent variable as stored in data (source intensity)}
\item{independent}{name of independent variable as stored in data (frequency)}
\item{group}{name of grouping variable, band (frequency band) in the case
radio astronomy}
\item{id}{name of variable as stored in data for identifying the subsets of
data for which a model needs to be developed. In radio astronomy this is the source
identifier named runcat}
\item{mcmc}{If TRUE then use Markov Chain Monte Carlo regression instead
of ordinary least squares}
\item{include_ext1}{If true then keep the datapoints in radio astronomy for which
the extract type equals 1. Extract type equals 1 data points are interpolations and
not real data. A model is likely to be more accurate if these datapoints are ommited
}
\item{minN}{minimum number of measurements required per group (frequency band)}
\item{do.res}{if true then also output the function residuals as a seperate dataframe}
\item{CX2}{Conficende interval for for Chi-square test, default = 0.95}
\item{do.log}{whether to take the log from the indpenedent and dependent variables}
\item{imidcol}{Name of column in which image id is stored}
\item{xtrsrc}{Name of column in which xtrsrc is stored}
\item{modeltime}{Whether to use constant of time series as a model (TRUE) or
to use one model per source based on the assuming of a exponential relationship between
source intensity and frequency (FALSE)}
}
\description{
This function takes the data and derives one models per source (the subset of data
identified by id). Here, the function converts the data to log space
and then fits a linear model. This is what is assumed to work in radio astronomogy.
Additionally the function tests whether the dependent values are normally distributed
per frequency band.
}
\value{
Dataframe with models and if do.res is set to TRUE it also includes the model residuals
}
\examples{
\dontrun{
V = getmodel(data=D,dependent="f_int",independent="freq_eff",group="band",id="runcat",
                 mcmc=FALSE,include_ext1=FALSE,minN=minN,doinR=doin.R,CX2=CX2)
}
}
\name{ctsky-package}
\alias{ctsky-package}
\alias{ctsky}
\docType{package}
\title{
\packageTitle{ctsky}
}
\description{
\packageDescription{ctsky}
}
\details{
\packageDESCRIPTION{ctsky}
\packageIndices{ctsky}
The package holds a set of functions I developed as part of the
path-finding project Compressing the sky into a large number of statistical
models.\cr
\cr
The package requires that MonetDB is installed with embedded R.\cr
\cr
The function \link{studyshell} is mainly a sanity check for all the code
as it reproduces all the performance evaluations The function \link{getmodel}
is probably the most important function: It generates the models and tests for
model assumptions.
}
\author{
\packageAuthor{ctsky}
}
\examples{
\dontrun{
library(MonetDB.R)
conn = dbConnect(MonetDB.R(),host="localhost", dbname="rsm", user="monetdb",
                            password="monetdb")
D = dbGetQuery(conn,paste("SELECT a1.runcat,a1.xtrsrc,i2.freq_eff,
                            i2.taustart_ts,x1.extract_type,x1.f_int,
                            x1.f_int_err,i2.band,i2.id
                            FROM assocxtrsource a1,
                            (SELECT t1.runcat
                            FROM
                            (SELECT a.runcat, i.band, count(*)
                            FROM assocxtrsource a, runningcatalog r,
                            extractedsource x, image i
                            WHERE a.runcat = r.id and a.xtrsrc = x.id and x.image = i.id
                            GROUP BY runcat, band having count(*) > 30) t1
                            GROUP BY t1.runcat) t2, runningcatalog r1,
                            extractedsource x1, image i2
                            WHERE
                            a1.runcat = r1.id and a1.runcat = t2.runcat and
                            a1.xtrsrc=x1.id and x1.image= i2.id
                            ORDER BY a1.runcat,a1.xtrsrc",sep=""))
dbDisconnect(conn)

D = D[order(D$runcat,D$band),]
studyshell(data=D,dependent="f_int",independent="freq_eff",group="band",id="runcat",
            resultpath= "~/CTS/results/",do.res=TRUE,CX2=0.95,minN = 30)

}
}
.. image:: https://zenodo.org/badge/47755645.svg
   :target: https://zenodo.org/badge/latestdoi/47755645

Compressing the Sky
====================

Path-finding project with CWI Database Architectures Group
on compress and handling astronomical data (Lofar)
Project duration: June 2015-December 2015

Guidance on folder structure:
-----------------------------
- ctsky_0.2.tar.gz in the root directory is the R-package source code
- /ctsky includes all files from which ctsky_0.2.tar.gz is compiled
- ctsky.Rcheck/ctsky-manual.pdf is the package manual
- due to file size i have not uploaded any example data
- embedded_sql.sh is an example of how to embedded the package in MonetDB SQL code
- dfspt-explore.R and lofar-explore.R are example files for exploring the astronomy data from within R based on calles to the ctsky-package

Notes on how to install MonetDB with R integration:
-------------------------------------------
- check MonetDB website for latest updates and instructions, these are the commands I used
- hg clone http://dev.monetdb.org/hg/MonetDB/ MonetDB
- cd MonetDB
- hg update Jul2015
- ./bootstrap
- ./configure —prefix=/some/install/dir —enable-rintegration=yes —enable-optimize
- make -j clean install
- into ~/.bashrc: export PATH=$PATH:/some/install/dir/bin

