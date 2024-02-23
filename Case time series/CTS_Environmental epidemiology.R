#A case study for applications in environmental epidemiology (09 August 2021, Antonio Gasparrini)

#library
library(dlnm) ; library(gnm) ; library(data.table) ; library(splines)

set.seed(13041975)
par(las=1)

#Simulating the original data
n <- 1601
dstart <- as.Date("2015-10-29")
dend   <- as.Date("2018-11-19")
date   <- seq(dstart, dend, by=1)
year   <- year(date)
month  <- month(date)
doy    <- yday(date)
dow    <- factor(wday(date))

fustart <- sample(seq(dstart, dend-10, by=1), n, replace=TRUE)
fuend   <- fustart + pmax(pmin(round(exp(rnorm(n, 5.1, 2))), dend-fustart), 10)
sum(fuend-fustart+1)

#While the follow-up distribution does not match perfectly the original study, the sampling parameters are set
#above to generate approximately the same number of total person-days, in this case, 363,901.
#Finally, we define some variables used to simulate the distribution of the environmental exposures and, later,
#the seasonal baseline risk. These variables are the cosine transformation of doy and quadratic splines of
#date with 5 degrees of freedom per year. In addition, we simulate 20 random smoke days occurring in the
#(Australian) summer. The code:
  
cosdoy <- cos(doy*2*pi / 366)
spldate <- bs(date, degree=2, int=TRUE, df=round(length(date)/365.25)*5)
smokeday <- date %in% sample(date[month %in% c(1,2,12)], 20)

#We are now ready to simulate the distribution of the three environmental stressors. In the original study,
#individual exposure series were reconstructed through the geo-location system of the smartphone by linkage
#with detailed spatio-temporal exposure maps. In order to simplify the simulation process, we derive here
#a single series for each stressor, assuming that all the 1601 subjects are exposed to the same levels on
#the same day. This does not affect the generality of the example, and in real-case settings, individual-level
#exposure series can nevertheless be used.
#The environmental exposures are created by assuming an underlying seasonal trend, represented by the
#cosine variable above, plus auto-correlated random normal deviations. Exponentiation is used to produce
#non-negative values of pollen (grains/m3) and pollution (PM2.5, 𝜇gr/m3), while temperature (∘C) is sampled
#directly. The code:

pollen  <- exp(cosdoy*2+2.5 + arima.sim(list(ar=0.5), length(date), sd=0.8))
pm      <- exp((-cosdoy)*1.6+2.5 + smokeday*3.2 +
            arima.sim(list(ar=0.6), length(cosdoy), sd=0.95))
tmean   <- cosdoy*6+15 + arima.sim(list(ar=0.6), length(cosdoy), sd=2.6)
envdata <- data.frame(date, pollen, pm, tmean)

#The variables are included in the dataframe envdata. 
#The definitions above provide a realistic distribution of the three exposures, 
#with pollen and temperature peaking in summer, 
#while PM2.5 shows higher wintertime levels but with isolated spikes in the summer corresponding to smoke days due to fires. 
#A visual representation is offered by the plots obtained through:
  
x11();par(mfrow=c(1,3))
plot(date, pollen, xlab="Date", ylab="Pollen", main="Pollen series", col=2,
     pch=19)
plot(date, pm, xlab="Date", ylab="PM2.5", main="Pollution series", col=2,
     pch=19)
plot(date, tmean, xlab="Date", ylab="Temperature", main="Temperature series",
     col=2, pch=19)

#The variables created above can now be used to define individual risk profiles of experiencing allergic symptoms. 
#These profiles will be simulated as risks associated to the three exposures on top of baseline trends.
#We first simulate the latter as a combination of shared underlying risks and individual-level deviations:

fortrend <- function(ind=TRUE) (cosdoy*1.6 + sin(doy*4*pi/366))/8+1 + if(ind)
  spldate %*% runif(ncol(spldate),-0.2,0.2) else 0

x11();plot(date, fortrend(ind=F), type="l", lwd=2, ylim=c(0.5,1.5), xlab="Date",
           ylab="OR", main="Shared baseline risk and indivdual deviations")
abline(h=1)
for(i in 1:7) lines(date, fortrend(ind=T), type="l", lty=2, col=i)


#The next step is the definition of the increase in risk due to exposure to the three environmental stressors.
#Specifically, we define non-linear relationships for pollen and temperature, with effects lagged up to 3 days,
#and a linear and unlagged association with PM2.5. First, we define the three functions to specify the three
#exposure-response risk shapes and the lag structure:
  
forpoll  <- function(x) 1.6 - 0.6*exp(-x/60)
forpm    <- function(x) exp(x/1000)
fortmean <- function(x) 1 + ifelse(x>15, 0.002*(x-15)^2, 0)
fwlag    <- function(lag) exp(-lag/1.5)

#These functions define relationships in the OR and lag scales, and can be represented graphically with:

x11();par(mfrow=c(1,4))
plot(0:200, forpoll(0:200), type="l", xlab="Pollen", ylab="OR",
     main="Exposure-response with pollen", col=2)
plot(0:100, forpm(0:100), type="l", xlab="PM2.5", ylab="OR",
     main="Exposure-response with PM2.5", col=2)
plot(0:30, fortmean(0:30), type="l", xlab="Temperature", ylab="OR",
     main="Exposure-response with temperature", col=2)
plot(0:50/10, fwlag(0:50/10), type="l", xlab="Lag", ylab="Weight",
     main="Lag structure", col=2)

#These shapes are similar to the associations estimated in the original study (Gasparrini 2021). 
#The lag structure is defined as weights, and can be used to represent a decreasing OR proportionally to time after
#the exposure occurred. As an example, we used the functions above to calculate the net OR in a given day
#after exposures to pollen of 50, 9, 135, and 93 grains/m3 in the same and past 3 days (lag 0–3):

exp(sum(log(forpoll(c(50,9, 135, 93))) * fwlag(0:3)))

#Simply, the expression above computes the log-OR for each exposure occurrence, which are then weighted
#depending on the lag and then summed and exponentiated to obtain the overall cumulative OR .
#We can now apply the same computation to the whole series for the three exposures, using first the function
#exphist to generate the matrix of lagged exposures, and then applying the expression for each row:
  
orpoll <- apply(exphist(pollen, lag=3), 1, function(x) exp(sum(log(forpoll(x)) * fwlag(0:3))))
orpm <- forpm(pm)

ortmean <- apply(exphist(tmean, lag=3), 1, function(x)
  exp(sum(log(fortmean(x)) * fwlag(0:3))))
orenv <- orpoll * orpm * ortmean

#Note that we assume a lag 0–3 for pollen and temperature, while we simply define a same-day association
#with no lag for PM2.5. The code above computes therefore the OR contribution for each exposure in each
#day, which are then multiplied to obtain the overall risk associated with all the three environmental stressors
#in the vector orenv. The series for temperature and all the exposures are graphically represented below:

x11();par(mfrow=c(1,2))
plot(date, ortmean, type="l", xlab="Date", ylab="OR",main="Risk of temperature")
plot(date, orenv, type="l", xlab="Date", ylab="OR",main="Risk of all environmental stressors")

#We have now all the information required for simulating the original data. These are created by looping in a
#list, producing the observations for each subject, and then binding them in a dataframe. Each of the blocks
#of code in the loop performs the following steps for each subject:
#1. Define the follow-up period and identify the related subset of the study period
#2. Create the total risk contribution in each follow-up day
#3. Sample the occurrence of respiratory symptoms within the follow-up period
#4. Put the information together in a dataframe, adding the series of environmental exposures

#Here is the R code:

dlist <- lapply(seq(n), function(i) {
  fudate <- seq(fustart[i], fuend[i], by=1)
  sub <- date %in% fudate
  ortot <- fortrend(ind=T)[sub] * orenv[sub] * (1 + wday(fudate) %in% c(2:6)*0.4)
  pbase <- plogis(-3.3 + 14/length(fudate) - 0.0015*length(fudate))
  sympt <- rbinom(sum(sub), 1, plogis(qlogis(pbase) + log(ortot)))
  data <- cbind(data.frame(id=paste0("sub",sprintf("%04d", i)), date=fudate,
                           year=year[sub], month=month[sub], dow=dow[sub], y=sympt), envdata[sub, -1])
  return(data.table(data))
})
data <- do.call(rbind, dlist)

#Analysis
dsub <- subset(data, id=="sub0036")

x11();par(mfrow=c(2,2))
plot(y~date, data=dsub, type="h", lty=2, ylim=c(0,2), yaxt="n", bty="l", xlab="",
     ylab="Day with \nsymptoms", mgp=c(2.2,0.7,0), lab=c(5,3,7))
points(y~date, data=subset(dsub,y>=1), pch=23, bg=3)
plot(pollen~date, data=dsub, type="l", lty=1, bty="l", col=2, xlab="",
     ylab="Pollen")
plot(pm~date, data=dsub, type="l", lty=1, bty="l", col=4, xlab="Date",
     ylab="PM2.5")
plot(tmean~date, data=dsub, type="l", lty=1, bty="l", col=3, xlab="Date",
     ylab="Temperature")


#We can now define the different terms to be included in the regression model. 
#First, we define a set of splines of time with approximately 8 degrees of freedom per year, 
# and subject/year/month strata indicators,
#to be used to model the shared seasonal trend and individual deviations, respectively.:
  
dftrend      <- round(as.numeric(diff(range(data$date))/365.25 * 8))
btrend       <- ns(data$date, knots=equalknots(data$date, dftrend-1))
data$stratum <- with(data, factor(paste(id, year, month, sep="-")))

#We now apply the function crossbasis() to parameterise distributed lag linear and 
#non-linear transformations of the environmental variables:
  
cbpoll  <- crossbasis(data$pollen, lag=3, argvar=list(knots=c(40,100)),
                     arglag=list(knots=1), group=data$id)
cbpm    <- crossbasis(data$pm, lag=3, arglag=list("integer"), group=data$id)
cbtmean <- crossbasis(data$tmean, lag=3, argvar=list(knots=1:2*10),
                      arglag=list(knots=1), group=data$id)

mod <- gnm(y ~ cbpoll + cbpm + cbtmean + btrend + dow, eliminate=stratum, data=data,
           family=binomial)

#The estimated coefficients and associated (co)variance matrix of the model can now be used to predict the
#association of the various terms with the risk of respiratory symptoms, using the function crosspred():

cppoll  <- crosspred(cbpoll , mod, at=0:20*10, cen=0)
cppm    <- crosspred(cbpm   , mod, at=0:20*5, cen=0)
cptmean <- crosspred(cbtmean, mod, cen=15, by=1.5)

#We can now represent graphically the association in both dimensions of exposure intensity and lag. 
#Specifically, the plots below represent the overall cumulative exposure-responses 
#(interpreted as the net associations accounting for the whole lag period), 
#the full bi-dimensional exposure-lag-responses for non-linear relationships of pollen and temperature, 
#and the lag-response corresponding to a 10𝜇gr/m3 increases in PM2.5. The code:
  
x11();par(mfrow=c(1,2))
plot(cppoll, "overall", xlab="Pollen", ylab="OR", col=2,
     main="Pollen: overall cumulative exposure-response", ylim=c(0.5,3))
plot(cppoll, xlab="Pollen", zlab="OR", main="Pollen: exposure-lag-response",
     cex.axis=0.8, col=2)

x11();par(mfrow=c(1,2))
plot(cppm, var=10, "overall", xlab="PM2.5", ylab="OR", col=4,
     main="PM2.5: overall cumulative exposure-response", ylim=c(0.95,1.20))
plot(cppm, var=10, ci="b", type="p", ylab="OR", col=4, pch=19, cex=1.7,
     xlab="Lag", main="PM2.5: lag-response", lab=c(3,5,7), ylim=c(0.995,1.015))

x11();par(mfrow=c(1,2))
plot(cptmean, "overall", xlab="Temperature", ylab="OR", col=3,
     main="Temperature: overall cumulative exposure-response", ylim=c(0.5,3))
plot(cptmean, xlab="Temperature", zlab="OR", ltheta=240, lphi=60, cex.axis=0.8,
     main="Temperature: exposure-lag-response", col=3)
