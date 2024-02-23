#A case study for applications in clinical epidemiology (09 August 2021, Antonio Gasparrini)
#librarry
pacman::p_load(dlnm,gnm,mgcv,pbs,data.table,scales)

set.seed(13041975)
par(las=1)

#Simulating the original data
n <- 3927
dstart <- as.Date("2007-01-01")
dend <- as.Date("2007-12-31")

#Then we generate the time variables across the follow-up period, namely date (calendar days), 
#time (a sequence of integers starting from 1), month (months in numbers), and doy (days of the year). 
#In addition,we randomly generate dob (date of birth) for each subject, 
#with age at start between 35 and 100 years old.

date <- seq(dstart, dend, by=1)
times <- seq(length(date))
month <- month(date)
doy <- yday(date)
dob <- sample(seq(dstart-round(100*365.25), dstart-round(35*365.25), by=1), n)

#These variables are used for simulating the temporal variation in the underlying risk of AMI, 
#with a cyclic seasonal trend and a long-term change by age modelled by a cosine function and polynomials, respectively.
#These effects are defined as a incident rate ratio (IRR), and created by the following code:
  
frrseas <- function(doy) (cos(doy*2*pi / 366) + 1) / 12 + 1
frrage  <- function(age) exp((age - 70) / 6)

#These temporal variations in risk along day of the year and age are represented in the graphs below:
x11();par(mfrow=c(1,2))
plot(1:365, frrseas(1:365), type="l", col=2, ylab="IRR", xlab="Day of the year",
     main="Simulated seasonal effect")
plot(35:90, frrage(35:90), type="l", col=2, ylab="IRR", xlab="Age",
     main="Simulated age effect")

#Now we create a function to define the IRR along the lag dimension. 
#In this case, this dimension represents the risk after a flu episode, 
#with the lag unit defined by day. Similarly, we illustrate the phenomenon graphically:

frrlag <- function(lag) exp(-(lag/10)) * 4 + 1
x11();plot(1:91, frrlag(1:91), type="l", col=2, ylab="IRR", xlab="Lag",
     main="Simulated lag effect")

#The graph indicates that, within a lag period of 3 months (1 to 91 days of lag) as in the original analyses,
#the risk is much increased in the first days after the flu episode, 
#but then it attenuates and tends to null after approximately one month.
#In the presence of multiple exposure episodes, lagged effects can cumulate in time, depending on the
#exposure profile of an individual. In this case, the risk at a given day is determined by the exposure history
#to flu, with potentially multiple flu episodes contributing at different lags for the same day. 
#For instance, the code below shows an example with the risk associated with four flu episodes in a 250-day period, 
#with the cumulated risk being the product of lag-specific contributions:

expprof <- as.numeric(seq(250) %in% c(15,100,110,160))
exphist <- exphist(expprof, lag=c(1,91), fill=0)
rrflu <- apply(exphist, 1, function(x) prod(frrlag(1:91)[x==1]))

x11();plot(seq(250), rrflu, type="l", col=2, ylab="Overall cumulative IRR", xlab="Days",
     main="Example of cumulated effect")
points(c(15,100,110,160), rep(1,4), pch=19)
legend("topright", c("Flu episodes","Cumulated IRR"), pch=c(19,NA), lty=c(NA,1),
       col=1:2, cex=0.8)


#We have now all the information required for simulating the original data. 
#These will consist of individual records with the following variables, with age measured in days:
#• id: the identifier of the subject
#• dob: date of birth
#• start: the age of the subject at the start of follow-up
#• end: the age of the subject at the end of follow-up
#• event: the age of the subject at the occurrence of the AMI event
#• flu*: multiple variables defining the age(s) of the subject at each flu episode

#The data are simulated by looping in a list, producing the observations for each subject, and then binding
#them in a dataframe. Each of the blocks of code in the loop performs the following steps for each subject:
#1. Sample the number of flu episode(s); define the risk of having a flu episode in each day; sample the
#   flu episodes and create an exposure profile
#2. Create the exposure history of flu for each day for a given lag period; compute the overall cumulative
#   AMI risk due to flu for each day
#3. Define the total AMI risk for each day, dependent on age, season, and flu; sample the unique AMI event
#4. Put the information together in a dataframe; add the flu episodes, setting them to NA if less than the
#   sampled maximum of 10
#Here is the R code (it takes less than a minute):

dlist <- lapply(seq(n), function(i) {
  nflu <- rpois(1,1) + 1
  expprof <- drop(rmultinom(1, nflu, frrseas(doy))) > 0 + 0
  exphist <- exphist(expprof, lag=c(1,91), fill=0)
  rrflu <- apply(exphist, 1, function(x) prod(frrlag(1:91)[x==1]))
  rrtot <- frrage(as.numeric((date-dob[i])/365.25)) * frrseas(doy) * rrflu
  devent <- date[drop(rmultinom(1, 1, rrtot))==1]
  data <- data.frame(id = paste0("sub", sprintf("%03d", i)), dob = dob[i],
                     start = as.numeric(dstart - dob[i]), end = as.numeric(dend - dob[i]),
                     event = as.numeric(devent - dob[i]))
  flu <- as.numeric(date[expprof == 1] - dob[i])
  for(j in seq(10)) data[paste0("flu", j)] <- if(j>nflu) NA else flu[j]
  return(data)
})
dataorig <- do.call(rbind, dlist)
dataorig

#Data expansion
sub <- dataorig[3,]

date    <- as.Date(sub$start:sub$end, origin=sub$dob)
datasub <- data.frame(
  id = sub$id,
  date = date,
  times = seq(length(date)),
  age = as.numeric(date-sub$dob)/365.25,
  y = as.numeric(date-sub$dob) %in% sub$event + 0,
  flu = as.numeric(date-sub$dob) %in% na.omit(as.numeric(sub[6:15])) + 0,
  month = month(date),
  doy = yday(date)
)
head(datasub)
datasub

exphistsub <- exphist(datasub$flu, lag=c(1,91), fill=0)
timeflu1   <- sub$flu1-sub$start+1
exphistsub[timeflu1 + 0:5, 1:10]

dlist <- lapply(seq(n), function(i) {
  sub <- dataorig[i,]
  date <- as.Date(sub$start:sub$end, origin=sub$dob)
  data <- data.frame(
    id = sub$id,
    date = date,
    times = seq(length(date)),
    age = as.numeric(date-sub$dob)/365.25,
    y = as.numeric(date-sub$dob) %in% sub$event + 0,
    flu = as.numeric(date-sub$dob) %in% na.omit(as.numeric(sub[6:15])) + 0,
    month = month(date),
    doy = yday(date)
  )
  exphist <- exphist(data$flu, lag=c(1,91), fill=0)
  return(data.table(cbind(data, exphist)))
})
data <- do.call(rbind, dlist)


x11();plot(unique(data$date), unique(data$date), ylim=c(0.5,5+0.5), yaxt="n",
     ylab="", xlab="Follow-up", frame.plot=F)
axis(2, at=5:1, labels=paste("Sub",5:1), lwd=0, las=1)
for(i in 5:1) {
  sub <- subset(data, id==unique(data$id)[i])
  flu <- sub$date[sub$flu==1]
  rect(flu+1, rep(i-0.3,length(flu)), flu+91, rep(i+0.3,length(flu)), border=NA,
       col=alpha("gold3",0.3))
  lines(sub$date, rep(i, nrow(sub)))
  points(sub$date[sub$y==1], i, pch=21, bg=2, cex=1.5)
}
View(data)
table(subset(data,id %in% c("sub001","sub002","sub003","sub004","sub005"))$flu,
      subset(data,id %in% c("sub001","sub002","sub003","sub004","sub005"))$id)


#Now, we replicate the main case time series analysis illustrated in the original article (Gasparrini 2021). 
#We first derive the terms to control for age and season using natural cubic and cyclic splines, respectively. 
#We use the wrapper function onebasis() that simplifies the prediction and plotting of these associations, 
#to be performed later. We call the function pbs from the package with the same name to generate the basis
#transformations for the cyclic splines. The code:
  
splage  <- onebasis(data$age, "ns", knots=quantile(data$age, c(1,3)*0.25))
splseas <- onebasis(data$doy, "pbs", df=3)


exphist <- data[,-c(1:8)]
cbspl <- crossbasis(exphist, lag=c(1,91), argvar=list("strata",breaks=0.5),
                    arglag=list("ns",knots=c(3,10,29)))

#We now have all the terms for fitting the fixed-effects Poisson regression using the function gnm(). The
#regression model includes all the predictors, and defines the conditional stratification through the argument
#eliminate. This is the code:
mspl <- gnm(y ~ cbspl+splage+splseas, data=data, family=poisson,
            eliminate=factor(id))

cpspl     <- crosspred(cbspl  , mspl, at=1)
cpsplage  <- crosspred(splage , mspl, cen=70   , to=90)
cpsplseas <- crosspred(splseas, mspl, cen=366/2, at=1:365)


x11();par(mfrow=c(1,3))
plot(cpsplage, col=2, ylab="IRR of AMI", xlab="Age", main="Risk by age",ylim=c(0,25))
mtext("Natural cubic splines with 3df", cex=0.8)
plot(cpsplseas, col=2, ylab="IRR of AMI", xlab="Day of the year",main="Risk by season", ylim=c(0.95,1.30))
mtext("Cyclic splines with 4df", cex=0.8)
plot(cpspl, var=1, col=2, ylab="IRR of AMI", xlab="Days after a flu episode",ylim=c(0,5), main="Risk by lag")
mtext("Natural cubic splines with 5df (DLM)", cex=0.8)

#An alternative parameterisation of cross-basis term can be used, 
#specifically using strata functions to represent the risk along lags. The code:
cbstr <- crossbasis(exphist, lag=c(1,91), argvar=list("strata",breaks=0.5),
                    arglag=list("strata",breaks=c(4,8,15,29)))

#We can now fit the alternative model:
mstr  <- gnm(y ~ cbstr+splage+splseas, data=data, family=poisson,eliminate=factor(id))
cpstr <- crosspred(cbstr, mstr, at=1)

x11();plot(cpstr, var=1, col=3, ylab="IRR of AMI", xlab="Days after a flu episode",ylim=c(0,5), main="Risk by lag")
mtext("Strata of lag periods (DLM)", cex=1)
lines(cpspl, var=1, lty=2)

resstr <- round(with(cpstr, t(rbind(matRRfit,matRRlow,matRRhigh))), 2)
colnames(resstr) <- c("IRR","low95%CI","high95%CI")
resstr[paste0("lag", c(1,4,8,15,29)),]
