#####################################################################################################
##### EGYPTIAN VULTURE POPULATION TREND IN BALKAN COUNTRIES #################
#####################################################################################################
## Velevski, M., Nikolov, S. C., Hallmann, B., Dobrev, V., Sidiropoulos, L., Saravia, V., Tsiakiris, R., Arkumarev, V., Galanaki, A., Kominos, T., Stara, K., Kret, E., Grubač, B., Lisičanec, E., Kastritis, T., Vavylis, D., Topi, M., Hoxha, B. & Oppel, S. 2015. Population decline and range contraction of the Egyptian Vulture Neophron percnopterus on the Balkan Peninsula. Bird Conservation International, 25(4): 440-450.
## article available at http://lifeneophron.eu/files/docs/1449818562_7.pdf
## SSM model based on Kery and Schaub 2012, Ch 5
## modified 31 January 2014 to include two separate trend estimates - from 2003-2013 (better data) and for longest time series
## modified 20 May 2014 to include loop over all 3 countries in one model
## excluded one estimate from Greece that was far out (45 pairs)

setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\EGVU_Balkans")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\submitted\\EGVU_Balkans")

bugs.dir<-("C:\\STEFFEN\\WinBUGS14")
bugs.dir<-("C:\\Program Files\\WinBUGS14")
library(R2WinBUGS)
library(R2jags)


################################## LOAD DATA  ######################################
EGVU<-read.table("Balkan_EV.csv", header=T, sep=",")
n.years<-dim(EGVU)[1]		## defines the number of years
#EGVU[31,4]<-NA			## remove the outlier of 45 pairs in Greece
#EGVU[c(26:30,32:34),4]<-NA			## remove the bad data for Greece

################################## SIMPLE STATE SPACE MODEL  ######################################

###################### WIN BUGS MODEL ####################################################


###############  Specify model in BUGS language ############
sink("EGVU_Balkantrend.txt")
cat("
model{

# Priors and constraints

for (s in 1:countries){			### start loop over every country
N.est[1,s] ~ dunif(0,500)   ## draw random value from a uniform distribution between 0 and 500 for initial population size 
mean.lambda[s] ~ dunif(0,10)	#Prior for mean growth rate
sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
sigma2.proc[s]<-pow(sigma.proc[s],2)
tau.proc[s]<-pow(sigma.proc[s],-2)
sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
sigma2.obs[s]<-pow(sigma.obs[s],2)
tau.obs[s]<-pow(sigma.obs[s],-2)
}


##### Likelihood function

for (s in 1:countries){			### start loop over every country
## State process for entire time series

for (t in 1:(T-1)){
lambda[t,s] ~ dnorm(mean.lambda[s], tau.proc[s])								# Distribution for random error of growth rate
N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										# Linear predictor (population size based on past pop size and change rate)
}														# run this loop over i= 40 years


## Observation process

for (t in 1:T){
y[t,s] ~ dnorm(N.est[t,s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
}														# run this loop over i= 40 years
}		## end country loop

## Derived parameters
for (t in 1:T){
pop.size[t]<-sum(N.est[t,])
}
mlam <- mean(lambda[1:(T-1),])  				# Arithmetic mean for whole time series
mlam10 <- mean(lambda[(T-10):(T-1),])   				# Arithmetic mean for 2003-2013
Balkan.lam<-mean(mean.lambda[])

for (s in 1:countries){
mlam10country[s] <- mean(lambda[(T-10):(T-1),s])   				# Arithmetic mean for 2003-2013 for each country
}

}														# close the model loop
", fill=TRUE)
sink()
###############  End of model in BUGS language ############




EV<-EGVU[EGVU$year>=1980,c(2,3,4)]
year<-seq(1980,2013,1)
n.years=length(year)
length(EV)
countries<-3


## Bundle data
bugs.data<-list(y=as.matrix(EV), T=n.years, countries=3)		## year is standardized so that covariate values are not too far away from zero

## Initial values
inits<- function() {list(sigma.proc=runif(countries,0,5),
mean.lambda=runif(countries,0.1,2),
sigma.obs=runif(countries,0,10),
N.est=matrix(data=(c(runif(countries,150,300),rep(NA,(n.years-1)*countries))), nrow=n.years,byrow=T))}


## Parameters to be estimated ('monitored') by WinBUGS
params<-c("lambda", "N.est","pop.size", "mean.lambda", "mlam", "Balkan.lam","mlam10","mlam10country")


## MCMC settings
ni<-5000000  	## number of iterations (draws per chain)
nt<-25		## thinning rate
nb<-2500000	## length of burn-in (number of iterations discarded at the start of each chain until convergence is reached)
nc<-3		## number of chains


## call WinBUGS from R   ###
## label model with name of country:
Balkan<- jags(data=bugs.data, inits=inits, parameters.to.save=params, model.file="EGVU_Balkantrend.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)


## summarize posteriors
sink("EGVU_trend_summaries2.txt")
print(Balkan,dig = 3)
sink()





## DRAW FIGURE AFTER RUNNING ALL THREE MODELS FOR GREECE, BULGARIA, AND MACEDONIA


EV.ssm<-Balkan$BUGSoutput
fitted<-lower<-upper<-Greece<-FYROM<-Bulgaria<-as.numeric()
for (i in 1:(n.years)){
fitted[i]<-quantile(EV.ssm$sims.list$pop.size[,i], 0.5)
lower[i]<-quantile(EV.ssm$sims.list$pop.size[,i], 0.025)
upper[i]<-quantile(EV.ssm$sims.list$pop.size[,i], 0.975)

Greece[i]<-quantile(EV.ssm$sims.list$N.est[,i,3], 0.5)
FYROM[i]<-quantile(EV.ssm$sims.list$N.est[,i,2], 0.5)
Bulgaria[i]<-quantile(EV.ssm$sims.list$N.est[,i,1], 0.5)
}
top<-max(upper)+7

pdf("EV_population_trend_all_Balkan.pdf", width=9, height=6)
postscript("Fig1_Balkan.eps", width=9, height=6)
#jpeg("Fig1_Balkan.jpg", width=9, height=6, units="in", res=600, quality=100)

library(Hmisc)
par(oma=c(0,0,0,0),mar=c(4,4,0,0), cex=1.2)
errbar(c(1980:2013),fitted, lower,upper,ylim=c(0,700), xlim=c(1980, 2013), ylab="Egyptian Vulture breeding pairs", xlab="Year", las=1, type='p', pch=16, main="",frame=FALSE, axes=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
axis(2, at=seq(0,700,100),labels=T, las=1, mgp=c(3,0.5,0), cex.axis=1.3)
axis(1, at=c(1980:2013),labels=T, , cex.axis=1.3)
points(x=c(1980:2013),fitted, type='l', col='black', lwd=2, lty=1)
points(x=c(1980:2013),Greece, type='l', col='grey38', lwd=2, lty=2)
points(x=c(1980:2013),FYROM, type='l', col='grey38', lwd=2, lty=3)
points(x=c(1980:2013),Bulgaria, type='l', col='grey38', lwd=1, lty=1)
#points(x=c(1980:2013),rowSums(EV), type='p', col='red', pch=4,lwd=2, lty=1)
legend(2013,700,lty=c(2,3,1), col="grey38", legend=c("Greece", "FYR of Macedonia", "Bulgaria"), lwd=c(2,2,1), bty='n',cex=1.2, xjust=1)



dev.off()
