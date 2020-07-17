##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS ON BALKAN PENINSULA
#
##########################################################################
# original model based on Lieury et al. 2015: Relative contribution of local demography and immigration in the recovery of a geographically-isolated population of the endangered Egyptian vulture. Biological Conservation 191: 349-356.

library(jagsUI)
library(tidyverse)
load("https:\\github.com\\steffenoppel\\vultures\\Balkan_EGVU_PVA_DATA.RData")
# load("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\EGVU_papers\\PVA_CaptiveRelease\\EGVU_IPM2020_output_v4_FINAL.RData")
# rm(list=setdiff(ls(), c("z.telemetry","INPUT")))
# head(INPUT)
# INPUT$nrep.terrvis<-NULL
# INPUT$rand.phi.offset<-NULL
# INPUT$J.fec.red<-NULL
# save.image("Balkan_EGVU_PVA_DATA.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPLANATION OF INPUT DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### the object 'INPUT' is a list with the following elements:
str(INPUT)
# y.terrvis: matrix of territory occupancy for 87 territories (rows) over 14 years (columns) with maximum number of adult birds observed per territory and year
# z.terrvis: matrix of initial values for territory occupancy for 87 territories (rows) over 14 years (columns) to initiate the estimation
# eff.terrvis: matrix of territory monitoring effort for 87 territories (rows) over 14 years (columns) scaled to the range -2 to +2
# f.obs: vector with one value per monitored territory specifying the first year of monitoring (1-13) 
# nsite.terrvis: number of observed territories for territory monitoring  (87)
# nprim.terrvis: number of observed primary occasions for territory monitoring (14)
# phase: vector of indicator (1 or 2) whether survival in a particular year was 'good' or 'bad'
# y.count: matrix with the number of adult breeding birds counted in each country (column, 4) in each year (row, 14)
# T.count: number of years with count data (14)
# countries: number of countries with count data (4)
# country.prop: vector of length 'countries' with proportion of the Balkan population in each country
# R.fec: vector of the number of monitored breeding pairs in each year for calculation of fecundity
# J.fec: vector of the number of fledglings produced by monitored breeding pairs in each year for calculation of fecundity
# y.telemetry: matrix of observed states (1-5) for 33 tracked juvenile birds (rows) over 118 months (columns):
# Observation states are as follows:
# 1 Tag ok, bird moving
# 2 Tag ok, bird not moving (dead)
# 3 Tag failed, bird observed alive
# 4 Dead bird recovered
# 5 No signal (=not seen)
# nind.telemetry: number of individual juveniles tracked with satellite telemetry (33)
# n.occasions.telemetry: number of months over which individuals were tracked (118)
# f.telemetry: vector with one value per tracked individual specifying the first month of tracking
# l.telemetry: vector with one value per tracked individual specifying the last month of tracking
# age.telemetry: matrix of individual age per month for 33 tracked juvenile birds (rows) over 119 months (columns)
# capt.telemetry: vector with one value per tracked individual specifying whether the bird was captive-bred (1) or wild (0)
# mig.telemetry: matrix specifying whether an individual migrated (1) or not (0) in a given month for 33 tracked juvenile birds (rows) over 119 months (columns)
# migprog.age: vector specifying at which age first autumn migration occurs in wild juveniles
# captprog.age: vector specifying at which age first autumn migration occurs in captive-released juveniles
# PROJECTION: number of years over which population is projected into the future (30)
# scen.capt.release: number of scenarios of captive-released birds (48)
# scen.imp.surv: number of scenarios of survival improvement in the wild (18)
# capt.release: matrix of the number of captive birds released every year (rows, 50 - only 30 are used if PROJECTION=30) for each scenario of captive release strategies (columns, 48)
# imp.surv: matrix of survival improvement (multiplier for estimated survival) for every year (rows, 30 - only 30 are used if PROJECTION=30) for each scenario of survival improvement (columns, 18)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INTEGRATED POPULATION MODEL FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("EGVU_IPM_Balkans.jags")
cat("
    model {
    #-------------------------------------------------
    # integrated population model for the Balkan Egyptian Vulture population
    # - age structured model with 6 age classes 
    # - age-specific probability to recruit at ages 4 or 5
    # - age-specific survival derived from tracking data for the first 3 years
    # - adult survival based on territory occupancy
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - FUTURE PROJECTION COMBINES SCENARIOS: baseline scenario (=do nothing, no chick removal)
    # - productivity supplemented by captive-bred juveniles (0-15 per year, for 10,20, or 30 years)
    # - improvement of survival 0-10% with time lag for full improvement to occur after 10 years
    #-------------------------------------------------
    
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    # Priors and constraints FOR FECUNDITY
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
    sigma.obs.count[s] ~ dunif(0.1,100)	#Prior for SD of observation process (variation in detectability)
    tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    # Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
    #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind.telemetry){
    for (t in f.telemetry[i]:(n.occasions.telemetry)){
    logit(phi.telemetry[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
    b.phi.capt*(capt.telemetry[i]) +     ### survival dependent on captive-release (captive-raised or other)
    b.phi.mig*(mig.telemetry[i,t]) +     ### survival dependent on first migration across the sea
    b.phi.age*(age.telemetry[i,t])     ### survival dependent on age (juvenile or other)
    } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi.telemetry ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi.telemetry/(1 - mean.phi.telemetry))    # logit transformed survival intercept
    
    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for changes in survival probability with age on logit scale
    b.phi.capt ~ dnorm(0, 0.001)          # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)               # Prior for  COST OF MIGRATION on survival probability on logit scale
    
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind.telemetry){
    for (t in f.telemetry[i]:(n.occasions.telemetry)){
    logit(p.obs.telemetry[i,t]) <- base.obs.telemetry
    logit(tag.fail.telemetry[i,t]) <- base.fail.telemetry
    logit(p.found.dead.telemetry[i,t]) <- base.recover.telemetry
    } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs.telemetry ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover.telemetry ~ dnorm(0, 0.001)            # Prior for intercept of recovery probability on logit scale
    
    
    # Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    ## Priors for detection probability
    lmu.p.terrvis <- log(mean.p.terrvis/(1 - mean.p.terrvis))    # logit transformed survival intercept
    mean.p.terrvis ~ dunif(0.5, 1)			
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    # Priors for survival
    for (nypterr in 1:2){   ## only 2 survival 'phases' - good or bad years
    mean.phi.terrvis[nypterr] ~ dunif(0.5, 1)   # informative prior for annual survival probability
    mean.rec.terrvis[nypterr] ~ dunif(0, 1)     # informative prior for annual recruitment probability after territory abandonment
    } 
    
    
    ### RANDOM OBSERVATION AND SURVIVAL EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
    rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    
    
    # Priors for population process and immigration
    
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(50, 100)   ##changed from JUV
    N1[1] ~ dunif(10, 35)
    N2[1] ~ dunif(5, 20)
    N3[1] ~ dunif(3, 15)
    N4[1] ~ dunif(0, 10)
    N5[1] ~ dunif(0, 5)
    N6[1] ~ dunif(200, 250)
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T.count){
    
    nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
    N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
    N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
    N4[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N4[tt-1]))                                                ### number of 5-year old survivors
    N6[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds
    
    } # tt
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (tlc in 1:T.count){
    Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
    for (s in 1:countries){			### start loop over every country
    Nterr.country[tlc,s] ~ dbin(country.prop[s],round(Nterr[tlc]))
    y.count[tlc,s] ~ dnorm(Nterr.country[tlc,s], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
    }
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    
    J.fec[tlc] ~ dpois(rho.fec[tlc])
    rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for juvenile survival from telemetry
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind.telemetry){
    
    for (t in f.telemetry[i]:(n.occasions.telemetry-1)){
    
    # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
    
    ps[1,i,t,1]<-1    ## dead birds stay dead
    ps[1,i,t,2]<-0
    ps[1,i,t,3]<-0
    
    ps[2,i,t,1]<-(1-phi.telemetry[i,t])
    ps[2,i,t,2]<-phi.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
    ps[2,i,t,3]<-phi.telemetry[i,t] * tag.fail.telemetry[i,t]
    
    ps[3,i,t,1]<-(1-phi.telemetry[i,t])
    ps[3,i,t,2]<-0
    ps[3,i,t,3]<-phi.telemetry[i,t]
    
    # Define probabilities of O(t) [last dim] given S(t)  [first dim]
    
    po[1,i,t,1]<-0
    po[1,i,t,2]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t]) * (1-p.found.dead.telemetry[i,t])
    po[1,i,t,3]<-0
    po[1,i,t,4]<-p.found.dead.telemetry[i,t]
    po[1,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t] * (1-p.found.dead.telemetry[i,t])
    
    po[2,i,t,1]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
    po[2,i,t,2]<-0
    po[2,i,t,3]<-0
    po[2,i,t,4]<-0
    po[2,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t]
    
    po[3,i,t,1]<-0
    po[3,i,t,2]<-0
    po[3,i,t,3]<-0
    po[3,i,t,4]<-0
    po[3,i,t,5]<-1
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind.telemetry){
    # Define latent state at first capture
    z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
    for (t in (f.telemetry[i]+1):n.occasions.telemetry){
    # State process: draw S(t) given S(t-1)
    z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
    } #t
    } #i
    
    
    # -------------------------------------------------        
    # 2.5. Likelihood for adult survival from territory monitoring
    # -------------------------------------------------
    
    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
    x.terrvis[ilterr,f.obsvis[ilterr]] <- 2 #firstobs[ilterr]
    
    for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
    z.terrvis[ilterr,klterr] ~ dbin(mean.phi.terrvis[phase[klterr]],x.terrvis[ilterr,klterr-1])
    mu.x[ilterr,klterr] ~ dbin(mean.rec.terrvis[phase[klterr]],(2-z.terrvis[ilterr,klterr]))
    x.terrvis[ilterr,klterr] <- mu.x[ilterr,klterr] + z.terrvis[ilterr,klterr]
    
    } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    
    for (iobs in 1:nsite.terrvis){
    for (kobs in f.obsvis[iobs]:nprim.terrvis){
    linpred.terrvis[iobs,kobs] <- lmu.p.terrvis+beta.obs.eff.terrvis*eff.terrvis[iobs,kobs] + rand.obs.terrvis[kobs]
    p.terrvis[iobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,kobs]))
    y.terrvis[iobs,kobs] ~ dbin(p.terrvis[iobs,kobs], x.terrvis[iobs,kobs])
    
    } 						 # close kobs loop over year
    } 							# close iobs loop over sites
    
    
    
    
    # -------------------------------------------------        
    # 3. DERIVED PARAMETERS
    # -------------------------------------------------
    
    ### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES
    
    ## for WILD BIRDS
    for (ageprog in 1:36){
    logit(phi.wild.telemetry[ageprog]) <- lp.mean.telemetry +      ### intercept for mean survival 
    b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
    b.phi.age*ageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.juv.telemetry<-phi.wild.telemetry[1]*phi.wild.telemetry[2]*phi.wild.telemetry[3]*phi.wild.telemetry[4]*phi.wild.telemetry[5]*phi.wild.telemetry[6]*phi.wild.telemetry[7]*phi.wild.telemetry[8]*phi.wild.telemetry[9]*phi.wild.telemetry[10]*phi.wild.telemetry[11]*phi.wild.telemetry[12]					### multiply monthly survival from age 1-12
    ann.phi.sec.telemetry<-phi.wild.telemetry[13]*phi.wild.telemetry[14]*phi.wild.telemetry[15]*phi.wild.telemetry[16]*phi.wild.telemetry[17]*phi.wild.telemetry[18]*phi.wild.telemetry[19]*phi.wild.telemetry[20]*phi.wild.telemetry[21]*phi.wild.telemetry[22]*phi.wild.telemetry[23]*phi.wild.telemetry[24]					### multiply monthly survival from age 13-24
    ann.phi.third.telemetry<-phi.wild.telemetry[25]*phi.wild.telemetry[26]*phi.wild.telemetry[27]*phi.wild.telemetry[28]*phi.wild.telemetry[29]*phi.wild.telemetry[30]*phi.wild.telemetry[31]*phi.wild.telemetry[32]*phi.wild.telemetry[33]*phi.wild.telemetry[34]*phi.wild.telemetry[35]*phi.wild.telemetry[36]					### multiply monthly survival from age 25-36
    
    ## for CAPTIVE-REARED DELAYED RELEASE BIRDS
    for (captageprog in 1:36){
    logit(phi.capt.telemetry[captageprog]) <- lp.mean.telemetry +      ### intercept for mean survival
    b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
    b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
    b.phi.age*captageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.capt.rel.first.year<-phi.capt.telemetry[8]*phi.capt.telemetry[9]*phi.capt.telemetry[10]*phi.capt.telemetry[11]*phi.capt.telemetry[12]*phi.capt.telemetry[13]*phi.capt.telemetry[14]*phi.capt.telemetry[15]*phi.capt.telemetry[16]*phi.capt.telemetry[17]*phi.capt.telemetry[18]*phi.capt.telemetry[19]*phi.capt.telemetry[20]*phi.capt.telemetry[21]*phi.capt.telemetry[22]*phi.capt.telemetry[23]*phi.capt.telemetry[24]					### first year for delayed-release bird is longer, from month 8 to 24 
    
    ### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    
    
    # Annual population growth rate
    for (ipop in 1:(T.count-1)){
    lambda.t[ipop] <- Nterr[ipop+1] / max(Nterr[ipop],1)
    loglambda.t[ipop]<-log(lambda.t[ipop])## for calculating geometric mean of overall population growth rate
    }
    
    #### OVERALL POPULATION GROWTH RATE  #########
    mean.lambda<-exp((1/(T.count-1))*sum(loglambda.t[1:(T.count-1)]))   # Geometric mean
    
    
    
    
    # -------------------------------------------------        
    # 4. PREDICTION INTO THE FUTURE
    # -------------------------------------------------
    ### INCLUDE SCENARIOS FOR N CAPTIVE BIRDS AND SURVIVAL IMPROVEMENT
    
    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){
    
    ### FUTURE FECUNDITY IS BASED ON mu.fec - removed the chick removal scenario
    
    for (fut in 1:PROJECTION){
    fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
    }
    
    
    
    # SPECIFY IMPROVEMENT OF SURVIVAL
    for (is in 1:scen.imp.surv){ 
    
    
    ## POPULATION PROCESS
    ### DOES NOT WORK WITH SAME ARRAY OF POP TRAJECTORY AS 3 dimensions required
    ## need to copy previous array elements
    nestlings.f[ncr,is,1] <- round(nestlings[T.count])   ##JUV[T.count]
    N1.f[ncr,is,1] <- N1[T.count]
    N2.f[ncr,is,1] <- N2[T.count]
    N3.f[ncr,is,1] <- N3[T.count]
    N4.f[ncr,is,1] <- N4[T.count]
    N5.f[ncr,is,1] <- N5[T.count]
    N6.f[ncr,is,1] <- N6[T.count]
    Nterr.f[ncr,is,1] <- Nterr[T.count]
    N2wild.f[ncr,is,1] <- 0
    N2released.f[ncr,is,1] <- 0
    
    for (fut in 2:PROJECTION){
    
    fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean.phi.terrvis[2],1) ### invalid parent error if survival>1
    
    ### probabilistic formulation
    nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as fecundity times number of territorial pairs
    N1.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### number of 1-year old survivors
    N2wild.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))                ### number of 2-year old wild survivors
    N2released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.capt.rel.first.year),1),round(capt.release[fut,ncr]))             ### number of 2-year old captive-released survivors 
    N2.f[ncr,is,fut] <-  N2wild.f[ncr,is,fut] + N2released.f[ncr,is,fut]       ### sum of the N2 cohort derived from wild and captive-bred juveniles 
    N3.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.third.telemetry),1),round(N2.f[ncr,is,fut-1]))                                                    ### number of 3-year old survivors
    N4.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round(N3.f[ncr,is,fut-1]))                                                       ### number of 4-year old survivors
    N5.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round(N4.f[ncr,is,fut-1]))                                                       ### number of 5-year old survivors
    N6.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round((N5.f[ncr,is,fut-1]+N6.f[ncr,is,fut-1])))                                   ### number of 6-year or older (adult) birds
    Nterr.f[ncr,is,fut] <- round((N4.f[ncr,is,fut] * 0.024) + (N5.f[ncr,is,fut] * 0.124) + (N6.f[ncr,is,fut]))
    
    
    } # fut
    
    for (fut2 in 1:(PROJECTION-1)){
    lambda.t.f[ncr,is,fut2] <- Nterr.f[ncr,is,fut2+1] / max(Nterr.f[ncr,is,fut2],1)
    loglambda.t.f[ncr,is,fut2]<-log(lambda.t.f[ncr,is,fut2])## for calculating geometric mean of overall population growth rate
    } # fut2
    
    
    #### FUTURE POPULATION GROWTH RATE  #########
    fut.lambda[ncr,is]<-exp((1/(PROJECTION-1))*sum(loglambda.t.f[ncr,is,1:(PROJECTION-1)]))   # Geometric mean
    
    } # end scenario of imp.surv
    
    } # end scenario of n capt released
    
    } # END of the JAGS model
    ",fill = TRUE)
sink()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP MCMC PARAMETERS AND RUN MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("mu.fec","lambda.t","b.phi.age","b.phi.capt","b.phi.mig","ann.phi.capt.rel.first.year",
           "ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry", "ann.phi.terrvis",     #"breed.prop4","breed.prop5",
           "mean.phi.terrvis","mean.lambda","fut.lambda","Nterr", "Nterr.f")


# MCMC settings
nc <- 3
nt <- 4
ni <- 50000
nb <- 10000

# specify initial values
initIPM <- function(){list(mean.p.terrvis=runif(1,0.5,1),
                           mean.phi.terrvis=runif(2,0.75, 1),
                           sigma.obs.count=runif(4,0.1,100),
                           mu.fec = runif(1,0,1),
                           z.telemetry = z.telemetry,                 ### matrix with initial values for the true states of tracked individuals
                           mean.phi.telemetry = runif(1, 0.9, 0.999),
                           base.obs.telemetry = rnorm(1,0, 0.001),
                           base.fail.telemetry = rnorm(1,0, 0.001),
                           base.recover.telemetry = rnorm(1,0, 0.001))}   


### RUN THE MODEL ###

EGVU.IPM <- jags(data=INPUT,
                       inits=initIPM,
                       parameters.to.save=paraIPM,
                       model.file="C:/STEFFEN/RSPB/Bulgaria/Analysis/PopulationModel/vultures/EGVU_IPM_Balkans.jags",    
                       n.chains=nc, n.thin=nt, n.burnin=nb, n.iter=ni,parallel=T)







  