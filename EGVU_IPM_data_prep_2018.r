# IPM for French Egyptian vulture populations
# specify model in bugs language
sink("NeoIPMi.bug")
cat("
model {
#-------------------------------------------------
# integrated population model base of the south-eastern population
# - age structured model with 6 age classes: 
# - status structured model with territorials and non-territorials
# - age-specific probability to recruit with 4 age classes (1-2 ; 3-4;5; 6+)
# - recovery
# - probability of ring loss
# - different resighting probability of birds having lost at least one ring
# - check if recapture probability of floaters do not increase with time
# - post breeding census, female-based assuming equal sex ratio & survival
# - survival parameters depend on supplementary feeding and NDVI
# - constant immigration rate and annual estimate of immigrant number
#-------------------------------------------------

# 1 Constrains and Priors: fixed age effect # see Kery and Schaub 7.7 #***********************************
# age-independent parameters 
   		 # Survival and resighting 
   for (tt in 1:(nyears-1)){
stdP[tt] <- (Placette[tt]-mean(Placette[]))/sd(Placette[]) 
# standardized number of vulture restaurant
stdNDVI[tt] <- (NDVI[tt]-mean(NDVI[]))/sd(NDVI[]) # standardized NDVI value                  
stdT[tt+1] <- (Time[tt+1]-mean(Time[]))/sd(Time[]) # standardized Time
                    
logit(sT[tt]) <- phiT.lim[tt]  # territorial survival + temporal variation
phiT.lim[tt] <- min(999, max(-999, phiT[tt])) #stabilize the logit
phiT[tt] <- mu.phiT + alpha*stdP[tt]  
# NDVI and conservation influence variation of territorial survival

pT[tt+1] <- mean.pT  # constant resighting probability of territorial
pLT[tt+1] <- mean.pLT  # different resighting probability after ring loss
      r[tt] <- mean.R 
      L[tt] <-mean.L
            } #t
mu.phiT <- log(mean.phiT / (1-mean.phiT))
    
    		# fecundity   
    for (tt in 1:nyears){  
log(fecSE[tt])<- mu.fec + eta*NtotSE[tt]  } #t
# Density feedback on Fecundity = Nb of fledging / estimated Nb of territorial 

    		# immigration
    for (tt in 1:(nyears-1)){ 
stdTd[tt] <- (Time[tt]-mean(Time[]))/sd(Time[])  
log(im[tt])<- mu.im + epsim[tt]     # Time trend of immigration rate
   epsim[tt] ~ dnorm(0,tauepsim)I(-10,10) } #t

# age-dependent parameters 
for (i in 1:nind){
   for (tt in f[i]:(nyears-1)){
logit(sF[i,tt]) <- phiF.lim[x[i,tt],tt]  
# age-specific floater survival + temporal variation
logit(pF[i,tt+1]) <- pF.lim[w[i,tt+1],tt+1]    
# age-specific floater resighting probability + temporal variation    
logit(pLF[i,tt+1]) <- pLF.lim[w[i,tt+1],tt+1]  
# age-specific floater resighting after tag loss + temporal variation    
a[i,tt] <- set[x[i,tt],tt]   } #t    
# age-specific recruitment probability between its birthday x-1 and x Kery&Schaub7
   } #i 
                 
  # for survival parameters and age at first recruitment
for (u in 1:4){
     for (tt in 1:(nyears-1)){
phiF.lim[u,tt] <- min(999, max(-999, phiF[u,tt])) #stabilize the logit
phiF[u,tt] <- mu.phiF[u] + alpha*stdP[tt] + beta*stdNDVI[tt]   
# Same restaurant effect between age class + NDVI variation
set[u,tt] <- mean.set[u]  
         } #t         
mu.phiF[u] <- log(mean.phiF[u] / (1-mean.phiF[u]))        
      } #u

  # for resighting probability
for (u in 1:3){
     for (tt in 1:(nyears-1)){
pF.lim[u,tt+1] <- min(999, max(-999, prF[u,tt+1])) #stabilize the logit
prF[u,tt+1] <- mu.pF[u] + kapa*stdT[tt+1]  # Trend of floater resighting
     
pLF.lim[u,tt+1] <- min(999, max(-999, prLF[u,tt+1])) #stabilize the logit
prLF[u,tt+1] <- mu.pLF[u] + kapa*stdT[tt+1]    
         } #t
      mu.pF[u] <- log(mean.pF[u] / (1-mean.pF[u]))        
      mu.pLF[u] <- log(mean.pLF[u] / (1-mean.pLF[u]))        
      } #u
  #************************************************
  # 2. Define the priors for the parameters     : non informative priors
  #************************************************
# Non informative priors
      mean.phiT ~ dunif(0, 1)         # Priors on survival of territorials
      mean.pT ~ dunif(0, 1)           # Priors on resighting of territorial
      mean.R ~ dunif(0, 1)            # Priors on recovery probability  
      mean.L ~ dunif(0, 1)            # Annual probability of ring loss  
      mean.pLT ~ dunif(0, 1)          # Lower resighting after ring-loss
      mean.phiF[1] ~ dunif(0, 1)      # Priors for age-spec. survival 
      mean.phiF[2] ~ dunif(0, 1)         
      mean.phiF[3] ~ dunif(0, 1)         
      mean.phiF[4] ~ dunif(0, 1)         
      mean.set[1] <- 0 # Priors for age-spec. recruitment, 0 for juveniles
      mean.set[2] ~ dunif(0, 1)      # Priors for age-spec. recruitment 
      mean.set[3] ~ dunif(0, 1)        
      mean.set[4] ~ dunif(0, 1)         
      mean.pF[1] ~ dunif(0, 1)       # Priors for age-spec. resighting rate 
      mean.pF[2] ~ dunif(0, 1)         
      mean.pF[3] ~ dunif(0, 1)         
      mean.pLF[1] <- mean.pF[1]   # Priors for age-spec. resighting rate 
      mean.pLF[2] <- mean.pF[2]          
      mean.pLF[3] ~ dunif(0, 1)   
      mu.fec ~ dunif(0, 2)           # Priors on fecundity
      mu.im ~ dunif(-5, 1)            # Priors on immigration rate

# Informative priors for parameters from Grande et al 2006    
# mean.phiT ~ dbeta(5, 1)  # territorial survival 0.83+/-0.02 SE
# mean.pT ~ dbeta(15, 1)   # territorial resighting 0.94+/-0.02
# mean.phiF[1] ~ dbeta(3, 1)    # age-spec. survival 0.73 +/- 0.02 
# mean.phiF[2] ~ dbeta(4, 1)    # 0.78 +/- 0.03 SE
# mean.phiF[3] ~ dbeta(4.8, 2.8)# 0.60 +/- 0.05 SE
# mean.phiF[4] ~ dbeta(3, 1)    # 0.75 +/- 0.02 SE
# mean.pF[1] ~ dbeta(1, 15)  # age-spec. resighting rate 0.05 +/- 0.02 SE
# mean.pF[2] ~ dbeta(1, 5)   #  0.15 +/- 0.02 SE
# mean.pF[3] ~ dbeta(1, 1)   #  0.55 +/- 0.05 SE
            
        # priors of covariable linear coefficient centered at 0 
            alpha ~ dnorm(0,0.01)I(-5,5) 
            eta  ~ dnorm(0,0.01)I(-5,5) 
            beta  ~ dnorm(0,0.01)I(-5,5) 
            kapa ~ dnorm(0,0.01)I(-5,5) 
            sigepsim ~ dunif(0,10) ; tauepsim <- pow(sigepsim,-2)    
          # Priors for initial population sizes
          #*****************************************
          njSE[1] ~ dnorm(7, 5)I(0,)  # number of 1years old
          n1SE[1] ~ dnorm(6, 5)I(0,)  # number of 1years old
          n2SE[1] ~ dnorm(5, 5)I(0,)  # number of 2years old
          n3SE[1] ~ dnorm(4, 5)I(0,)  # number of 3years old
          n4SE[1] ~ dnorm(3, 5)I(0,)  # number of 4years old   
          n5SE[1] ~ dnorm(5, 5)I(0,)  # number of 5years old   
          
          NtotSE[1] ~ dnorm(18, 1)I(0,)    # number of territorial
                     
  # 3. Derived parameters
  #********************************************
         # Population growth rate r[t], Mean population growth rate
for(tt in 1:(nyears-1)) {
   lambdaSE[tt] <- NtotSE[tt+1]/NtotSE[tt]
   loglaSE[tt] <- log(lambdaSE[tt])} #tt
 mean.lambdaSE <- exp((1/(nyears-1))*sum(loglaSE[1:(nyears-1)]))
          
         # Annual survival parameters
for (tt in 1:(nyears-1)){
      surT[tt] <- exp(phiT[tt])/(1+exp(phiT[tt])) } #t # Territorial survival      
for (tt in 1:(nyears-1)){        
   for (u in 1:4){
      surF[u,tt] <- exp(phiF[u,tt])/(1+exp(phiF[u,tt])) } #u } #t# Floaters survival   

         # Annual resighting parameters
for (tt in 1:(nyears-1)){        
      recF[3,tt+1] <- exp(prF[3,tt+1])/(1+exp(prF[3,tt+1]))     } #t
      
   # Unstandardized coefficients
         NDVIF <- beta/sd(NDVI[])
         Food <- alpha/sd(Placette[])
         pTime <- kapa/sd(Time[])
 mean.im <- exp(mu.im)                   # mean immigration rate 
  
  # 4. Likelihoods of the single data sets
  #-------------------------------------------------
# 4.1. Likelihood for population count data (state-space model)
   # 4.1.1 System process: female based matrix model 
                                                     
   for (ti in 2:nyears){
        # Floaters compartment     
      n1SE[ti] ~ dbin(surF[1,ti-1],njSE[ti-1])  # number of 1yr old
      n2SE[ti] ~ dbin(surF[1,ti-1],n1SE[ti-1])  # number of 2yr old
 
  sn3[ti-1]<-surF[2,ti-1]*(1-set[2,ti-1])# product of survival and recruitment
  sn4[ti-1]<-surF[2,ti-1]*(1-set[2,ti-1])
  sn5[ti-1]<-surF[3,ti-1]*(1-set[3,ti-1])
  sn6[ti-1]<-surF[4,ti-1]*(1-set[4,ti-1])  
  
      n3SE[ti] ~ dbin(sn3[ti-1],n2SE[ti-1])     # number of 3yr old
      n4SE[ti] ~ dbin(sn4[ti-1],n3SE[ti-1])     # number of 4yr old   
      n5SEa[ti] ~ dbin(sn5[ti-1],n4SE[ti-1])    # number of 5yr old 
      n5SEb[ti] ~ dbin(sn6[ti-1],n5SE[ti-1])   # number of 6yr old +       
      n5SE[ti] <- n5SEa[ti] + n5SEb[ti]       # total number of 5yr old  +

        # Territorials compartment
  st3[ti-1]<-surF[2,ti-1]*set[2,ti-1]       
  st4[ti-1]<-surF[2,ti-1]*set[2,ti-1]
  st5[ti-1]<-surF[3,ti-1]*set[3,ti-1]
  st6[ti-1]<-surF[4,ti-1]*set[4,ti-1]
      t3SE[ti] ~ dbin(st3[ti-1],n2SE[ti-1])  # new territorial of 3yro
      t4SE[ti] ~ dbin(st4[ti-1],n3SE[ti-1])  # new territorial of 4yr old   
      t5SEa[ti] ~ dbin(st5[ti-1],n4SE[ti-1])  # new territorial of 5yr old   
      t5SEb[ti] ~ dbin(st6[ti-1],n5SE[ti-1])   # new territorial of 6yr old 
     floSE[ti] <- n1SE[ti] + n2SE[ti] + n3SE[ti] + n4SE[ti] + n5SE[ti] 
# Annual number of floaters
      newSE[ti] <- t3SE[ti] + t4SE[ti] + t5SEa[ti] + t5SEb[ti] 
# Annual new territorial recruits
totSE[ti] ~ dbin(sT[ti-1],NtotSE[ti-1]) # territorials alive from the last year
      imSE[ti] <- im[ti-1] * NtotSE[ti-1]
IMSE[ti] ~ dpois(imSE[ti])  # Annual number of immigrants
      
NtotSE[ti]<- t3SE[ti] +t4SE[ti] +t5SEa[ti] +t5SEb[ti] +IMSE[ti] +totSE[ti]        
# total population of  territorials
      
        # Reproduction process
      juvSE[ti] <- fecSE[ti] *0.5 * NtotSE[ti]
      njSE[ti] ~ dpois(juvSE[ti])          # number of fledglings          
       # Observation process      
  NTSE[ti] ~ dpois(NtotSE[ti])       } #t # Poisson observation error of census count      

# 4.2 Likelihood for capture-recapture data: multi state model 
# -------------------------------------------------
# Parameters:
# sF: survival probability of floaters
# sT: survival probability of territorials
# a(i): probability to start breeding at age i 
# pF: recapture probability of floaters
# pT: recapture probability of territorials  
# r: probability to recover a recently dead animals
# L: annual probability to lose one colored plastic ring
# -------------------------------------------------
# States (S) = biological meaningful state from identifying processes
# 1 fledgling with all coloured rings
# 2 floater with all coloured rings
# 3 territorial with all coloured rings
# 4 recently dead and recovered
# 5 floater with one ring lost
# 6 territorial with one ring lost
# 7 floater with two ring lost
# 8 territorial with two rings lost
# 9 dead or recently dead and not recovered
# 10 floater with all ring lost
# 11 territorial with all ring lost

# Observations (O) = observation process, data we have
# 1 seen as juvenile = ringed
# 2 seen as floater
# 3 seen as territorial
# 4 recovered dead
# 5 identify as floater despite one ring lost
# 6 identify as territorial despite one ring lost
# 7 identify as floater despite two ring lost
# 8 identify as territorial despite two rings lost
# 9 not seen or recovered
# 10 identify as floater despite all ring lost
# 11 identify as territorial despite all ring lost
# -------------------------------------------------
# define state-transition and observation matrices   
for (i in 1:nind){  # for each bird, i.e. we translate the likelihood of each individuals life history : such individual effect in necessary for age-variation of survival

	# Transition probability between t+1 and t
   for (t in f[i]:(nyears-1)){                  
# the first occasion indicate the ringing year, i.e. the start of age
# Loop over time: each t is one year more for the individual
# Define probabilities of state S(t+1) given S(t) : Transition matrix between each state
     # First index = states at time t, last index = states at time t+1
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- sF[i,t] * (1-a[i,t]) * (1-L[t])     
# Transition probability during year t while individuals were aged at x
      ps[1,i,t,3] <- sF[i,t] * a[i,t] * (1-L[t]) 
      ps[1,i,t,4] <- (1-sF[i,t])*r[t]
      ps[1,i,t,5] <- sF[i,t] * (1-a[i,t]) * L[t] *(1-L[t]) 
      ps[1,i,t,6] <- sF[i,t] * a[i,t] * L[t]  *(1-L[t])
      ps[1,i,t,7] <- sF[i,t] * (1-a[i,t]) * L[t] * L[t] *(1-L[t])
      ps[1,i,t,8] <- sF[i,t] * a[i,t] * L[t] * L[t] *(1-L[t])
      ps[1,i,t,9] <- (1-sF[i,t])*(1-r[t])
      ps[1,i,t,10] <- sF[i,t] * (1-a[i,t]) * L[t] * L[t] *L[t]
      ps[1,i,t,11] <- sF[i,t] * a[i,t] * L[t] * L[t] *L[t]
                        
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- sF[i,t] * (1-a[i,t]) * (1-L[t])  
      ps[2,i,t,3] <- sF[i,t] * a[i,t]  * (1-L[t])  
      ps[2,i,t,4] <- (1-sF[i,t])*r[t]
      ps[2,i,t,5] <- sF[i,t] * (1-a[i,t]) * L[t] *(1-L[t]) 
      ps[2,i,t,6] <- sF[i,t] * a[i,t]  * L[t] *(1-L[t]) 
      ps[2,i,t,7] <- sF[i,t] * (1-a[i,t]) * L[t] * L[t] *(1-L[t])  
      ps[2,i,t,8] <- sF[i,t] * a[i,t]  * L[t] * L[t] *(1-L[t])  
      ps[2,i,t,9] <- (1-sF[i,t])*(1-r[t])
      ps[2,i,t,10] <- sF[i,t] * (1-a[i,t]) * L[t] * L[t] * L[t]  
      ps[2,i,t,11] <- sF[i,t] * a[i,t]  * L[t] * L[t] * L[t]  

      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- sT[t] * (1-L[t])        # not have lost any ring
      ps[3,i,t,4] <- (1-sT[t])*r[t]
      ps[3,i,t,5] <- 0
      ps[3,i,t,6] <- sT[t] * L[t]*(1-L[t])   # have lost one ring only
      ps[3,i,t,7] <- 0
      ps[3,i,t,8] <- sT[t] * L[t] * L[t] *(1-L[t])
      ps[3,i,t,9] <- (1-sT[t])*(1-r[t])      # have lost two ring only
      ps[3,i,t,10] <- 0
      ps[3,i,t,11] <- sT[t] * L[t] * L[t] * L[t] # have lost all ring

      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- 0
      ps[4,i,t,6] <- 0
      ps[4,i,t,7] <- 0
      ps[4,i,t,8] <- 0
      ps[4,i,t,9] <- 1
      ps[4,i,t,10] <- 0
      ps[4,i,t,11] <- 0
      
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0 
      ps[5,i,t,3] <- 0 
      ps[5,i,t,4] <- (1-sF[i,t])*r[t]
      ps[5,i,t,5] <- sF[i,t] * (1-a[i,t]) * (1-L[t]) 
      ps[5,i,t,6] <- sF[i,t] * a[i,t] * (1-L[t])
      ps[5,i,t,7] <- sF[i,t] * (1-a[i,t]) * L[t] * (1-L[t])  
      ps[5,i,t,8] <- sF[i,t] * a[i,t]  * L[t] * (1-L[t]) 
      ps[5,i,t,9] <- (1-sF[i,t])*(1-r[t])
      ps[5,i,t,10] <- sF[i,t] * (1-a[i,t]) * L[t] * L[t]   
      ps[5,i,t,11] <- sF[i,t] * a[i,t]  * L[t] * L[t]  
      
      ps[6,i,t,1] <- 0
      ps[6,i,t,2] <- 0
      ps[6,i,t,3] <- 0
      ps[6,i,t,4] <- (1-sT[t])*r[t]
      ps[6,i,t,5] <- 0
      ps[6,i,t,6] <- sT[t] * (1-L[t])
      ps[6,i,t,7] <- 0
      ps[6,i,t,8] <- sT[t] * L[t] * (1-L[t])
      ps[6,i,t,9] <- (1-sT[t])*(1-r[t])
      ps[6,i,t,10] <- 0
      ps[6,i,t,11] <- sT[t] * L[t] * L[t] 

      ps[7,i,t,1] <- 0
      ps[7,i,t,2] <- 0 
      ps[7,i,t,3] <- 0 
      ps[7,i,t,4] <- (1-sF[i,t])*r[t]
      ps[7,i,t,5] <- 0
      ps[7,i,t,6] <- 0
      ps[7,i,t,7] <- sF[i,t] * (1-a[i,t]) * (1-L[t])   
      ps[7,i,t,8] <- sF[i,t] * a[i,t] * (1-L[t])  
      ps[7,i,t,9] <- (1-sF[i,t])*(1-r[t])
      ps[7,i,t,10] <- sF[i,t] * (1-a[i,t]) * L[t]   
      ps[7,i,t,11] <- sF[i,t] * a[i,t] * L[t]  
       
      ps[8,i,t,1] <- 0
      ps[8,i,t,2] <- 0
      ps[8,i,t,3] <- 0
      ps[8,i,t,4] <- (1-sT[t])*r[t]
      ps[8,i,t,5] <- 0
      ps[8,i,t,6] <- 0
      ps[8,i,t,7] <- 0
      ps[8,i,t,8] <- sT[t] * (1-L[t])
      ps[8,i,t,9] <- (1-sT[t])*(1-r[t])
      ps[8,i,t,10] <- 0
      ps[8,i,t,11] <- sT[t] * L[t] 

      ps[9,i,t,1] <- 0
      ps[9,i,t,2] <- 0
      ps[9,i,t,3] <- 0
      ps[9,i,t,4] <- 0
      ps[9,i,t,5] <- 0
      ps[9,i,t,6] <- 0
      ps[9,i,t,7] <- 0
      ps[9,i,t,8] <- 0
      ps[9,i,t,9] <- 1
      ps[9,i,t,10] <- 0
      ps[9,i,t,11] <- 0

      ps[10,i,t,1] <- 0
      ps[10,i,t,2] <- 0 
      ps[10,i,t,3] <- 0 
      ps[10,i,t,4] <- (1-sF[i,t])*r[t]
      ps[10,i,t,5] <- 0
      ps[10,i,t,6] <- 0
      ps[10,i,t,7] <- 0  
      ps[10,i,t,8] <- 0 
      ps[10,i,t,9] <- (1-sF[i,t])*(1-r[t])
      ps[10,i,t,10] <- sF[i,t] * (1-a[i,t])   
      ps[10,i,t,11] <- sF[i,t] * a[i,t]  
      ps[11,i,t,1] <- 0
      ps[11,i,t,2] <- 0
      ps[11,i,t,3] <- 0
      ps[11,i,t,4] <- (1-sT[t])*r[t]
      ps[11,i,t,5] <- 0
      ps[11,i,t,6] <- 0
      ps[11,i,t,7] <- 0
      ps[11,i,t,8] <- 0
      ps[11,i,t,9] <- (1-sT[t])*(1-r[t])
      ps[11,i,t,10] <- 0
      ps[11,i,t,11] <- sT[t]  

      # Define probabilities of O(t+1) given S(t+1)  : Transition matrix between state and observation, i.e. matrix of recapture probability 
# First index = states at time t+1, last index = observations at time t+1
# Observation process occurs at the beginning of the next year!
      po[1,i,t+1,1] <- 0
      po[1,i,t+1,2] <- 0
      po[1,i,t+1,3] <- 0
      po[1,i,t+1,4] <- 0
      po[1,i,t+1,5] <- 0
      po[1,i,t+1,6] <- 0
      po[1,i,t+1,7] <- 0
      po[1,i,t+1,8] <- 0
      po[1,i,t+1,9] <- 1  # as in Kery and Schaub’s Book; all juvenile are not mark and we are not interested by this observation transition t+1o estimate survival
      po[1,i,t+1,10] <- 0
      po[1,i,t+1,11] <- 0

# Observation probability at the beginning of the t+1 year when birds was aged x+1

 
      po[2,i,t+1,1] <- 0
      po[2,i,t+1,2] <- pF[i,t+1]      
      po[2,i,t+1,3] <- 0
      po[2,i,t+1,4] <- 0
      po[2,i,t+1,5] <- 0
      po[2,i,t+1,6] <- 0
      po[2,i,t+1,7] <- 0
      po[2,i,t+1,8] <- 0
      po[2,i,t+1,9] <- 1-pF[i,t+1]
      po[2,i,t+1,10] <- 0
      po[2,i,t+1,11] <- 0
      
po[3,i,t+1,1] <- 0
      po[3,i,t+1,2] <- 0
      po[3,i,t+1,3] <- pT[t+1]
      po[3,i,t+1,4] <- 0
      po[3,i,t+1,5] <- 0
      po[3,i,t+1,6] <- 0
      po[3,i,t+1,7] <- 0
      po[3,i,t+1,8] <- 0
      po[3,i,t+1,9] <-  1-pT[t+1]
      po[3,i,t+1,10] <- 0
      po[3,i,t+1,11] <- 0

      po[4,i,t+1,1] <- 0
      po[4,i,t+1,2] <- 0
      po[4,i,t+1,3] <- 0
      po[4,i,t+1,4] <- 1
      po[4,i,t+1,5] <- 0
      po[4,i,t+1,6] <- 0
      po[4,i,t+1,7] <- 0
      po[4,i,t+1,8] <- 0
      po[4,i,t+1,9] <- 0
      po[4,i,t+1,10] <- 0
      po[4,i,t+1,11] <- 0
      po[5,i,t+1,1] <- 0
      po[5,i,t+1,2] <- 0
      po[5,i,t+1,3] <- 0
      po[5,i,t+1,4] <- 0
      po[5,i,t+1,5] <- pLF[i,t+1]      
      po[5,i,t+1,6] <- 0
      po[5,i,t+1,7] <- 0
      po[5,i,t+1,8] <- 0
      po[5,i,t+1,9] <- 1-pLF[i,t+1]
      po[5,i,t+1,10] <- 0
      po[5,i,t+1,11] <- 0
 
      po[6,i,t+1,1] <- 0
      po[6,i,t+1,2] <- 0
      po[6,i,t+1,3] <- 0
      po[6,i,t+1,4] <- 0
      po[6,i,t+1,5] <- 0
      po[6,i,t+1,6] <- pLT[t+1]
      po[6,i,t+1,7] <- 0
      po[6,i,t+1,8] <- 0
      po[6,i,t+1,9] <- 1-pLT[t+1]
      po[6,i,t+1,10] <- 0
      po[6,i,t+1,11] <- 0

      po[7,i,t+1,1] <- 0
      po[7,i,t+1,2] <- 0
      po[7,i,t+1,3] <- 0
      po[7,i,t+1,4] <- 0
      po[7,i,t+1,5] <- 0
      po[7,i,t+1,6] <- 0
      po[7,i,t+1,7] <- pLF[i,t+1]
      po[7,i,t+1,8] <- 0
     po[7,i,t+1,9] <- 1-pLF[i,t+1]
      po[7,i,t+1,10] <- 0
      po[7,i,t+1,11] <- 0

      po[8,i,t+1,1] <- 0
      po[8,i,t+1,2] <- 0
      po[8,i,t+1,3] <- 0
      po[8,i,t+1,4] <- 0
      po[8,i,t+1,5] <- 0
      po[8,i,t+1,6] <- 0
      po[8,i,t+1,7] <- 0
      po[8,i,t+1,8] <- pLT[t+1]
      po[8,i,t+1,9] <- 1-pLT[t+1]
      po[8,i,t+1,10] <- 0
      po[8,i,t+1,11] <- 0

      po[9,i,t+1,1] <- 0
      po[9,i,t+1,2] <- 0
      po[9,i,t+1,3] <- 0
      po[9,i,t+1,4] <- 0
      po[9,i,t+1,5] <- 0
      po[9,i,t+1,6] <- 0
      po[9,i,t+1,7] <- 0
      po[9,i,t+1,8] <- 0
      po[9,i,t+1,9] <- 1
      po[9,i,t+1,10] <- 0
      po[9,i,t+1,11] <- 0
      po[10,i,t+1,1] <- 0
      po[10,i,t+1,2] <- 0
      po[10,i,t+1,3] <- 0
      po[10,i,t+1,4] <- 0
      po[10,i,t+1,5] <- 0
      po[10,i,t+1,6] <- 0
      po[10,i,t+1,7] <- 0
      po[10,i,t+1,8] <- 0
      po[10,i,t+1,9] <- 1
      po[10,i,t+1,10] <- 0
      po[10,i,t+1,11] <- 0
      po[11,i,t+1,1] <- 0
      po[11,i,t+1,2] <- 0
      po[11,i,t+1,3] <- 0
      po[11,i,t+1,4] <- 0
      po[11,i,t+1,5] <- 0
      po[11,i,t+1,6] <- 0
      po[11,i,t+1,7] <- 0
      po[11,i,t+1,8] <- 0
      po[11,i,t+1,9] <- 1
      po[11,i,t+1,10] <- 0
      po[11,i,t+1,11] <- 0} #t } #i
 

 
# Likelihood
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- MSmat[i,f[i]]
   for (t in (f[i]+1):nyears){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      MSmat[i,t] ~ dcat(po[z[i,t], i, t,])
      } #t
   } #i
        
# 4.3. Likelihood for fecundity: Poisson regression from the number of surveyed brood
# -------------------------------------------------
for (ttt in 1:nyears){
   NFSE[ttt] ~ dpois(rho[ttt])
   rho[ttt] <- NTSE[ttt]*fecSE[ttt]   } #

} # END of the BUGS model
",fill = TRUE)
sink()
    
# bundle data
dataIPM <- list(MSmat = MSmat, NTSE = NTSE, NFSE = NFSE, nyears = nyears, nind = nind, f = f, x = x, w = w,Placette = Placette,NDVI = NDVI, Time=Time, z = known.state.ms(MSmat,9))

# initial values
initIPM <- function(){list(mean.phiF = runif(4,0.6,0.9), mean.phiT = runif(1,0.8,1), mean.pF = runif(3,0,0.5), mean.pT = runif(1,0.8,1), mean.R = runif(1,0,0.5), mean.L = runif(1,0,1), mean.pLT = runif(1,0,1), mu.fec = runif(1,0,2),mu.im = runif(1,-5,0.5), njSE = c(NA,rep(7,(nyears-1))), n1SE = c(NA,rep(6,(nyears-1))), n2SE = rep(5,(nyears)), n3SE = rep(4,(nyears)), n4SE = rep(3,(nyears)), n5SEa = c(NA,rep(2,(nyears-1))),IMSE = c(NA,rep(1,(nyears-1))),  5SEb = c(NA,rep(2,(nyears-1))), totSE = c(NA,rep(12,(nyears-1))), t3SE = c(NA,rep(2,(nyears-1))), t4SE = c(NA,rep(2,(nyears-1))), t5SEa = c(NA,rep(1,(nyears-1))),t5SEb = c(NA,rep(2,(nyears-1))), alpha = rnorm(1), beta = rnorm(1),eta = rnorm(1),kapa = rnorm(1),z = ms.init.z(MSmat,f))}

# MCMC settings
nc <- 2
nt <- 10
ni <- 100000
nb <- 50000     

# parameters monitored
paraIPM <-c("mean.phiF","mean.phiT","mean.pF","mean.pT","mean.set","mean.R" ,"mean.L","mean.pLT","mean.pLF","mean.im","im","NtotSE","IMSE", "newSE","floSE","lambdaSE","mean.lambdaSE","surF","surT","fecSE","Food", "NDVIF","eta","recF","pTime")

# call winbugs from r (brt 12000 s)
NeoIPMi <- bugs(dataIPM, inits= initIPM, parameters=paraIPM, "NeoIPMi.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
