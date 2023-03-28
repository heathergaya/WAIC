library(runjags)

calc.waic <- function(x){
  vars <- grep("loglike.waic", colnames(x$mcmc[[1]])) #find the output that relates to loglike
  like <- as.matrix(x$mcmc[,vars,]) 
  fbar <- colMeans(exp(like)) #mean log-likelihood 
  Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
  WAIC<- -2*sum(log(fbar))+2*Pw
  return(WAIC)
}
calc.jointlike <- function(x){
  vars <- grep("loglike.new", colnames(x$mcmc[[1]])) #find the output that relates to loglike
  like <- as.matrix(x$mcmc[,vars,]) 
  fbar <- colMeans(exp(like)) #mean likelihood 
  Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
  WAIC_ish<- -2*sum(log(fbar))+2*Pw
  return(WAIC_ish)
}

calc.jointlike2 <- function(x){
  vars <- grep("loglike.new2", colnames(x$mcmc[[1]])) #find the output that relates to loglike
  like <- as.matrix(x$mcmc[,vars,]) 
  fbar <- colMeans(exp(like)) #mean likelihood 
  Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
  WAIC_ish<- -2*sum(log(fbar))+2*Pw
  return(WAIC_ish)
}

calc.Dsel <- function(x){
  vars <- grep("y.new", colnames(x$mcmc[[1]]))
  y.save <- as.matrix(x$mcmc[,vars,]) 
  y.mean <- apply(y.save,2,mean)
  y.real <- c(jd.birds$y)
  sum.1 <- sum((y.real-y.mean)^2) #want to minimize this
  sum.2 <- sum(apply(y.save,2,var))
  D_sel = sum.1+sum.2
  return(D_sel)
}


modelstring.birds_ac = "
  model 
{
for (i in 1:n.sites){ 
  #first year
  log(lambda[i,1]) <- beta0 + beta1*tree[i,1]*on[1]
  N[i,1] ~ dpois(lambda[i,1])
  log_N[i,1] <- logdensity.pois(N[i,1], lambda[i,1])
  psi[i,1] <- 1 #make JAGS happy
   
  for (t in 2:n.years){
  psi[i,t] <- exp(b0 + b1*t + b2*tree[i,t]*on[1]) #simple time trend; same for all sites
  lambda[i,t] <- lambda[i,t-1]*psi[i,t] 
  N[i,t] ~ dpois(lambda[i,t]) 
  log_N[i,t] <- logdensity.pois(N[i,t], lambda[i,t])

  } #end t

  for (t in 1:n.years){
    for (j in 1:n.visit){
      logit(p[i,j,t]) <- alpha0 + alpha1*wind[i,j,t]*on[2]
      y[i,j,t] ~ dbin(p[i,j,t], N[i,t])
      y.new[i,j,t] ~ dbin(p[i,j,t], N[i,t])
      loglike.waic[i,j,t] <- logdensity.bin(y[i,j,t], p[i,j,t], N[i,t])
    } #end j
  loglike.new1[i,t] <- sum(loglike.waic[i,,t]) +log_N[i,t] #sum all, despite autocorrelation
  loglike.new2[i,t] <- mean(loglike.waic[i,,t]) +log_N[i,t] #only take the average detection for each site
  } #end t 
  
} #end i

b0 ~ dnorm(0,1) 
b1 ~ dnorm(0,1) #time trend
b2 ~ dnorm(0,1) #tree 
beta0 ~ dunif(-3,3)
beta1 ~ dunif(-2,2)
alpha0 ~ dunif(-2,2)
alpha1 ~ dunif(-2,2)
}
"

set.seed(20) 

tree <- beta.0 <- beta.1 <- wind <- detect1 <- detect2 <- detect1b <- detect2b <- abund.bird <- abund.bird2 <-  list()
#a produces total prob detection = (.22, .67)
# b produces (.63 to .97)
alpha0_a <- -2
alpha1_a <- -.5
alpha0_b <- -.5
alpha1_b <- -.5
n.sites <- 50
n.years <- 5
b0 <- -1
b1 <- .3 
b1_b <- .5 #not realistic but good for testing
b2 <- .1
params <- c("loglike.waic", "loglike.new1", "loglike.new2", "y.new")
tree <- runif(50*n.years)
tree <- scale(tree); attr(tree,"scaled:scale")<-NULL;attr(tree,"scaled:center")<-NULL
tree.mat <- array(tree, dim = c(n.sites, n.years)) #sites, years 
wind <- runif(50*4*n.years, 0,15) # fake wind speeds
wind <- scale(wind)
attr(wind,"scaled:center")<-NULL; attr(wind,"scaled:scale")<-NULL
wind.mat <- array(wind, dim = c(n.sites, 4, n.years)) #sites, visits, years 

for(i in 1:100){
print(paste("I just started run #", i, sep = ""))
beta.0[[i]] <- runif(1)
beta.1[[i]] <- runif(1, min = -2, max = 2)

abund.bird[[i]] <- abund.bird2[[i]] <- lambda <- lambda2 <- psi <- psi2 <- N <- N2 <- matrix(NA, nrow = n.sites, ncol = n.years)
abund.bird[[i]][,1] <- N[,1] <- abund.bird2[[i]][,1] <- N2[,1] <-  rpois(n.sites, exp(beta.0[[i]]+beta.1[[i]]*tree.mat[,1]))
lambda[,1] <- lambda2[,1] <- exp(beta.0[[i]]+beta.1[[i]]*tree.mat[,1])

for(t in 2:n.years){
  #small trend
  psi[,t] <- exp(b0+b1*t+b2*tree.mat[,t])
  lambda[,t] <- lambda[,t-1]*psi[,t]
  N[,t] <- rpois(n.sites, lambda[,t])
  abund.bird[[i]][,t] <- N[,t]
  
  #large time trend
  psi2[,t] <- exp(b0+b1_b*t+b2*tree.mat[,t])
  lambda2[,t] <- lambda2[,t-1]*psi2[,t]
  N2[,t] <- rpois(n.sites, lambda2[,t])
  abund.bird2[[i]][,t] <- N2[,t]
}

det <- det2 <- detb <- det2b <- array(NA, c(n.sites, 4, n.years))
for(k in 1:nrow(det)){
  for(t in 1:n.years){
    det[k,,t] <- rbinom(4, abund.bird[[i]][k,t], plogis(alpha0_a+alpha1_a*wind.mat[k,,t]))
    det2[k,,t] <- rbinom(4, abund.bird[[i]][k,t], plogis(alpha0_b+alpha1_b*wind.mat[k,,t]))
    detb[k,,t] <- rbinom(4, abund.bird2[[i]][k,t], plogis(alpha0_a+alpha1_a*wind.mat[k,,t]))
    det2b[k,,t] <- rbinom(4, abund.bird2[[i]][k,t], plogis(alpha0_b+alpha1_b*wind.mat[k,,t]))
  }
}
detect1[[i]] <- det
detect2[[i]] <- det2
detect1b[[i]] <- detb
detect2b[[i]] <- det2b

WAIC_rank1 <- jointlike_rank1 <- jointlike2_rank1<- Dsel_rank1 <- WAIC_rank2 <- jointlike_rank2 <- jointlike2_rank2 <- Dsel_rank2<- rep(NA, 300)
WAIC_w1 <- jointlike_w1 <- jointlike2_w1<- Dsel_w1 <- WAIC_w2 <- jointlike_w2 <- jointlike2_w2 <- Dsel_w2 <- rep(NA, 300)
top_waic1 <- top_waic2 <- top_jointlike1 <- top_jointlike2 <- top_jointlike21 <- top_jointlike22 <- top_Dsel1 <- top_Dsel2 <- rep(NA, 300)

#### Small Time Trend ########
##### imperfect detection ####
#### Params ####
jd.birds <- list(y = det,
                 wind = wind.mat,
                 tree = tree.mat, 
                 n.sites = n.sites,
                 n.visit = 4, n.years = n.years)
N.init <- N
ji.birbs <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1),
  N = N.init,
  b0 = b0
)}

### global model/correct model ####
jd.birds$on <- c(1,1)
jags.birds1 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### no covariates model ####
jd.birds$on <- c(0,0)
jags.birds2 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### tree, no wind ####
jd.birds$on <- c(1,0)
jags.birds3 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### wind, no tree ####
jd.birds$on <- c(0,1)
jags.birds4 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")


mymodels <- mget(ls()[grep("jags.birds", ls())]) #grabs all models in environment
WAIC <- jointlike <- jointlike2 <-  Dsel <- data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC[k,2] <-  calc.waic(mymodels[[k]])
  jointlike[k,2] <- calc.jointlike(mymodels[[k]])
  jointlike2[k,2] <- calc.jointlike2(mymodels[[k]])
  Dsel[k,2] <- calc.Dsel(mymodels[[k]])
}

WAIC$deltaWAIC <- WAIC$loglike-min(WAIC$loglike)
WAIC$rel_like <- exp(-.5*WAIC$deltaWAIC)
WAIC$weight <- WAIC$rel_like/sum(WAIC$rel_like)
WAIC <- WAIC[order(-WAIC$weight),]  
WAIC_rank1[i] <- which(WAIC$modname == 1)
WAIC_w1[i] <- WAIC$weight[WAIC_rank1[i]]
top_waic1[i] <- WAIC[1,1]
jointlike$deltajointlike <- jointlike$loglike-min(jointlike$loglike)
jointlike$rel_like <- exp(-.5*jointlike$deltajointlike)
jointlike$weight <- jointlike$rel_like/sum(jointlike$rel_like)
jointlike <- jointlike[order(-jointlike$weight),]  
jointlike_rank1[i] <- which(jointlike$modname == 1)
jointlike_w1[i] <- jointlike$weight[jointlike_rank1[i]]
top_jointlike1[i] <- jointlike[1,1]

Dsel$deltaDsel <- Dsel$loglike-min(Dsel$loglike)
Dsel$rel_like <- exp(-.5*Dsel$deltaDsel)
Dsel$weight <- Dsel$rel_like/sum(Dsel$rel_like)
Dsel <- Dsel[order(-Dsel$weight),]  
Dsel_rank1[i] <- which(Dsel$modname == 1)
Dsel_w1[i] <- Dsel$weight[Dsel_rank1[i]]
top_Dsel1[i] <- Dsel[1,1]

jointlike2$deltajointlike2 <- jointlike2$loglike-min(jointlike2$loglike)
jointlike2$rel_like <- exp(-.5*jointlike2$deltajointlike2)
jointlike2$weight <- jointlike2$rel_like/sum(jointlike2$rel_like)
jointlike2 <- jointlike2[order(-jointlike2$weight),]  
jointlike2_rank1[i] <- which(jointlike2$modname == 1)
jointlike2_w1[i] <- jointlike2$weight[jointlike2_rank1[i]]
top_jointlike21[i] <- jointlike2[1,1]

##### High detection ####
##### Params2 ####
print(paste("I just started High Detection models for run #", i, sep = ""))

jd.birds <- list(y = det2,
                 wind = wind.mat,
                 tree = tree.mat, 
                 n.sites = n.sites,
                 n.visit = 4, n.years = n.years)

### global model/correct model ####
jd.birds$on <- c(1,1)
jags.birds1 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### no covariates model ####
jd.birds$on <- c(0,0)
jags.birds2 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### tree, no wind ####
jd.birds$on <- c(1,0)
jags.birds3 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### wind, no tree ####
jd.birds$on <- c(0,1)
jags.birds4 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")


#### WAIC and Joint Likelihood Tables for "perfect" Detection Data ####
mymodels2 <- mget(ls()[grep("jags.birds", ls())])
WAIC2 <- jointlike2 <- jointlike22 <-  Dsel2 <-data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC2[k,2] <-  calc.waic(mymodels2[[k]])
  jointlike2[k,2] <- calc.jointlike(mymodels2[[k]])
  jointlike22[k,2] <- calc.jointlike2(mymodels2[[k]])
  Dsel2[k,2] <- calc.Dsel(mymodels2[[k]])
}

WAIC2$deltaWAIC <- WAIC2$loglike-min(WAIC2$loglike)
WAIC2$rel_like <- exp(-.5*WAIC2$deltaWAIC)
WAIC2$weight <- WAIC2$rel_like/sum(WAIC2$rel_like)
WAIC2 <- WAIC2[order(-WAIC2$weight),]  
WAIC_rank2[i] <- which(WAIC2$modname == 1)
WAIC_w2[i] <- WAIC2$weight[WAIC_rank2[i]]
top_waic2[i] <- WAIC2[1,1]

jointlike2$deltajointlike <- jointlike2$loglike-min(jointlike2$loglike)
jointlike2$rel_like <- exp(-.5*jointlike2$deltajointlike)
jointlike2$weight <- jointlike2$rel_like/sum(jointlike2$rel_like)
jointlike2 <- jointlike2[order(-jointlike2$weight),]  
jointlike_rank2[i] <- which(jointlike2$modname == 1)
jointlike_w2[i] <- jointlike2$weight[jointlike_rank2[i]]
top_jointlike2[i] <- jointlike2[1,1]

jointlike22$deltajointlike2 <- jointlike22$loglike-min(jointlike22$loglike)
jointlike22$rel_like <- exp(-.5*jointlike22$deltajointlike2)
jointlike22$weight <- jointlike22$rel_like/sum(jointlike22$rel_like)
jointlike22 <- jointlike22[order(-jointlike22$weight),]  
jointlike2_rank2[i] <- which(jointlike22$modname == 1)
jointlike2_w2[i] <- jointlike22$weight[jointlike2_rank2[i]]
top_jointlike22[i] <- jointlike22[1,1]

Dsel2$deltaDsel <- Dsel2$loglike-min(Dsel2$loglike)
Dsel2$rel_like <- exp(-.5*Dsel2$deltaDsel)
Dsel2$weight <- Dsel2$rel_like/sum(Dsel2$rel_like)
Dsel2 <- Dsel2[order(-Dsel2$weight),]  
Dsel_rank2[i] <- which(Dsel2$modname == 1)
Dsel_w2[i] <- Dsel2$weight[Dsel_rank2[i]]
top_Dsel2[i] <- Dsel2[1,1]

#### End of small time trend models. 
##### Large Time Trend #######

WAIC_rank1b <- jointlike_rank1b <- jointlike2_rank1b <- Dsel_rank1b <- WAIC_rank2b <- jointlike_rank2b <- jointlike2_rank2b <- Dsel_rank2b <- rep(NA, 300)
WAIC_w1b <- jointlike_w1b <- jointlike2_w1b <- Dsel_w1b <- WAIC_w2b <- jointlike_w2b <- jointlike2_w2b <- Dsel_w2b <- rep(NA, 300)
top_waic1b <- top_waic2b <- top_jointlike1b <- top_jointlike2b <- top_jointlike21b <- top_jointlike22b <- top_Dsel1b <- top_Dsel2b <- rep(NA, 300)

##### imperfect detection ####
#### Params ####
jd.birds <- list(y = detb,
                 wind = wind.mat,
                 tree = tree.mat, 
                 n.sites = n.sites,
                 n.visit = 4, n.years = n.years)
N.init <- N2
ji.birbs <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1),
  N = N.init,
  b0 = b0
)}

### global model/correct model ####
jd.birds$on <- c(1,1)
jags.birds1 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### no covariates model ####
jd.birds$on <- c(0,0)
jags.birds2 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### tree, no wind ####
jd.birds$on <- c(1,0)
jags.birds3 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### wind, no tree ####
jd.birds$on <- c(0,1)
jags.birds4 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")


mymodels <- mget(ls()[grep("jags.birds", ls())]) #grabs all models in environment
WAIC <- jointlike <- jointlike2 <-  Dsel <- data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC[k,2] <-  calc.waic(mymodels[[k]])
  jointlike[k,2] <- calc.jointlike(mymodels[[k]])
  jointlike2[k,2] <- calc.jointlike2(mymodels[[k]])
  Dsel[k,2] <- calc.Dsel(mymodels[[k]])
}

WAIC$deltaWAIC <- WAIC$loglike-min(WAIC$loglike)
WAIC$rel_like <- exp(-.5*WAIC$deltaWAIC)
WAIC$weight <- WAIC$rel_like/sum(WAIC$rel_like)
WAIC <- WAIC[order(-WAIC$weight),]  
WAIC_rank1b[i] <- which(WAIC$modname == 1)
WAIC_w1b[i] <- WAIC$weight[WAIC_rank1b[i]]
top_waic1b[i] <- WAIC[1,1]
jointlike$deltajointlike <- jointlike$loglike-min(jointlike$loglike)
jointlike$rel_like <- exp(-.5*jointlike$deltajointlike)
jointlike$weight <- jointlike$rel_like/sum(jointlike$rel_like)
jointlike <- jointlike[order(-jointlike$weight),]  
jointlike_rank1b[i] <- which(jointlike$modname == 1)
jointlike_w1b[i] <- jointlike$weight[jointlike_rank1b[i]]
top_jointlike1b[i] <- jointlike[1,1]

Dsel$deltaDsel <- Dsel$loglike-min(Dsel$loglike)
Dsel$rel_like <- exp(-.5*Dsel$deltaDsel)
Dsel$weight <- Dsel$rel_like/sum(Dsel$rel_like)
Dsel <- Dsel[order(-Dsel$weight),]  
Dsel_rank1b[i] <- which(Dsel$modname == 1)
Dsel_w1b[i] <- Dsel$weight[Dsel_rank1b[i]]
top_Dsel1b[i] <- Dsel[1,1]

jointlike2$deltajointlike2 <- jointlike2$loglike-min(jointlike2$loglike)
jointlike2$rel_like <- exp(-.5*jointlike2$deltajointlike2)
jointlike2$weight <- jointlike2$rel_like/sum(jointlike2$rel_like)
jointlike2 <- jointlike2[order(-jointlike2$weight),]  
jointlike2_rank1b[i] <- which(jointlike2$modname == 1)
jointlike2_w1b[i] <- jointlike2$weight[jointlike2_rank1b[i]]
top_jointlike21b[i] <- jointlike2[1,1]

##### High detection ####
##### Params2 ####
print(paste("I just started High Detection models for run #", i, sep = ""))

jd.birds <- list(y = det2b,
                 wind = wind.mat,
                 tree = tree.mat, 
                 n.sites = n.sites,
                 n.visit = 4, n.years = n.years)

### global model/correct model ####
jd.birds$on <- c(1,1)
jags.birds1 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### no covariates model ####
jd.birds$on <- c(0,0)
jags.birds2 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### tree, no wind ####
jd.birds$on <- c(1,0)
jags.birds3 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")

##### wind, no tree ####
jd.birds$on <- c(0,1)
jags.birds4 <- run.jags(model = modelstring.birds_ac, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")


#### WAIC and Joint Likelihood Tables for "perfect" Detection Data ####
mymodels2 <- mget(ls()[grep("jags.birds", ls())])
WAIC2 <- jointlike2 <- jointlike22 <-  Dsel2 <-data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC2[k,2] <-  calc.waic(mymodels2[[k]])
  jointlike2[k,2] <- calc.jointlike(mymodels2[[k]])
  jointlike22[k,2] <- calc.jointlike2(mymodels2[[k]])
  Dsel2[k,2] <- calc.Dsel(mymodels2[[k]])
}

WAIC2$deltaWAIC <- WAIC2$loglike-min(WAIC2$loglike)
WAIC2$rel_like <- exp(-.5*WAIC2$deltaWAIC)
WAIC2$weight <- WAIC2$rel_like/sum(WAIC2$rel_like)
WAIC2 <- WAIC2[order(-WAIC2$weight),]  
WAIC_rank2b[i] <- which(WAIC2$modname == 1)
WAIC_w2b[i] <- WAIC2$weight[WAIC_rank2b[i]]
top_waic2b[i] <- WAIC2[1,1]

jointlike2$deltajointlike <- jointlike2$loglike-min(jointlike2$loglike)
jointlike2$rel_like <- exp(-.5*jointlike2$deltajointlike)
jointlike2$weight <- jointlike2$rel_like/sum(jointlike2$rel_like)
jointlike2 <- jointlike2[order(-jointlike2$weight),]  
jointlike_rank2b[i] <- which(jointlike2$modname == 1)
jointlike_w2b[i] <- jointlike2$weight[jointlike_rank2b[i]]
top_jointlike2b[i] <- jointlike2[1,1]

jointlike22$deltajointlike2 <- jointlike22$loglike-min(jointlike22$loglike)
jointlike22$rel_like <- exp(-.5*jointlike22$deltajointlike2)
jointlike22$weight <- jointlike22$rel_like/sum(jointlike22$rel_like)
jointlike22 <- jointlike22[order(-jointlike22$weight),]  
jointlike2_rank2b[i] <- which(jointlike22$modname == 1)
jointlike2_w2b[i] <- jointlike22$weight[jointlike2_rank2b[i]]
top_jointlike22b[i] <- jointlike22[1,1]

Dsel2$deltaDsel <- Dsel2$loglike-min(Dsel2$loglike)
Dsel2$rel_like <- exp(-.5*Dsel2$deltaDsel)
Dsel2$weight <- Dsel2$rel_like/sum(Dsel2$rel_like)
Dsel2 <- Dsel2[order(-Dsel2$weight),]  
Dsel_rank2b[i] <- which(Dsel2$modname == 1)
Dsel_w2b[i] <- Dsel2$weight[Dsel_rank2b[i]]
top_Dsel2b[i] <- Dsel2[1,1]

} #end test

##### Results ########
dput(abund.bird, "abund.birds_weaktrend.txt")
dput(abund.bird2, "abund.birds_strongtrend.txt")
dput(beta.0, "beta0_markov.txt")
dput(beta.1, "beta1_markov.txt")
dput(detect1, "detect1_weaktrend.txt")
dput(detect2, "detect2_weaktrend.txt")
dput(detect1b, "detect1_strongtrend.txt")
dput(detect2b, "detect2_strongtrend.txt")
dput(tree.mat, "tree_markov.txt")
dput(wind.mat, "wind_markov.txt")

results_weaktrend <- data.frame(n.sim = 1:100, n.sites = n.sites, 
                      WAIC_rank1 = WAIC_rank1, WAIC_rank2 = WAIC_rank2, 
                      WAIC_W1 = WAIC_w1, WAIC_W2 = WAIC_w2, 
                      jointlike_rank1 = jointlike_rank1, jointlike_rank2 =  jointlike_rank2,
                      jointlike_w1 = jointlike_w1, jointlike_w2 = jointlike_w2,
                      Dsel_rank1 = Dsel_rank1, Dsel_rank2 =  Dsel_rank2,
                      Dsel_w1 = Dsel_w1, Dsel_w2 = Dsel_w2,
                      jointlike2rank1 = jointlike2_rank1, jointlike2rank2 = jointlike2_rank2,
                      jointlike2w1 = jointlike2w1, jointlike2w2 = jointlike2w2,
                      topwaic_1 = top_waic1, topwaic_2 = top_waic2, 
                      topjoint_1 = top_jointlike1, topjoint_2 = top_jointlike2,
                      topDsel_1 = top_Dsel1, topDsel_2 = top_Dsel2,
                      topjointlike2_1 = top_jointlike21, topjointlike2_2 = top_jointlike22 )

dput(results_weaktrend, "results_weaktrend.txt")

results_strongtrend <- data.frame(n.sim = 1:100, n.sites = n.sites, 
                                  WAIC_rank1b = WAIC_rank1b, WAIC_rank2b = WAIC_rank2b, 
                                  WAIC_W1b = WAIC_w1b, WAIC_W2b = WAIC_w2b, 
                                  jointlike_rank1b = jointlike_rank1b, jointlike_rank2b =  jointlike_rank2b,
                                  jointlike_w1b = jointlike_w1b, jointlike_w2b = jointlike_w2b,
                                  Dsel_rank1b = Dsel_rank1b, Dsel_rank2b =  Dsel_rank2b,
                                  Dsel_w1b = Dsel_w1b, Dsel_w2b = Dsel_w2b,
                                  jointlike2rank1b = jointlike2_rank1b, jointlike2rank2b = jointlike2_rank2b,
                                  jointlike2w1b = jointlike2w1b, jointlike2w2b = jointlike2w2b,
                                  topwaic_1b = top_waic1b, topwaic_2b = top_waic2b, 
                                  topjoint_1b = top_jointlike1b, topjoint_2b = top_jointlike2b,
                                  topDsel_1b = top_Dsel1b, topDsel_2b = top_Dsel2b,
                                  topjointlike2_1b = top_jointlike21b, topjointlike2_2b = top_jointlike22b )

dput(results_strongtrend, "results_strongtrend.txt")




######## unused  ####
# modelstring.birds_dm = "
#   model 
# {
# for (i in 1:n.sites){ 
# #first year
#   log(lambda[i,1]) <- beta0 + beta1*tree[i]*on[1]
#   N[i,1] ~ dpois(lambda[i,1])
#   log_N[i,1] <- logdensity.pois(N[i,1], lambda[i,1])
#   for (t in 2:n.years){
#     lambda[i,t] <- lambda[i,t-1]*phi + lambda[i,t-1]*gamma
#     
#     S[i,t] ~ dbin(phi, N[i,t-1]) #Survival
#     log_S[i,t] <- logdensity.bin(S[i,t], phi, N[i,t-1])
#     
#     G[i,t] ~ dpois(lambda[i,t-1]*gamma) #recruitment
#     log_G[i,t] <- logdensity.pois(G[i,t], lambda[i,t-1]*gamma)
#     
#     N[i,t] <- S[i,t] + G[i,t]
#     log_N[i,t] <- log_S[i,t]+log_G[i,t]
# 
#   } #end t
# 
#   for (t in 1:n.years){
#     for (j in 1:n.visit){
#       logit(p[i,j,t]) <- alpha0 + alpha1*wind[i,j,t]*on[2]
#       y[i,j,t] ~ dbin(p[i,j,t], N[i,t])
#       y.new[i,j,t] ~ dbin(p[i,j,t], N[i,t])
#       loglike.waic[i,j,t] <- logdensity.bin(y[i,j,t], p[i,j,t], N[i,t])
#     } #end j
#   loglike.new1[i,t] <- sum(loglike.waic[i,,t]) +log_N[i,t] #sum all, despite autocorrelation
#   loglike.new2[i,t] <- mean(loglike.waic[i,,t]) +log_N[i,t] #only take the average detection for each site
#   } #end t 
#   
# }
# 
# phi ~ dunif(0,1) # survival
# gamma.0 ~ dnorm(0,1) #recruitment
# log(gamma) <- gamma.0
# beta0 ~ dunif(-3,3)
# beta1 ~ dunif(-2,2)
# alpha0 ~ dunif(-2,2)
# alpha1 ~ dunif(-2,2)
# }
# "
# 
# set.seed(20) 
# 
# tree <- beta.0 <- beta.1 <- wind <- detect1 <- detect2 <- abund.bird <-  list()
# #a produces total prob detection = (.22, .67)
# # b produces (.63 to .97)
# alpha0_a <- -2
# alpha1_a <- -.5
# alpha0_b <- -.5
# alpha1_b <- -.5
# n.sites <- 50
# n.years <- 5
# 
# params <- c("loglike.waic", "loglike.new1", "loglike.new2", "y.new")
# tree <- runif(50)
# tree <- scale(tree); attr(tree,"scaled:scale")<-NULL;attr(tree,"scaled:center")<-NULL
# wind <- runif(50*4*5, 0,15) # fake wind speeds
# wind <- scale(wind)
# attr(wind,"scaled:center")<-NULL; attr(wind,"scaled:scale")<-NULL
# wind.mat <- array(wind, dim = c(n.sites, 4, n.years)) #sites, visits, years 
# 
# i <- 1
# print(paste("I just started run #", i, sep = ""))
# beta.0[[i]] <- runif(1)
# beta.1[[i]] <- runif(1, min = -2, max = 2)
# phi <- .8
# gamma0 <- -.2
# gamma <- exp(gamma0) #.8
# abund.bird[[i]] <- lambda <- S <- G <- N<- matrix(NA, nrow = n.sites, ncol = n.years)
# abund.bird[[i]][,1] <- N[,1] <-  rpois(n.sites, exp(beta.0[[i]]+beta.1[[i]]*tree))
# lambda[,1] <- exp(beta.0[[i]]+beta.1[[i]]*tree)
# for(t in 2:n.years){
#   lambda[,t] <- lambda[,t-1]*phi + lambda[,t-1]*gamma
#   S[,t] <- rbinom(n.sites, N[,t-1], phi)
#   G[,t] <- rpois(n.sites, lambda[i,t-1]*gamma)
#   
#   N[,t] <- S[,t] + G[,t]
#   abund.bird[[i]][,t] <- N[,t]
# }
# 
# det <- det2 <- array(NA, c(n.sites, 4, n.years))
# for(k in 1:nrow(det)){
#   for(t in 1:n.years){
#     det[k,,t] <- rbinom(4, abund.bird[[i]][k,t], plogis(alpha0_a+alpha1_a*wind.mat[k,,t]))
#     det2[k,,t] <- rbinom(4, abund.bird[[i]][k,t], plogis(alpha0_b+alpha1_b*wind.mat[k,,t]))
#   }
# }
# detect1[[i]] <- det
# detect2[[i]] <- det2
# 
# ##### imperfect detection ####
# #### Params ####
# jd.birds <- list(y = det,
#                  wind = wind.mat,
#                  tree = c(tree), 
#                  n.sites = n.sites,
#                  n.visit = 4, n.years = n.years)
# N.init <- array(NA, c(n.sites, n.years))
# N.init[,1] <- N[,1]
# G.init <- G
# S.init <- S
# ji.birbs <- function(){list(
#   alpha0 = runif(1, 0, 1),
#   alpha1= runif(1),
#   beta0 = runif(1),
#   beta1 = runif(1, -1, 1),
#   N = N.init,
#   gamma.0 = gamma0,
#   phi = phi,
#   G = G,
#   S = S
# )}