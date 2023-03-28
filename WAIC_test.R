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

calc.marglike <- function(x){
  vars <- grep("logmarglike", colnames(x$mcmc[[1]])) #find the output that relates to loglike
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


modelstring.birds = "
  model 
{
for (i in 1:n.sites){ 
  log(lambda[i]) <- beta0 + beta1*tree[i]*on[1]
  N[i] ~ dpois(lambda[i])
  log_N[i] <- logdensity.pois(N[i], lambda[i])
  for (t in 1:n.visit){
  logit(p[i,t]) <- alpha0 + alpha1*wind[i,t]*on[2]
  y[i,t] ~ dbin(p[i,t], N[i])
  y.new[i,t] ~ dbin(p[i,t], N[i])
  loglike.waic[i,t] <- logdensity.bin(y[i,t], p[i,t], N[i])
  }
  loglike.new[i] <- sum(loglike.waic[i,])+log_N[i]
  
  
}

beta0 ~ dunif(-3,3)
beta1 ~ dunif(-2,2)
alpha0 ~ dunif(-2,2)
alpha1 ~ dunif(-2,2)
}
"

modelstring.marg = "
model {

for (i in 1:n.sites){ 
  log(lambda[i]) <- beta0 + beta1*tree[i]*on[1]
    for(t in 1:n.visit){
      logit(p[i,t]) <- alpha0 + alpha1*wind[i,t]*on[2]
      } #end t
    for(q in 1:nPossibleN) {
      PrN[i,q] <- dpois(possibleN[q], lambda[i])
        for (t in 1:n.visit){
          logcondlike.y.visit[i,t,q] <- logdensity.bin(y[i,t], p[i,t], possibleN[q])
            } #end t again
      condlike.y[i,q] <- exp(sum(logcondlike.y.visit[i,,q]))
    } #end q
  logmarglike[i] <- log(PrN[i,] %*% condlike.y[i,])
  zeros[i] ~ dpois(-logmarglike[i])
}

beta0 ~ dunif(-3,3)
beta1 ~ dunif(-2,2)
alpha0 ~ dunif(-2,2)
alpha1 ~ dunif(-2,2)
}"


set.seed(20) 

tree <- beta.0 <- beta.1 <- wind <- detect1 <- detect2 <- abund.bird <-  list()
#a produces total prob detection = (.22, .67)
# b produces (.63 to .97)
alpha0_a <- -2
alpha1_a <- -.5
alpha0_b <- -.5
alpha1_b <- -.5
n.sites <- rep(c(15,25,50), each = 100)
WAIC_rank1 <- jointlike_rank1 <- marglike_rank1<- Dsel_rank1 <- WAIC_rank2 <- jointlike_rank2 <- marglike_rank2<- Dsel_rank2<- rep(NA, 300)
WAIC_w1 <- jointlike_w1 <- marg_w1<- Dsel_w1 <- WAIC_w2 <- jointlike_w2 <- marg_w2 <- Dsel_w2 <- rep(NA, 300)
top_waic1 <- top_waic2 <- top_jointlike1 <- top_jointlike2 <- top_marg1 <- top_marg2 <- top_Dsel1 <- top_Dsel2 <- rep(NA, 300)
params <- c("loglike.waic", "loglike.new", "y.new")
tree <- runif(50)
tree <- scale(tree); attr(tree,"scaled:scale")<-NULL;attr(tree,"scaled:center")<-NULL
wind <- runif(50*4, 0,15) # fake wind speeds
wind <- scale(wind)
attr(wind,"scaled:center")<-NULL; attr(wind,"scaled:scale")<-NULL
wind.mat <- matrix(wind, ncol = 4, byrow = 1)

for (i in 275:300){
print(paste("I just started run #", i, sep = ""))
beta.0[[i]] <- runif(1)
beta.1[[i]] <- runif(1, min = -2, max = 2)
abund.bird[[i]] <- rpois(n.sites[i], exp(beta.0[[i]]+beta.1[[i]]*tree[1:n.sites[i]]))
detect1[[i]] <- rbinom(n.sites[i]*4, rep(abund.bird[[i]], each = 4), plogis(alpha0_a+alpha1_a*wind[1:(n.sites[i]*4)]))
detect1[[i]] <- matrix(detect1[[i]], ncol = 4, byrow = T)
detect2[[i]] <- rbinom(n.sites[i]*4, rep(abund.bird[[i]], each = 4), plogis(alpha0_b+alpha1_b*wind[1:(n.sites[i]*4)]))
detect2[[i]] <- matrix(detect2[[i]], ncol = 4, byrow = T)

birds1 <- data.frame(Site = seq(1:n.sites[[i]]), Visit1 = detect1[[i]][,1], Visit2 = detect1[[i]][,2], Visit3 = detect1[[i]][,3], Visit4 = detect1[[i]][,4], Wind1 = wind.mat[1:n.sites[i],1], Wind2 = wind.mat[1:n.sites[i],2], Wind3 = wind.mat[1:n.sites[i],3], Wind4 = wind.mat[1:n.sites[i],4], PercentCover = tree[1:n.sites[i]])

birds2<- data.frame(Site = seq(1:n.sites[[i]]), Visit1 = detect2[[i]][,1], Visit2 = detect2[[i]][,2], Visit3 = detect2[[i]][,3], Visit4 = detect2[[i]][,4], Wind1 = wind.mat[1:n.sites[i],1], Wind2 = wind.mat[1:n.sites[i],2], Wind3 = wind.mat[1:n.sites[i],3], Wind4 = wind.mat[1:n.sites[i],4], PercentCover = tree[1:n.sites[i]])

##### imperfect detection ####
#### Params ####
jd.birds <- list(y = as.matrix(birds1[,c("Visit1", "Visit2", "Visit3", "Visit4")]),
                 wind = as.matrix(birds1[,c("Wind1", "Wind2", "Wind3", "Wind4")]),
                 tree = birds1[,"PercentCover"], 
                 n.sites = n.sites[[i]],
                 n.visit = 4)
jd.marg <-  list(y = as.matrix(birds1[,c("Visit1", "Visit2", "Visit3", "Visit4")]),
                 wind = as.matrix(birds1[,c("Wind1", "Wind2", "Wind3", "Wind4")]),
                 tree = birds1[,"PercentCover"], 
                 n.sites = n.sites[[i]],
                 n.visit = 4,
                 nPossibleN =length(1:(max(abund.bird[[i]]) + 30)), 
                 possibleN= 1:(max(abund.bird[[i]]) + 30), 
                 zeros = rep(0, n.sites[i]))
ji.birbs <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1),
  N = apply(jd.birds$y, 1, max)+2
)}
ji.marg <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1)
)}


### global model/correct model ####
jd.birds$on <- c(1,1)
jags.birds1 <- run.jags(model = modelstring.birds, 
                       inits = ji.birbs,
                       monitor = params, 
                       data = jd.birds, n.chains = 3, burnin =  5000,
                       sample = 2000, method = "parallel")
jd.marg$on <- c(1,1)
jags.marg1 <- run.jags(model = modelstring.marg, 
                        inits = ji.marg,
                        monitor = "logmarglike", 
                        data = jd.marg, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
##### no covariates model ####
jd.birds$on <- c(0,0)
jags.birds2 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg$on <- c(0,0)
jags.marg2 <- run.jags(model = modelstring.marg, 
                       inits = ji.marg,
                       monitor = "logmarglike", 
                       data = jd.marg, n.chains = 3, burnin =  5000,
                       sample = 2000, method = "parallel")

##### tree, no wind ####
jd.birds$on <- c(1,0)
jags.birds3 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg$on <- c(1,0)
jags.marg3 <- run.jags(model = modelstring.marg, 
                       inits = ji.marg,
                       monitor = "logmarglike", 
                       data = jd.marg, n.chains = 3, burnin =  5000,
                       sample = 2000, method = "parallel")

##### wind, no tree ####
jd.birds$on <- c(0,1)
jags.birds4 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs,
                        monitor = params, 
                        data = jd.birds, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg$on <- c(0,1)
jags.marg4 <- run.jags(model = modelstring.marg, 
                       inits = ji.marg,
                       monitor = "logmarglike", 
                       data = jd.marg, n.chains = 3, burnin =  5000,
                       sample = 2000, method = "parallel")

#### WAIC and Joint Likelihood Tables for Imperfect Detection Data ####
mymodels <- mget(ls()[grep("jags.birds", ls())]) #grabs all models in environment
margmodels <- mget(ls()[grep("jags.marg", ls())])
WAIC <- jointlike <- marglike <-  Dsel <- data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC[k,2] <-  calc.waic(mymodels[[k]])
  jointlike[k,2] <- calc.jointlike(mymodels[[k]])
  marglike[k,2] <- calc.marglike(margmodels[[k]])
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

marglike$deltamarglike <- marglike$loglike-min(marglike$loglike)
marglike$rel_like <- exp(-.5*marglike$deltamarglike)
marglike$weight <- marglike$rel_like/sum(marglike$rel_like)
marglike <- marglike[order(-marglike$weight),]  
marglike_rank1[i] <- which(marglike$modname == 1)
marg_w1[i] <- marglike$weight[marglike_rank1[i]]
top_marg1[i] <- marglike[1,1]

##### Perfect detection ####
##### Params2 ####
print(paste("I just started Perfect Detection models for run #", i, sep = ""))

jd.birds2 <- list(y = as.matrix(birds2[,c("Visit1", "Visit2", "Visit3", "Visit4")]),
                 wind = as.matrix(birds2[,c("Wind1", "Wind2", "Wind3", "Wind4")]),
                 tree = birds2[,"PercentCover"], 
                 n.sites = n.sites[i],
                 n.visit = 4)
jd.marg2 <-  list(y = as.matrix(birds2[,c("Visit1", "Visit2", "Visit3", "Visit4")]),
                 wind = as.matrix(birds2[,c("Wind1", "Wind2", "Wind3", "Wind4")]),
                 tree = birds2[,"PercentCover"], 
                 n.sites = n.sites[[i]],
                 n.visit = 4,
                 nPossibleN =length(1:(max(abund.bird[[i]]) + 30)), 
                 possibleN= 1:(max(abund.bird[[i]]) + 30), 
                 zeros = rep(0, n.sites[i]))

ji.birbs2 <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1),
  N = apply(jd.birds2$y, 1, max)+2
)}
ji.marg2 <- function(){list(
  alpha0 = runif(1, 0, 1),
  alpha1= runif(1),
  beta0 = runif(1),
  beta1 = runif(1, -1, 1)
)}

### global model/correct model ####
jd.birds2$on <- c(1,1)
jags.like_1 <-  run.jags(model = modelstring.birds, 
                         inits = ji.birbs2,
                         monitor = params, 
                         data = jd.birds2, n.chains = 3, burnin =  5000,
                         sample = 2000, method = "parallel")
jd.marg2$on <- c(1,1)
j.margs1 <- run.jags(model = modelstring.marg, 
                       inits = ji.marg2,
                       monitor = "logmarglike", 
                       data = jd.marg2, n.chains = 3, burnin =  5000,
                       sample = 2000, method = "parallel")

##### no covariates ####
jd.birds2$on <- c(0,0)
jags.like_2 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs2,
                        monitor = params, 
                        data = jd.birds2, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg2$on <- c(0,0)
j.margs2 <- run.jags(model = modelstring.marg, 
                     inits = ji.marg2,
                     monitor = "logmarglike", 
                     data = jd.marg2, n.chains = 3, burnin =  5000,
                     sample = 2000, method = "parallel")
##### tree, no wind ####
jd.birds2$on <- c(1,0)
jags.like_3 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs2,
                        monitor = params, 
                        data = jd.birds2, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg2$on <- c(1,0)
j.margs3 <- run.jags(model = modelstring.marg, 
                     inits = ji.marg2,
                     monitor = "logmarglike", 
                     data = jd.marg2, n.chains = 3, burnin =  5000,
                     sample = 2000, method = "parallel")

##### no tree, wind ####
jd.birds2$on <- c(0,1)
jags.like_4 <- run.jags(model = modelstring.birds, 
                        inits = ji.birbs2,
                        monitor = params, 
                        data = jd.birds2, n.chains = 3, burnin =  5000,
                        sample = 2000, method = "parallel")
jd.marg2$on <- c(0,1)
j.margs4 <- run.jags(model = modelstring.marg, 
                     inits = ji.marg2,
                     monitor = "logmarglike", 
                     data = jd.marg2, n.chains = 3, burnin =  5000,
                     sample = 2000, method = "parallel")

#### WAIC and Joint Likelihood Tables for Imperfect Detection Data ####
mymodels2 <- mget(ls()[grep("jags.like_", ls())])
margmodels2 <- mget(ls()[grep("j.marg", ls())])
WAIC2 <- jointlike2 <- marglike2 <-  Dsel2 <-data.frame(modname = 1:4, loglike = rep(NA, 4))
for(k in 1:4){
  WAIC2[k,2] <-  calc.waic(mymodels2[[k]])
  jointlike2[k,2] <- calc.jointlike(mymodels2[[k]])
  marglike2[k,2] <- calc.marglike(margmodels2[[k]])
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

Dsel2$deltaDsel <- Dsel2$loglike-min(Dsel2$loglike)
Dsel2$rel_like <- exp(-.5*Dsel2$deltaDsel)
Dsel2$weight <- Dsel2$rel_like/sum(Dsel2$rel_like)
Dsel2 <- Dsel2[order(-Dsel2$weight),]  
Dsel_rank2[i] <- which(Dsel2$modname == 1)
Dsel_w2[i] <- Dsel2$weight[Dsel_rank2[i]]
top_Dsel2[i] <- Dsel2[1,1]

marglike2$deltamarglike <- marglike2$loglike-min(marglike2$loglike)
marglike2$rel_like <- exp(-.5*marglike2$deltamarglike)
marglike2$weight <- marglike2$rel_like/sum(marglike2$rel_like)
marglike2 <- marglike2[order(-marglike2$weight),]  
marglike_rank2[i] <- which(marglike2$modname == 1)
marg_w2[i] <- marglike2$weight[marglike_rank2[i]]
top_marg2[i] <- marglike2[1,1]
}

##### Results ########
dput(abund.bird, "abund.birds.txt")
dput(beta.0, "beta0.txt")
dput(beta.1, "beta1.txt")
dput(detect1, "detect1.txt")
dput(detect2, "detect2.txt")
dput(tree, "tree.txt")
dput(wind, "wind.txt")



results <- data.frame(n.sim = 1:300, n.sites = n.sites, 
                      WAIC_rank1 = WAIC_rank1, WAIC_rank2 = WAIC_rank2, 
                      WAIC_W1 = WAIC_w1, WAIC_W2 = WAIC_w2, 
                      jointlike_rank1 = jointlike_rank1, jointlike_rank2 =  jointlike_rank2,
                      jointlike_w1 = jointlike_w1, jointlike_w2 = jointlike_w2,
                      Dsel_rank1 = Dsel_rank1, Dsel_rank2 =  Dsel_rank2,
                      Dsel_w1 = Dsel_w1, Dsel_w2 = Dsel_w2,
                      marg_rank1 = marglike_rank1, marg_rank2 = marglike_rank2,
                      Marg_w1 = marg_w1, Marg_w2 = marg_w2,
                      topwaic_1 = top_waic1, topwaic_2 = top_waic2, 
                      topjoint_1 = top_jointlike1, topjoint_2 = top_jointlike2,
                      topDsel_1 = top_Dsel1, topDsel_2 = top_Dsel2,
                      topmarg_1 = top_marg1, topmarg_2 = top_marg2 )

#dput(results2, "results.txt")

###### Results #######
#setwd("~/Desktop/U_Georgia/Chandler_Meetings/WAIC_Musing/WAIC")
results <- dget("results.txt")
library(dplyr)
results %>% 
  summarize(WAIC_rank1= sum(WAIC_rank1 == 1),
            WAICj_rank1= sum(jointlike_rank1 == 1),
            Marg_rank1 = sum(marg_rank1 == 1),
            Dsel_rank1 = sum(Dsel_rank1 == 1),
            WAIC_rank2= sum(WAIC_rank2 == 1),
            WAICj_rank2= sum(jointlike_rank2 == 1),
            Marg_rank2 = sum(marg_rank2 == 1),
            Dsel_rank2 = sum(Dsel_rank2 == 1)) 

results %>% 
  summarize(WAIC_rank1= sum(WAIC_rank1 == 1) + sum(WAIC_rank2 == 1),
      WAICj_rank1= sum(jointlike_rank1 == 1)+sum(jointlike_rank2 == 1),
      Marg_rank1 = sum(marg_rank1 == 1) + sum(marg_rank2 == 1),
      Dsel_rank1 = sum(Dsel_rank1 == 1) + sum(Dsel_rank2 == 1))


ranks <-   results %>% 
  group_by(n.sites) %>%
  summarize(WAIC_rank1= sum(WAIC_rank1 == 1),
            WAICj_rank1= sum(jointlike_rank1 == 1),
            Marg_rank1 = sum(marg_rank1 == 1),
            Dsel_rank1 = sum(Dsel_rank1 == 1),
            WAIC_rank2= sum(WAIC_rank2 == 1),
            WAICj_rank2= sum(jointlike_rank2 == 1),
            Marg_rank2 = sum(marg_rank2 == 1),
            Dsel_rank2 = sum(Dsel_rank2 == 1)) 
ranks

weights <-  results %>% 
  group_by(n.sites) %>%
  summarize(WAIC_w1= mean(WAIC_W1),
            WAIC_sd1 = sd(WAIC_W1),
            WAICj_w1= mean(jointlike_w1),
            WAICj_sd1 = sd(jointlike_w1),
            Marg_weight1 = mean(Marg_w1),
            Marg_std1 = sd(Marg_w1),
            Dsel_w1= mean(Dsel_w1),
            Dsel_sd1 = sd(Dsel_w1),
            WAIC_w2= mean(WAIC_W2),
            WAIC_sd2 = sd(WAIC_W2),
            WAICj_w2= mean(jointlike_w2),
            WAICj_sd2 = sd(jointlike_w2),
            Marg_weight2 = mean(Marg_w2),
            Marg_std2 = sd(Marg_w2),
            Dsel_w2= mean(Dsel_w2),
            Dsel_sd2 = sd(Dsel_w2)) 

weights

gg.results <- data.frame(sites = paste(rep(results$n.sites,8), " sites", sep = ""), rank = c(results$topwaic_1, results$topwaic_2, results$topjoint_1, results$topjoint_2, results$topmarg_1, results$topmarg_2, results$topDsel_1, results$topDsel_2), Method = rep(c("WAIC", "WAICj", "Marg", "Dsel"), each = nrow(results)*2), detection = rep(c(rep("Low Detection", nrow(results)), rep("High Detection", nrow(results))), 4)) 
head(gg.results)

library(ggplot2)
gg.results$rank <- as.factor(gg.results$rank)
test <- gg.results %>%
  group_by(sites, Method, detection, rank, .drop = F) %>%
  count(rank, .drop = F)
test$n <- (test$n)/100
test <- test[!is.na(test$rank),]
(pp <- ggplot(test, aes(x = rank, group = Method, fill = Method))+
  #geom_histogram(aes(y = ..density..), binwidth = .8, bins = 4, alpha = .25, position = "dodge")+
  geom_col(aes(y = n), alpha = .85, position = "dodge")+
  facet_wrap(~detection+sites, nrow = 2)+
  theme_classic()+
  xlab("Top Ranked Model")+
  ylab("Frequency")+
  ylim(0,1)
  #+theme(strip.text.x = element_blank())
  )


library(egg)
pp<- tag_facet(pp, open = "", close = "", 
          tag_pool = c("A15", "A25", "A50", "B15", "B25", "B50", hjust = -.45, fontface = 2, size = 5))
pp
#ggsave("Ranking_WAIC.jpeg", pp, dpi = 400)

par(mfrow = c(2,4))
weights <- data.frame(n.sites = results$n.sites, WAICj = results$jointlike_w1, WAIC = results$WAIC_W1, marg = results$Marg_w1, Dsel = results$Dsel_w1)
boxplot(WAIC ~ n.sites, weights, ylim = c(0,1), main = "WAIC  \nLow Det")
boxplot(WAICj ~ n.sites, weights,ylim = c(0,1), main = "WAICj  \nLow Det")
boxplot(marg ~ n.sites, weights,ylim = c(0,1), main = "marg  \nLow Det")
boxplot(Dsel ~ n.sites, weights,ylim = c(0,1), main = "Dsel  \nLow Det")
weights2 <- data.frame(n.sites = results$n.sites, WAICj = results$jointlike_w2, WAIC = results$WAIC_W2, marg = results$Marg_w2, Dsel = results$Dsel_w2)
boxplot(WAIC ~ n.sites, weights2, ylim = c(0,1), main = "WAIC \nHigh Det")
boxplot(WAICj ~ n.sites, weights2,ylim = c(0,1), main = "WAICj \nHigh Det")
boxplot(marg ~ n.sites, weights2,ylim = c(0,1), main = "marg \nHigh Det")
boxplot(Dsel ~ n.sites, weights2,ylim = c(0,1), main = "Dsel \nHigh Det")



