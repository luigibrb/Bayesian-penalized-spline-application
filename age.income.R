library(R2WinBUGS)
library(R2jags)
set.seed(1974)
# set.seed(3241)


data(age.income, package = "SemiPar")
n<-dim(age.income)[1]
covariate<-age.income$age
response<-age.income$log.income
X<-cbind(rep(1,n),covariate)


#### Inizializzazione ####

# plot(covariate, response)

num.knots<-20
knots<-quantile(unique(covariate),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])

inits.b<-rep(0,num.knots)
inits<-function(){
  list("beta"=c(0,0),"b"=inits.b,"taub"=0.01,"taueps"=0.01)}

parameters<-list("lambda","sigmab","sigmaeps","beta","b","ystar")

Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

data<-list("response"=response,"X"=X,"Z"=Z,"n"=n,"num.knots"=num.knots)



#### BUGS Fitting ####

age.income.fit<- bugs(data = data, inits = inits, parameters.to.save = parameters, model.file = "model.cp.txt",
                 n.chains = 1, n.iter = 100000, n.burnin = 10000, n.thin = 1,
                 debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE)


age.income.fit$DIC

age.income.fit$summary

plot(covariate, response)

# points(x=covariate, y=age.income.fit$median$beta[1] + age.income.fit$median$beta[2]*covariate + Z %*% age.income.fit$median$b, col="red")

age.income.prediction <- cbind(covariate,age.income.fit$median$ystar)

lines(age.income.prediction, col = "red")

#### BUGS Diagnostic ####

# str(age.income.fit)
plot(age.income.fit$sims.array[,1,"beta[1]"], type = "l", main = "Traceplot of Beta[1] simulated values")
plot(age.income.fit$sims.array[,1,"beta[2]"], type = "l")
# plot(age.income.fit$sims.list$beta[,1], type = "l")
# plot(age.income.fit$sims.list$beta[,2], type = "l")
acf(age.income.fit$sims.array[,,"beta[1]"])
acf(age.income.fit$sims.array[,,"beta[2]"])


##### SECOND BUGS FIT WITH DIFFERENT STARTING POINTS #######

inits.b2<-rep(10,num.knots)
inits2<-function(){
  list("beta"=c(10,10),"b"=inits.b2,"taub"=10,"taueps"=10)}


age.income.fit2<- bugs(data = data, inits = inits2, parameters.to.save = parameters, model.file = "model.cp.txt",
                      n.chains = 1, n.iter = 10000, n.burnin = 0, n.thin = 1,
                      debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE)


age.income.fit2$DIC

age.income.fit2$summary

plot(covariate, response, main = "Age-Income Predictor (Bugs, Second Simulation)")

# points(x=covariate, y=age.income.fit$median$beta[1] + age.income.fit$median$beta[2]*covariate + Z %*% age.income.fit$median$b, col="red")

age.income.prediction2 <- cbind(covariate,age.income.fit2$median$ystar)

lines(age.income.prediction2, col = "red")

#### BUGS 2 Diagnostic ####

# str(age.income.fit)
# head(names(age.income.fit$sims.array[1,1,]), 30)
plot(age.income.fit$sims.array[,1,"beta[1]"], type = "l")
plot(age.income.fit$sims.array[,1,"beta[2]"], type = "l")
# plot(age.income.fit$sims.list$beta[,1], type = "l")
# plot(age.income.fit$sims.list$beta[,2], type = "l")
acf(age.income.fit$sims.array[,,"beta[1]"])
acf(age.income.fit$sims.array[,,"beta[2]"])


#### JAGS #####
parameters.jags <- c("lambda","sigmab","sigmaeps","beta","b","ystar")
age.income.jags <- jags(data, inits, parameters.jags, model.file = "model.cp.txt",
                       n.chains = 1, n.iter = 100000, n.burnin = 10000, n.thin = 1)

plot(age.income.jags$BUGSoutput$sims.array[,1,"beta[1]"], type = "l")
acf(age.income.jags$BUGSoutput$sims.array[,1,"beta[1]"])
plot(age.income.jags$BUGSoutput$sims.array[,1,"beta[2]"], type = "l")
acf(age.income.jags$BUGSoutput$sims.array[,1,"beta[2]"])
plot(age.income.jags$BUGSoutput$sims.array[,1,"b[1]"], type = "l")
acf(age.income.jags$BUGSoutput$sims.array[,1,"b[1]"])

age.income.jags$BUGSoutput$DIC


# age.income.fit$median
# abs((age.income.jags$BUGSoutput$median$beta - age.income.fit$median$beta) / age.income.fit$median$beta)

age.income.jags$BUGSoutput$summary

# age.income.jags.pred <- age.income.jags$BUGSoutput$summary[1:22, c(3,5,7)]

plot(covariate, response)

age.income.prediction.jags <- cbind(covariate,age.income.jags$BUGSoutput$median$ystar)

lines(age.income.prediction.jags, col = "blue", lty = 2)

# age.income.prediction.jags <- age.income.jags$BUGSoutput$summary[27:231, c(3,5,7)]
# 
# age.income.pred_low <- cbind(covariate, age.income.prediction.jags[,1])
# age.income.pred_med <- cbind(covariate, age.income.prediction.jags[,2])
# age.income.pred_up <- cbind(covariate, age.income.prediction.jags[,3])
# 
# lines(age.income.pred_low, col = "green")
# lines(age.income.pred_med, col = "red")
# lines(age.income.pred_up, col = "green")
# 
# age.income.jags.pred <- age.income.jags$BUGSoutput$summary[27:231, c(3,5,7)]
# 
# age.income.jags.pred_low <- cbind(covariate, age.income.jags.pred[,1])
# age.income.jags.pred_med <- cbind(covariate, age.income.jags.pred[,2])
# age.income.jags.pred_up <- cbind(covariate, age.income.jags.pred[,3])
# 
# plot(covariate, response)
# 
# lines(age.income.jags.pred_low, col = "green")
# lines(age.income.jags.pred_med, col = "red")
# lines(age.income.jags.pred_up, col = "green")

age.income.jags2 <- jags(data, inits2, parameters.jags, model.file = "model.cp.txt",
                        n.chains = 1, n.iter = 100000, n.burnin = 0, n.thin = 1)
plot(age.income.jags2$BUGSoutput$sims.array[,1,"beta[1]"], type = "l", main = "Traceplot of Beta[1] simulated values")
acf(age.income.jags2$BUGSoutput$sims.array[,1,"beta[1]"])
plot(age.income.jags2$BUGSoutput$sims.array[,1,"beta[2]"], type = "l")
acf(age.income.jags2$BUGSoutput$sims.array[,1,"beta[2]"])
plot(age.income.jags2$BUGSoutput$sims.array[,1,"b[1]"], type = "l")
acf(age.income.jags2$BUGSoutput$sims.array[,1,"b[1]"])

age.income.jags2$BUGSoutput$DIC


# age.income.fit$median
# abs((age.income.jags$BUGSoutput$median$beta - age.income.fit$median$beta) / age.income.fit$median$beta)

age.income.jags2$BUGSoutput$summary

# age.income.jags.pred <- age.income.jags$BUGSoutput$summary[1:22, c(3,5,7)]

plot(covariate, response, main = "Age-Income Predictor (Jags, Second Simulation)")

age.income.prediction.jags2 <- cbind(covariate,age.income.jags2$BUGSoutput$median$ystar)

lines(age.income.prediction.jags2, col = "red")
