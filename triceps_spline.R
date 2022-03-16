library(R2WinBUGS)
library(R2jags)


data("triceps", package = "MultiKink")


# tri.age.plot <- ggplot(triceps, aes(x=age, y=triceps)) +
#   geom_point(alpha=0.55, color="black") +
#   theme_minimal()
# tri.age.plot
# 
# plot(triceps$age, triceps$triceps)
# a <- lm(triceps$triceps~triceps$age)
# abline(a = 6.1972, b = 0.2158)
# 
# plot(triceps$age[triceps$age<=10], triceps$triceps[triceps$age<=10])
# lm(triceps$triceps[triceps$age<=10]~triceps$age[triceps$age<=10])
# abline(a = 8.5871, b = -0.2762, col = "red")
# 
# plot(triceps$age[triceps$age>10 & triceps$age<=20], triceps$triceps[triceps$age>10 & triceps$age<=20])
# lm(triceps$triceps[triceps$age>10 & triceps$age<=20]~triceps$age[triceps$age>10 & triceps$age<=20])
# abline(a = -1.3202, b = 0.6936, col = "green")
# 
# plot(triceps$age[triceps$age>20 & triceps$age<=30], triceps$triceps[triceps$age>20 & triceps$age<=30])
# lm(triceps$triceps[triceps$age>20 & triceps$age<=30]~triceps$age[triceps$age>20 & triceps$age<=30])
# abline(a = 9.1757, b = 0.1629, col = "blue")


#### Triceps data ####
# data("triceps")
# n<-dim(triceps)[1]
# covariate<-triceps$age
# response<-triceps$triceps
# X<-cbind(rep(1,n),covariate)

#### Triceps data sampling ####
n<-dim(triceps)[1]
set.seed(1974)
sample.tri <- sample(1:n, 200)     # <<<----- cambia qui la dimensione del campione
n<-length(sample.tri)
covariate<-triceps$age[sample.tri]
response<-triceps$triceps[sample.tri]
ord <- order(covariate)
covariate<-covariate[ord]
response<-response[ord]
X<-cbind(rep(1,200),covariate)


#### Inizializzazione ####

plot(covariate, response)

num.knots<-20
knots<-quantile(unique(covariate),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])
# knots<- seq(0, 50, length.out = num.knots+1)
# knots<-knots[-1]

inits.b<-rep(0,num.knots)
inits<-function(){
  list("beta"=c(0,0),"b"=inits.b,"taub"=0.01,"taueps"=0.01)}

# inits.b2<- rep(1,num.knots)                             # seconda inizializzazione (meno precisa)
# inits2<-function(){
#   list(beta=c(1,1),b=inits.b2,taub=0.01,taueps=0.01)}

# parameters<-list("lambda","sigmab","sigmaeps","beta","b","ystar")



Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

data<-list("response"=response,"X"=X,"Z"=Z,"n"=n,"num.knots"=num.knots)
parameters<-list("lambda","sigmab","sigmaeps","beta","b","ystar")



### Fitting ####
triceps.fit <- bugs(data = data, inits = inits, parameters.to.save = parameters, model.file = "model.cp.txt",
                 n.chains = 1, n.iter = 100000, n.burnin = 10000, n.thin = 1,
                 debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE)

triceps.fit

triceps.fit$DIC

triceps.fit$median

triceps.fit$summary

# points(x=covariate, y=triceps.fit$median$beta[1] + triceps.fit$median$beta[2]*covariate + Z %*% triceps.fit$median$b, col="red")
plot(covariate, response, main = "Triceps Predictor (20 knots, Bugs Function)")

lines(x=covariate, y=triceps.fit$median$ystar, col = "green")

# str(triceps.fit)
# head(names(triceps.fit$sims.array[1,1,]), 30)
plot(triceps.fit$sims.array[,1,"beta[1]"], type = "l")
plot(triceps.fit$sims.array[,1,"beta[2]"], type = "l")
# plot(triceps.fit$sims.list$beta[,1], type = "l")
# plot(triceps.fit$sims.list$beta[,2], type = "l")
acf(triceps.fit$sims.array[,1,"beta[1]"])
acf(triceps.fit$sims.array[,1,"beta[2]"])



### JAGS ####


parameters.jags <- c("lambda","sigmab","sigmaeps","beta","b","ystar")
triceps.fit.jags <- jags(data, inits, parameters.jags, model.file = "model.cp.txt",
                    n.chains=1, n.iter=100000, n.burnin=10000)


plot(triceps.fit.jags$BUGSoutput$sims.array[,1,"beta[1]"], type = "l", main = "Traceplot of Beta[1] simulated values (Jags)")
acf(triceps.fit.jags$BUGSoutput$sims.array[,1,"beta[1]"], main = "Autocorrelation of Beta[1] simulated values (Jags)")
plot(triceps.fit.jags$BUGSoutput$sims.array[,1,"beta[2]"], type = "l", main = "Traceplot of Beta[2] simulated values (Jags)")
acf(triceps.fit.jags$BUGSoutput$sims.array[,1,"beta[2]"], main = "Autocorrelation of Beta[1] simulated values (Jags)")
plot(triceps.fit.jags$BUGSoutput$sims.array[,1,"b[1]"], type = "l", main = "Traceplot of b[1] simulated values (Jags)")
acf(triceps.fit.jags$BUGSoutput$sims.array[,1,"b[1]"], main = "Autocorrelation of b[1] simulated values (Jags)")

triceps.fit.jags$BUGSoutput$DIC

triceps.jags.prediction <- cbind(covariate, triceps.fit.jags$BUGSoutput$median$ystar)

plot(covariate, response, main = "Triceps Predictor (Jags)")

lines(triceps.jags.prediction, col = "red")





#### MODELLO A 10 NODI #####


num.knots<-10
knots<-quantile(unique(covariate),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])

inits.b<-rep(0,num.knots)
inits<-function(){
  list("beta"=c(0,0),"b"=inits.b,"taub"=0.01,"taueps"=0.01)}



Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

data<-list("response"=response,"X"=X,"Z"=Z,"n"=n,"num.knots"=num.knots)



#### Fitting ####
triceps.fit.10 <- bugs(data = data, inits = inits, parameters.to.save = parameters, model.file = "model.10.txt",
                    n.chains = 1, n.iter = 100000, n.burnin = 10000, n.thin = 1,
                    debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE)

triceps.fit.10

triceps.fit.10$DIC

triceps.fit.10$median

triceps.fit.10$summary


plot(covariate, response, main = "Triceps Predictor (Bugs, 10 knots)")

lines(x=covariate, y=triceps.fit.10$median$ystar, col = "blue")


plot(triceps.fit.jags.10$sims.array[,1,"beta[1]"], type = "l")
plot(triceps.fit.jags.10$sims.array[,1,"beta[2]"], type = "l")
acf(triceps.fit.jags.10$sims.array[,1,"beta[1]"])
acf(triceps.fit.jags.10$sims.array[,1,"beta[2]"])
