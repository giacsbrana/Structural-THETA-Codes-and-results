
# "Structural Theta: codes and Results (Part I)"

## The R codes

# This files provides the codes that allow to reproduce all results as in Table 3 to 6 of Sbrana and Silvestrini (IJF, 2024).
# 
# If you have any question please email to : giacomo.sbrana@neoma-bs.fr
# 
# ## Table 3: Monte Carlo results 
# 
# Below we report the code that allows to reproduce the results as in Table 3 using the Structural THETA.
# 
# This code generates two sets of processes (DGP1 and DGP2). You can select the one you wish.

library("forecast")

THETAstructural<-function(y,steps){
  
  likelihood<-function(parameters){q=abs(parameters[1]);co=abs(parameters[2])
  state<-p<-rep(0,length(y)+1)
  state[1]=y[1]
  p[1]=10000
  k<-inn<-rep(0,length(y))
  k[1]=p[1]/(p[1]+1);
  sigmae=0
  
  p[2]=p[1]-p[1]*k[1]+q;
  state[2]=co+state[1]+k[1]*inn[1]
  
  for(t in (2):length(y)){
    
    k[t]=p[t]/(p[t]+1);
    inn[t]=y[t]-state[t]
    sigmae=sigmae+inn[t]^2/(p[t]+1)
    p[t+1]=p[t]-p[t]*k[t]+q;
    state[t+1]=co+state[t]+k[t]*inn[t]
  }
  
  sum(log(p+1))+length(y)*log(sigmae/length(y))}
  
  results<-optim(c(.01,.001),likelihood)
  
  q=abs(results[[1]][1]);co=abs(results[[1]][2])
  state<-p<-rep(0,length(y)+1)
  state[1]=y[1]
  p[1]=10000
  k<-inn<-rep(0,length(y))
  k[1]=p[1]/(p[1]+1);
  sigmae=0
  
  p[2]=p[1]-p[1]*k[1]+q;
  state[2]=co+state[1]+k[1]*inn[1]
  
  for(t in (2):length(y)){
    
    k[t]=p[t]/(p[t]+1);
    inn[t]=y[t]-state[t]
    sigmae=sigmae+inn[t]^2/(p[t]+1)
    p[t+1]=p[t]-p[t]*k[t]+q;
    state[t+1]=co+state[t]+k[t]*inn[t]
  }
  
  PointF<-rep(0,steps);PointF[1]=state[length(y)+1];
  for( t in 2:steps){PointF[t]=co*(t-1)+state[length(y)+1]}
  
  
  H<-sigmae/(length(y)-1);P<-tail(p,1)*H;Eta<-q*H;
  Inter<-c();Inter[1]=P;
  for(j in 2:steps){Inter[j]=Inter[j-1]+Eta};
  Interv<-c();
  for(j in 1:steps){Interv[j]=Inter[j]+H};
  prob1=.85;prob2=.95; 
  lower1<-PointF-qnorm((1+prob1)/2)*sqrt(Interv);
  #for(d in 1:steps){if(lower1[d]<0){lower1[d]=0}};
  upper1<-PointF+qnorm((1+prob1)/2)*sqrt(Interv);
  lower2<-PointF-qnorm((1+prob2)/2)*sqrt(Interv);
  #for(d in 1:steps){if(lower2[d]<0){lower2[d]=0}};
  upper2<-PointF+qnorm((1+prob2)/2)*sqrt(Interv);
  
  list(mean=PointF,lower=cbind(lower1,lower2),upper=cbind(upper1,upper2))
  
}


AUTO<-function(y,steps){forecast::forecast(auto.arima(y),h=steps,level = c(.85,.95))}
THETA<-function(y,steps){forecast::forecast(thetaf(y,h=steps,level = c(.85,.95)))}
ETS<-function(y,steps){forecast::forecast(ets(y),h=steps)}
competitors<-list(THETAstructural,ETS,AUTO)

####################### data Monte Carlo ##################################################

steps<-6
freq<-frequ<-1
times=10000

set.seed(1)

DGP1=list()
for(u in 1:times){
  n=60+steps
  eta=runif(1,.1,.5)*rnorm(n)
  mu=y=trend=c()
  alpha=runif(1,.1,5);beta=runif(1,.01,.1)
  ar=runif(1,.7,.99)
  delta=runif(1,-.7,-.1)
  trend[1]=alpha+beta
  mu[1]=eta[1];y[1]=trend[1]+mu[1]
  for(t in 2:n){mu[t]=ar*mu[t-1]+eta[t]+delta*eta[t-1]
  trend[t]=alpha+beta*t
  y[t]=trend[t]+mu[t] 
  }
  DGP1[[u]]=y
}


DGP2=list()
for(u in 1:times){
  n=60+steps
  eta=runif(1,.1,.5)*rnorm(n)
  mu=y=c()
  ar=runif(1,.1,.7)
  delta=runif(1,-.7,0)
  co=runif(1,.01,.1);
  trend=somma=c();trend[1]=co
  mu[1]=eta[1];y[1]=co;somma[1]=mu[1]
  for(t in 2:n){mu[t]=ar*mu[t-1]+eta[t]+delta*eta[t-1]
  trend[t]=trend[t-1]+co
  somma[t]=mu[t]+somma[t-1]
  y[t]=somma[t]+trend[t]
  
  }
  DGP2[[u]]=y
}


######################## Select the dataset to use! ##############################
DATA=DGP1
goods=length(DATA)

#############################################################################################

rmsse<-function(y,yf,steps){OutSample=tail(y,length(yf));sqrt(mean((OutSample[1:steps]-yf[1:steps])^2)/mean(diff(head(y,(length(y)-length(yf))),frequ)^2))}
mase<-function(y,yf,steps){OutSample=tail(y,length(yf));(mean(abs(OutSample[1:steps]-yf[1:steps])))/mean(abs(diff(head(y,(length(y)-length(yf))),frequ)))}


mis<-function(yout,lower,upper,prob,yin){
  a<-b<-d<-rep(0,length(yout));
  for(t in 1:length(yout)){
    if(yout[t]<lower[t]){a[t]=(2/(1-prob))*(lower[t]-yout[t])}
    if(yout[t]>upper[t]){b[t]=(2/(1-prob))*(yout[t]-upper[t])}
    d[t]=upper[t]-lower[t]}
  mean(a+b+d)/mean(abs(diff(yin,freq)))
}

par(mfrow=c(1,3))

resMASE<-resRMSSE<-resSMAPE<-matrix(NA,steps,length(competitors))
resCOVER85<-resCOVER95<-resMIS85<-resMIS95<-matrix(NA,1,length(competitors))
results1MASE<-resultsMASE<-results1RMSSE<-resultsRMSSE<-results1SMAPE<-resultsSMAPE<-forec<-list()
resultsCOVER85<-resultsMIS85<-resultsCOVER95<-resultsMIS95<-matrix(NA,goods,length(competitors))
results1COVER85<-results1MIS85<-results1COVER95<-results1MIS95<-list()

for(h in 1:length(DATA)){
  
  y=ts(DATA[[h]],frequency = freq);
  for(j in 1:length(competitors)){
    forec[[j]]<-competitors[[j]](head(y,(length(y)-steps)),steps)
    for(s in 1:steps){
      resRMSSE[s,j]<-rmsse(y,forec[[j]]$mean,s)
      resMASE[s,j]<-mase(y,forec[[j]]$mean,s) 
    }
    
    resMIS85[1,j]<-mis(tail(y,steps),forec[[j]]$lower[,1],forec[[j]]$upper[,1],.85,head(y,length(y)-steps));
    resMIS95[1,j]<-mis(tail(y,steps),forec[[j]]$lower[,2],forec[[j]]$upper[,2],.95,head(y,length(y)-steps));
    
  }
  results1RMSSE[[1]]=resRMSSE
  results1MASE[[1]]=resMASE
  
  resultsRMSSE[[h]]<-Reduce("+", results1RMSSE)
  resultsMASE[[h]]<-Reduce("+", results1MASE)
  resultsMIS85[h,]<-resMIS85
  resultsMIS95[h,]<-resMIS95
  data_frame <- data.frame(col1 = 1:steps)  
  
}

print(round(Reduce("+", resultsRMSSE)/h,3))
print(rbind(round(colMeans(resultsMIS85,na.rm = TRUE),3),round(colMeans(resultsMIS95,na.rm = TRUE),3)))


## Table 4, 5 and 6: The M4 competition (Point forecast)

#This code allows to reproduce the results as in Table 4,5 and 6 of Sbrana & Silvestrini (IJF, 2024).

#The results refer to Yearly data used in the M4 competition comparing:
  
#  The Structural (MSOE) THETA, SES, Holt, Damped, Comb, THETA classic, AutoTHETA, DOTM, RWDAR, naive, naive seasonal, naive2.

#One can test additional methods by adding them among the competitors as in the Benchmarks function below.

#One can change the frequency of data. For example, we provide the code to use All M4 data, which allows to reproduce Table 6.

#We used part of the code already available on github. However, we adapted it to the specific M4 competition. 


library("Mcomp")
library("M4comp2018")
library("forecTheta")

############################# All M4 data #####################################
##### This code allows testing the performance over all M4 data ###############
###############################################################################
# Mdata<-M4
# In<-Out<-list();for(i in 1:length(Mdata)){In[[i]]=Mdata[[i]]$x;Out[[i]]=Mdata[[i]]$xx}
# replic<-length(In)
# data_train = In
# data_test=Out

############################# Yearly data #####################################
##### This code allows testing the performance over as specific frequency #####
###############################################################################
############## You may change "Yearly" with any other frequency you wish ######

Mdata<-Filter(function(l) l$period=="Yearly",M4)
In<-Out<-list();for(i in 1:length(Mdata)){In[[i]]=Mdata[[i]]$x;Out[[i]]=Mdata[[i]]$xx}
replic<-length(In)
data_train = In
data_test=Out
fh=6
###### Dynamic Optimized Theta Model ##########################################
#(forecTheta library) 
######## https://www.sciencedirect.com/science/article/pii/S0169207016300243#!


DOTM<-function(y,steps){s=freq
dotm(ts(y,frequency = s),h=steps,level=c(80,95))
#dotm(y,h=steps)$mean
}


################################# Structural THETA code ###############################

THETAstructural=function(y,steps){frq=freq
prob1=.80;prob2=.95;
if(frq!=1){
  s=freq;
  for(t in 1:length(y)){if(y[t]!=0){y<-y[t:length(y)];break}}
  w<-rep(1/(2*s),s+1);w[2:s]<-1/s
  cma<-matrix(NA,length(y),1);
  for(g in 1:(length(y)-s)){cma[g+s/2]<-sum(w*y[g:(g+s)])};
  residuals<-y/cma
  sfactors<-c();for(seas in 1:s){
    sfactors[seas]<-mean(na.omit(residuals[seq(seas,length(y)-s+seas,by=s)]))}
  sfactors<-sfactors*s/sum(sfactors)
  sfactout<-rep(sfactors,length(y)+steps)[(length(y)+1):(length(y)+steps)]
  y<-y/rep(sfactors,ceiling(length(y)/s))[1:length(y)]}

likelihood<-function(parameters){q=abs(parameters[1]);co=abs(parameters[2])
state<-p<-rep(0,length(y)+1)
state[1]=y[1]
p[1]=10000
k<-inn<-rep(0,length(y))
k[1]=p[1]/(p[1]+1);
sigmae=0

p[2]=p[1]-p[1]*k[1]+q;
state[2]=co+state[1]+k[1]*inn[1]

for(t in (2):length(y)){
  
  k[t]=p[t]/(p[t]+1);
  inn[t]=y[t]-state[t]
  sigmae=sigmae+inn[t]^2/(p[t]+1)
  p[t+1]=p[t]-p[t]*k[t]+q;
  state[t+1]=co+state[t]+k[t]*inn[t]
}

sum(log(p+1))+length(y)*log(sigmae/length(y))}

results<-optim(c(.01,1),likelihood)

q=abs(results[[1]][1]);co=abs(results[[1]][2])
state<-p<-rep(0,length(y)+1)
state[1]=y[1]
p[1]=10000
k<-inn<-rep(0,length(y))
k[1]=p[1]/(p[1]+1);

p[2]=p[1]-p[1]*k[1]+q;
state[2]=co+state[1]+k[1]*inn[1]

for(t in (2):length(y)){
  
  k[t]=p[t]/(p[t]+1);
  inn[t]=y[t]-state[t]
  p[t+1]=p[t]-p[t]*k[t]+q;
  state[t+1]=co+state[t]+k[t]*inn[t]
}

PointF<-rep(0,steps);PointF[1]=state[length(y)+1];
for( t in 2:steps){PointF[t]=co*(t-1)+state[length(y)+1]}

if(frq!=1){PointF=PointF*sfactout}


PointF

}

# This code forecast with the RWDAR as in Sbrana & Silvestrini (IJF, 2023)

rwdar<-function(y,h){s=freq
if(s!=1){w<-rep(1/(2*s),s+1);w[2:s]<-1/s
cma<-matrix(NA,length(y),1);
for(g in 1:(length(y)-s)){cma[g+s/2]<-sum(w*y[g:(g+s)])};
residuals<-y/cma
sfactors<-c();for(seas in 1:s){
  sfactors[seas]<-mean(na.omit(residuals[seq(seas,length(y)-s+seas,by=s)]))}
sfactors<-sfactors*s/sum(sfactors)
sfactout<-rep(sfactors,length(y)+h)[(length(y)+1):(length(y)+h)]
y<-y/rep(sfactors,ceiling(length(y)/s))[1:length(y)]}else{sfactout=rep(1,h)}

obs<-length(y);  lambda <- beta <- c();   v <- rep(0,obs)
lambda[1]=y[1]
beta[1]<-0
ConcLogLikelihood <- function(Kg){q<-1-exp(-abs(Kg[1]));p<-1-exp(-abs(Kg[2]));co<-abs(Kg[3])
k1<- -(2*sqrt(q))/(sqrt(q)*(-1+p)-sqrt(4+q* (1+p)^2))
k2<-(p*(2+q+q*p-sqrt(q)*sqrt(4+q *(1+p)^2)))/(2+2 *q*p)
for (t in 2:(obs)) {
  v[t]<-y[t]-lambda[t-1]-beta[t-1]
  lambda[t]<-co+lambda[t-1]+k1*v[t]
  beta[t]<-p*beta[t-1]+k2*v[t]
}
sum(v^2)
}
result<-optim(c(.8,.4,2),ConcLogLikelihood)
q<-1-exp(-abs(result[[1]][[1]])); p<-1-exp(-abs(result[[1]][[2]]));co<-abs(result[[1]][[3]])
if(p==1){p=.99}
k1<- -(2*sqrt(q))/(sqrt(q)*(-1+p)-sqrt(4+q* (1+p)^2))
k2<-(p*(2+q+q*p-sqrt(q)*sqrt(4+q *(1+p)^2)))/(2+2 *q*p)
for (t in 2:(obs)) {
  v[t]<-y[t]-lambda[t-1]-beta[t-1]
  lambda[t]<-co+lambda[t-1]+k1*v[t]
  beta[t]<-p*beta[t-1]+k2*v[t]
}
T<-Tmatp<-diag(1,2);T[2,2]<-p;Z<-c(1,1)

Forecast<-c()
for(i in 1:h){
  Forecast[i]<-Z%*%Tmatp%*%rbind(lambda[obs],beta[obs]);
  Forecast[i]<-Forecast[i]+(i-1)*co
  Tmatp<-T%*%Tmatp
}
Forecast*sfactout
}

##################################################################################

SeasonalityTest <- function(input, ppy){
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )
    
    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }
  return(test_seasonal)
}

Theta.models.fit <- function(input, fh, theta, curve, model, seasonality , plot=FALSE, positive=TRUE){
  #Check if the inputs are valid
  if (theta<1){ theta <- 1}
  if (fh<1){ fh <- 1}
  #Estimate theta line weights
  outtest <- naive(input, h=fh)$mean
  wses <- (1/theta) ; wlrl <- (1-wses)
  #Estimate seasonaly adjusted time series
  ppy <- frequency(input)
  if (seasonality=="N"){
    des_input <- input ; SIout <- rep(1, fh) ; SIin <- rep(1, length(input))
  }else if (seasonality=="A"){
    Dec <- decompose(input, type="additive")
    des_input <- input-Dec$seasonal 
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    Dec <- decompose(input, type="multiplicative")
    des_input <- input/Dec$seasonal 
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }
  
  #Estimate theta line zero
  observations <- length(des_input)
  xs <- c(1:observations)
  xf = xff <- c((observations+1):(observations+fh))
  dat=data.frame(des_input=des_input, xs=xs)
  newdf <- data.frame(xs = xff)
  
  if (curve=="Exp"){
    estimate <- lm(log(des_input)~xs)
    thetaline0In <- exp(predict(estimate))+input-input
    thetaline0Out <- exp(predict(estimate, newdf))+outtest-outtest
  }else{
    estimate <- lm(des_input ~ poly(xs, 1, raw=TRUE))
    thetaline0In <- predict(estimate)+des_input-des_input
    thetaline0Out <- predict(estimate, newdf)+outtest-outtest
  }
  
  #Estimete Theta line (theta)
  if (model=="A"){
    thetalineT <- theta*des_input+(1-theta)*thetaline0In
  }else{
    thetalineT <- (des_input^theta)*(thetaline0In^(1-theta))
  }
  
  #forecasting TL2
  sesmodel <- ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  
  #Theta forecasts
  if (model=="A"){
    forecastsIn <- as.numeric(thetaline2In*wses)+as.numeric(thetaline0In*wlrl)+des_input-des_input
    forecastsOut <- as.numeric(thetaline2Out*wses)+as.numeric(thetaline0Out*wlrl)+outtest-outtest
  }else{
    forecastsIn <- ((as.numeric(thetaline2In)^(1/theta))*(as.numeric(thetaline0In)^(1-(1/theta))))+des_input-des_input
    forecastsOut <- ((as.numeric(thetaline2Out)^(1/theta))*(as.numeric(thetaline0Out)^(1-(1/theta))))+outtest-outtest
  }
  
  #Seasonal adjustments
  if (seasonality=="A"){
    forecastsIn <- forecastsIn+SIin
    forecastsOut <- forecastsOut+SIout
  }else{
    forecastsIn <- forecastsIn*SIin
    forecastsOut <- forecastsOut*SIout
  }
  
  #Zero forecasts become positive
  if (positive==T){
    for (i in 1:length(forecastsOut)){
      if (forecastsOut[i]<0){ forecastsOut[i] <- 0 }
    }
  }
  
  if (plot==TRUE){
    united <- cbind(input,forecastsOut)
    for (ik in 1:(observations+fh)){ united[ik,1] = sum(united[ik,2],united[ik,1], na.rm = TRUE) }
    plot(united[,1],col="black",type="l",main=paste("Model:",model,",Curve:",curve,",Theta:",theta),xlab="Time",ylab="Values",
         ylim=c(min(united[,1])*0.85,max(united[,1])*1.15))
    lines(forecastsIn, col="green") ; lines(forecastsOut, col="green")
    lines(thetaline2In, col="blue") ; lines(thetaline2Out, col="blue")
    lines(thetaline0In, col="red") ; lines(thetaline0Out, col="red")
  }
  
  output=list(fitted=forecastsIn,mean=forecastsOut,
              fitted0=thetaline0In,mean0=thetaline0Out,
              fitted2=thetaline2In,mean2=thetaline2Out,
              model=paste(seasonality,model,curve,
                          c(round(theta,2)),
                          round(sesmodel$model$par[1],3),
                          round(sesmodel$model$par[2],3)))
  
  return(output)
}

AutoTheta<- function(input, fh, positive=TRUE){
  
  if (min(input)>0){
    molist <- c("M","A") ; trlist <- c("Lrl","Exp")
  }else{
    molist <- c("A") ; trlist <- c("Lrl")
  }
  #Scale
  base <- mean(input) ; input <- input/base
  
  #Check seasonality & Create list of models
  ppy <- frequency(input) ; ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input, ppy) }
  if (ST==T){
    
    selist <- c("M","A")
    listnames <- c()
    for (i in 1:length(selist)){
      for (ii in 1:length(molist)){
        for (iii in 1:length(trlist)){
          listnames <- c(listnames,paste(selist[i], molist[ii], trlist[iii]))
        }
      }
    }
    
  }else{
    
    listnames <- c()
    for (ii in 1:length(molist)){
      for (iii in 1:length(trlist)){
        listnames <- c(listnames, paste("N", molist[ii], trlist[iii]))
      }
    }
    
  }
  #Exclude instable models
  excluded <- c("N M Lrl", "A M Lrl", "A M Exp", "M M Lrl")
  listnames <- listnames[!(listnames %in% excluded)]
  modellist <- NULL
  for (i in 1:length(listnames)){
    modellist[length(modellist)+1] <- list(c(substr(listnames,1,1)[i], substr(listnames,3,3)[i],
                                             substr(listnames,5,7)[i]))
  }
  
  #Start validation
  errorsin <- c() ; models <- NULL
  
  #With this function determine opt theta per case
  optfun <- function(x, input, fh, curve, model, seasonality){
    mean(abs(Theta.models.fit(input=input, fh, theta=x, curve, model, seasonality , plot=FALSE)$fitted-input))
  }
  
  for (j in 1:length(listnames)){
    optTheta <- suppressWarnings(optimize(optfun, c(1:3), 
                                          input=input, fh=fh, curve=modellist[[j]][3], model=modellist[[j]][2], 
                                          seasonality=modellist[[j]][1])$minimum)
    
    fortheta <- Theta.models.fit(input=input, fh=fh, theta=optTheta, curve=modellist[[j]][3], model=modellist[[j]][2], 
                                 seasonality=modellist[[j]][1], plot=F)
    models[length(models)+1] <- list(fortheta)
    errorsin <- c(errorsin, mean(abs(input-fortheta$fitted)))
  }
  
  #Select model and export
  selected.model <- models[[which.min(errorsin)]]
  description <- selected.model$model
  
  #Estimate Prediction Intervals
  frc <- selected.model$mean*base
  fitted <- selected.model$fitted*base
  residuals_t <- as.numeric(input*base-fitted)
  
  if (frequency(input)==1){
    m <- 12
  }else if (frequency(input)==4){
    m <- 4
  }else{
    m <- 1
  }
  
  pisl <- frc-1.960*sd(residuals_t)*sqrt(1+m*(c(1:fh)-1))
  pisu <- frc+1.960*sd(residuals_t)*sqrt(1+m*(c(1:fh)-1))
  if (positive==T){
    pisl[pisl<0] <- 0 ; pisu[pisu<0] <- 0
  }
  output <- list(fitted=fitted, mean=frc, description=description, piu=pisu, pil=pisl) 
  
  return(output)
}

#################################################################################

smape_cal <- function(outsample, forecasts){
  #Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts){
  #Used to estimate MASE
  frq <- frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
}

naive_seasonal <- function(input, fh){
  #Used to estimate Seasonal Naive
  frcy <- frequency(input)
  frcst <- naive(input, h=fh)$mean 
  if (frcy>1){ 
    frcst <- head(rep(as.numeric(tail(input,frcy)), fh), fh) + frcst - frcst
  }
  return(frcst)
}

Theta.classic <- function(input, fh){
  #Used to estimate Theta classic
  
  #Set parameters
  wses <- wlrl<-0.5 ; theta <- 2
  #Estimate theta line (0)
  observations <- length(input)
  xt <- c(1:observations)
  xf <- c((observations+1):(observations+fh))
  train <- data.frame(input=input, xt=xt)
  test <- data.frame(xt = xf)
  
  estimate <- lm(input ~ poly(xt, 1, raw=TRUE))
  thetaline0In <- as.numeric(predict(estimate))
  thetaline0Out <- as.numeric(predict(estimate,test))
  
  #Estimate theta line (2)
  thetalineT <- theta*input+(1-theta)*thetaline0In
  sesmodel <- ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  
  #Theta forecasts
  forecastsIn <- (thetaline2In*wses)+(thetaline0In*wlrl)
  forecastsOut <- (thetaline2Out*wses)+(thetaline0Out*wlrl)
  
  #Zero forecasts become positive
  for (i in 1:length(forecastsOut)){
    if (forecastsOut[i]<0){ forecastsOut[i]<-0 }
  }
  
  output=list(fitted = forecastsIn, mean = forecastsOut,
              fitted0 = thetaline0In, mean0 = thetaline0Out,
              fitted2 = thetaline2In, mean2 = thetaline2Out)
  
  return(output)
}


Benchmarks <- function(input, fh){
  #Used to estimate the statistical benchmarks of the M4 competition
  
  #Estimate seasonally adjusted time series
  ppy <- frequency(input) ; ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==TRUE){
    Dec <- decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    des_input <- input ; SIout <- rep(1, fh)
  }
  
  #f1 <- naive(input, h=fh)$mean #Naive
  #f2 <- naive_seasonal(input, fh=fh) #Seasonal Naive
  #f3 <- naive(des_input, h=fh)$mean*SIout #Naive2
  f4 <- ses(des_input, h=fh)$mean*SIout #Ses
  f5 <- holt(des_input, h=fh, damped=F)$mean*SIout #Holt
  f6 <- holt(des_input, h=fh, damped=T)$mean*SIout #Damped
  f7 <- Theta.classic(input=des_input, fh=fh)$mean*SIout #Theta
  f8 <- (f4+f5+f6)/3 #Comb
  f9 <- THETAstructural(input,fh)
  f10 <-DOTM(input,fh)$mean
  f11 <-AutoTheta(input,fh)$mean
  f12 <-rwdar(input,fh)
  
  return(list(f4,f5,f6,f7,f8,f9,f10,f11,f12))
  
}

Names_benchmarks <- c("SES", "Holt", "Damped", "Theta", "Comb", "MSOE Theta","DOTM","AutoTheta","RWDAR")

Total_smape=Total_mase <- array(NA,dim = c(length(Names_benchmarks), fh, length(data_train)))

# # PARALLEL
# system.time(expr =
#               forecasts <- mclapply(X=data_train, FUN = Benchmarks, fh=fh, mc.cores=20)
# )

forecasts <- lapply(X = data_train, FUN = Benchmarks, fh=fh)

for (i in 1:length(forecasts)) {
  insample <- data_train[[i]]
  outsample <- data_test[[i]]
  forecast = forecasts[[i]]
  
  #sMAPE
  for (j in 1:length(Names_benchmarks)){
    Total_smape[j,,i] <- smape_cal(outsample, forecast[[j]]) #j the # of the benchmark
  }
  #MASE
  for (j in 1:length(Names_benchmarks)){
    Total_mase[j,,i] <- mase_cal(insample, outsample, forecast[[j]]) #j the # of the benchmark
  }
  
}

print("########### sMAPE ###############")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i], round(mean(Total_smape[i,,]), 3)))
}
print("########### MASE ################")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i], round(mean(Total_mase[i,,]), 3)))
}

####

Total_mase_m=matrix(NA, length(Names_benchmarks), fh)
for (i in 1:length(Names_benchmarks)){
  for (j in 1:fh){
    Total_mase_m[i,j] = mean(Total_mase[i,j,])
  }
}
Total_mase_m=t(Total_mase_m)
colnames(Total_mase_m)=Names_benchmarks

overalResults <- matrix(c(colMeans(Total_mase_m),apply(Total_mase_m,2,median)),
                        ncol(Total_mase_m), 2, dimnames=list(colnames(Total_mase_m),c("Mean","Median")))
round(overalResults,5)

boxplot(Total_mase_m)
points(colMeans(Total_mase_m),col="red",pch=16)

#options(download.file.method="wininet")
#install.packages("tsutils")
#help(print.nemenyi)
library(tsutils)
nemenyi(Total_mase_m, plottype="mcb")