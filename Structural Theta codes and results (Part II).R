
# Structural Theta: codes and Results (Part II)"


## The R codes

# This files provides the codes that allow to reproduce results as in Table 7 to 9 of Sbrana and Silvestrini (IJF, 2024).
# 
# If you have any question please email to : giacomo.sbrana@neoma-bs.fr
# 
## Table 7: M4 competition: Prediction interval results 
# 
# Below we report the code that allows to reproduce the results as in Table 7 using the Structural THETA only. 
# 
# You may also change method and obtain the results for your choice.

library("M4comp2018")

####################### Yearly data ##################################################

Mdata<-Filter(function(l) l$period=="Yearly",M4)
In<-Out<-list();for(i in 1:length(Mdata)){In[[i]]=Mdata[[i]]$x;Out[[i]]=Mdata[[i]]$xx}
replic<-length(In)
steps<-6
freq<-frequ<-1

####################### Quarterly data ###############################################

# Mdata<-Filter(function(l) l$period=="Quarterly",M4)
# In<-Out<-list();for(i in 1:length(Mdata)){In[[i]]=Mdata[[i]]$x;Out[[i]]=Mdata[[i]]$xx}
# replic<-length(In)
# steps<-fh<-8
# freq<-frequ<-frq<-4

##########################################
THETAstructural<-function(y,steps){
  
  if(freq!=1){
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
  
  if(freq!=1){PointF=PointF*sfactout}
  
  
  H<-sigmae/(length(y)-1);P<-tail(p,1)*H;Eta<-q*H;
  Inter<-c();Inter[1]=P;
  for(j in 2:steps){Inter[j]=Inter[j-1]+Eta};
  Interv<-c();
  for(j in 1:steps){Interv[j]=Inter[j]+H};
  prob1=.85;prob2=.95; 
  lower1<-PointF-qnorm((1+prob1)/2)*sqrt(Interv);
  for(d in 1:steps){if(lower1[d]<0){lower1[d]=0}};
  upper1<-PointF+qnorm((1+prob1)/2)*sqrt(Interv);
  lower2<-PointF-qnorm((1+prob2)/2)*sqrt(Interv);
  for(d in 1:steps){if(lower2[d]<0){lower2[d]=0}};
  upper2<-PointF+qnorm((1+prob2)/2)*sqrt(Interv);
  
  list(mean=PointF,lower=cbind(lower1,lower2),upper=cbind(upper1,upper2))
  
}

#########################################


msis<-function(yout,lower,upper,prob,yin){
  a<-b<-d<-rep(0,length(yout));
  for(t in 1:length(yout)){
    if(yout[t]<lower[t]){a[t]=(2/(1-prob))*(lower[t]-yout[t])}
    if(yout[t]>upper[t]){b[t]=(2/(1-prob))*(yout[t]-upper[t])}
    d[t]=upper[t]-lower[t]}
  mean(a+b+d)/mean(abs(diff(yin,freq)))
}


resultsMIS95<-matrix(NA,length(In),1)

for(h in 1:length(In)){
  
  y=ts(c(In[[h]],Out[[h]]),frequency = freq);  
  forec=THETAstructural(head(y,(length(y)-steps)),steps)
  resultsMIS95[h,]<-msis(tail(y,steps),forec$lower[,2],forec$upper[,2],.95,head(y,length(y)-steps));
  
}


print(colMeans(resultsMIS95,na.rm = TRUE))




## Table 8: M3 competition: SMAPE results 

# Below we report the code that allows to reproduce the results as in Table 8 using the Structural THETA. 

# You may also change method and obtain the results for it.

library("Mcomp")

################################### MSOE (Structural) THETA #############################

## Warning : This code provides only point forecasts (not prediction intervals) #########
######## This is done since the M3 was focusing only on point forecasts #################
###### The full code can be found using this link: https://rpubs.com/giac76/1167511 #####

THETAstructural<-function(y,steps,freq){
  
  if(freq!=1){
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
  
  if(freq!=1){PointF=PointF*sfactout}
  
  
  H<-sigmae/(length(y)-1);P<-tail(p,1)*H;Eta<-q*H;
  Inter<-c();Inter[1]=P;
  for(j in 2:steps){Inter[j]=Inter[j-1]+Eta};
  Interv<-c();
  for(j in 1:steps){Interv[j]=Inter[j]+H};
  prob1=.85;prob2=.95; 
  lower1<-PointF-qnorm((1+prob1)/2)*sqrt(Interv);
  for(d in 1:steps){if(lower1[d]<0){lower1[d]=0}};
  upper1<-PointF+qnorm((1+prob1)/2)*sqrt(Interv);
  lower2<-PointF-qnorm((1+prob2)/2)*sqrt(Interv);
  for(d in 1:steps){if(lower2[d]<0){lower2[d]=0}};
  upper2<-PointF+qnorm((1+prob2)/2)*sqrt(Interv);
  
  PointF
  
}

######################### For M3 Data  ##########################

M3y<-subset(M3,"Yearly")
M3q<-subset(M3,"Quarterly")
M3m<-subset(M3,"Monthly")
M3o<-subset(M3,"Other")
In<-Out<-list();
for(i in 1:length(M3y)){In[[i]]=M3y[[i]]$x;Out[[i]]=M3y[[i]]$xx}
for(i in 1:length(M3q)){In[[length(M3y)+i]]=M3q[[i]]$x;Out[[length(M3y)+i]]=M3q[[i]]$xx}
for(i in 1:length(M3m)){In[[length(M3y)+length(M3q)+i]]=M3m[[i]]$x;Out[[length(M3y)+length(M3q)+i]]=M3m[[i]]$xx}
for(i in 1:length(M3o)){In[[length(M3y)+length(M3q)+length(M3m)+i]]=M3o[[i]]$x;Out[[length(M3y)+length(M3q)+length(M3m)+i]]=M3o[[i]]$xx}

# This is to chech that Naive 
#naiV=function(y,step, freq){naive(y,steps)$mean}

smape1<-function(DataOut,Forc,h){200*abs(DataOut[h]-Forc[h])/(DataOut[h]+Forc[h])}

################# Yearly data ###########################

steps<-6
freq<-1

SMAPE1y<-matrix(NA,length(M3y),6)
for(j in 1:length(M3y)){
  FO=THETAstructural(In[[j]],steps,freq)
  for(w in 1:steps){
    SMAPE1y[j,w]=smape1(Out[[j]][1:w],FO[1:w],w)
  }
}

################# Quarterly data###########################

steps<-8
freq<-4

SMAPE1q<-matrix(NA,length(M3q),steps)
for(j in 1:length(M3q)){
  FO=THETAstructural(In[[length(M3y)+j]],steps,freq)
  for(w in 1:steps){
    SMAPE1q[j,w]=smape1(Out[[length(M3y)+j]][1:w],FO[1:w],w)
  }
}

################# Monthly data###########################

steps<-18
freq<-12


SMAPE1m<-matrix(NA,length(M3m),steps)
for(j in 1:length(M3m)){
  FO=THETAstructural(In[[length(M3y)+length(M3q)+j]],steps,freq)
  for(w in 1:steps){
    SMAPE1m[j,w]=smape1(Out[[length(M3y)+length(M3q)+j]][1:w],FO[1:w],w)
  }
}

################# Other data###########################

steps<-8
freq<-1

SMAPE1o<-matrix(NA,length(M3o),steps)
for(j in 1:length(M3o)){
  FO=THETAstructural(In[[length(M3y)+length(M3q)+length(M3m)+j]],steps,freq)
  for(w in 1:steps){
    SMAPE1o[j,w]=smape1(Out[[length(M3y)+length(M3q)+length(M3m)+j]][1:w],FO[1:w],w)
  }
}

################# All data###########################

SMAPE<-matrix(NA,length(In),18)
for(j in 1:length(In)){
  if(j<=length(M3y)){steps=6;freq=1}
  if(j>length(M3y)&j<=(length(M3y)+length(M3q))){steps=8;freq=4}
  if(j>(length(M3y)+length(M3q))&j<=(length(M3y)+length(M3q)+length(M3m))){steps=18;freq=12}
  if(j>(length(M3y)+length(M3q)+length(M3m))){steps=8;freq=1}
  FO=THETAstructural(In[[j]],steps,freq)
  for(w in 1:steps){
    SMAPE[j,w]=smape1(Out[[j]][1:w],FO[1:w],w)
  }
}

resAve=c();j=1
for(i in c(4,6,8,12,15,18)){resAve[j]=mean(SMAPE[,1:i],na.rm = T);j=j+1};

round(c(mean(SMAPE1y[,1:4]),mean(SMAPE1y[,1:6])),2)
round(c(mean(SMAPE1q[,1:4]),mean(SMAPE1q[,1:6]),mean(SMAPE1q[,1:8])),2)
round(c(mean(SMAPE1m[,1:4]),mean(SMAPE1m[,1:6]),mean(SMAPE1m[,1:8]),mean(SMAPE1m[,1:12]),mean(SMAPE1m[,1:15]),mean(SMAPE1m[,1:18])),2)
round(c(mean(SMAPE1o[,1:4]),mean(SMAPE1o[,1:6]),mean(SMAPE1o[,1:8])),2)
print(round(resAve,2))



## Table 9: M5 competition: 

# Below we report the code that allows to reproduce the results as in Table 9 using the Structural THETA. 

# You may also change method and obtain the results for your choice.

aa1=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_insample1.txt"),header = FALSE))
aa2=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_insample2.txt"),header = FALSE))
aa3=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_insample3.txt"),header = FALSE))
aa4=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_insample4.txt"),header = FALSE))
aa5=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_insample5.txt"),header = FALSE))
outsample=as.matrix(read.table(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/M5_outsample.txt"),header = FALSE))
W=read.csv(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/WeightsM5.csv"),header = TRUE)
sales=rbind(aa1,aa2,aa3,aa4,aa5)
who=read.csv(url("https://raw.githubusercontent.com/giacsbrana/M5_dataset/main/ItemDeptCatStoreState.csv"))
steps=28

StrucThetaM5=function(y,steps){prob1=.5;prob2=.67;prob3=.95;prob4=.99;freq=7;

s=freq;h=steps;w<-rep(1/s,s);cma<-matrix(NA,length(y),1);
for(g in 4:(length(y)-3)){cma[g]<-sum(w*y[(g-3):(g+3)])};
residuals<-y-cma;sfactors<-c();for(seas in 1:s){
  sfactors[seas]<-mean(na.omit(residuals[seq(seas,length(y)-s+seas,by=s)]))}
sfactout<-rep(sfactors,length(y)+h)[(length(y)+1):(length(y)+h)]
y<-y-rep(sfactors,ceiling(length(y)/s))[1:length(y)]

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

PointF=PointF+sfactout


H=sigmae/(length(y)-1);
Inter<-c();for(j in 1:steps){Inter[j]=H*(((2+q+sqrt(q^2+4*q))/2)+q*(j-1))}
fo=PointF;Interv=Inter;

lower0<-upper0<-fo;lower50<-fo-qnorm((1+prob1)/2)*sqrt(Interv);
lower67<-fo-qnorm((1+prob2)/2)*sqrt(Interv);lower95<-fo-qnorm((1+prob3)/2)*sqrt(Interv);
lower99<-fo-qnorm((1+prob4)/2)*sqrt(Interv);upper50<-fo+qnorm((1+prob1)/2)*sqrt(Interv);
upper67<-fo+qnorm((1+prob2)/2)*sqrt(Interv);upper95<-fo+qnorm((1+prob3)/2)*sqrt(Interv);
upper99<-fo+qnorm((1+prob4)/2)*sqrt(Interv);list(mean=fo,lower=cbind(lower0,lower50,lower67,lower95,lower99),
                                                 upper=cbind(upper0,upper50,upper67,upper95,upper99)) 
}


steps=28

spl=function(y,act,metodo,steps){#y=ts(y,frequency = 7)
  #n=forecast(auto.arima(ts(y,frequency = 7)),h=steps,level = c(0,50,67,95,99))
  #n=naive(y,steps,level = c(0,50,67,95,99))
  n=metodo
  spl<-matrix(NA,9,steps)
  u=c(0.75,0.835,0.975,0.995)
  for(g in 1:4){
    for(s in 1:steps){
      if(n$upper[s,g+1]<= act[s]){spl[g,s]=u[g]*(act[s]-n$upper[s,g+1])}else{spl[g,s]=(1-u[g])*(n$upper[s,g+1]-act[s])}
    }
  }
  u=c(0.25,0.165,0.025,0.005)
  for(g in 1:4){
    for(s in 1:steps){n$lower[,g+1][n$lower[,g+1]<0]=0
    if(n$lower[s,g+1]<= act[s]){spl[g+4,s]=u[g]*(act[s]-n$lower[s,g+1])}else{spl[g+4,s]=(1-u[g])*(n$lower[s,g+1]-act[s])}
    }
  }
  for(s in 1:steps){
    if(n$upper[s,1]<= act[s]){spl[9,s]=.5*(act[s]-n$upper[s,1])}else{spl[9,s]=(1-.5)*(n$upper[s,1]-act[s])}
  }
  mean(rowMeans(spl)/mean(abs(diff(y))))
}


sale <- sales

stBU=matrix(0,nrow(sale),steps)


####################### Level 12 ###################

WRMSSE12<-0
SPL12<-0
for(good in 1:nrow(sale)){
  tims=sale[good,]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  act=outsample[good,]
  st=StrucThetaM5(y,steps)
  stBU[good,]=st$mean
  metodo=as.numeric(stBU[good,])
  err<-act-metodo
  SPL12<-SPL12+spl(y,act,st,steps)*W$weight[good]
  WRMSSE12<-WRMSSE12+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[good]
  
}

c(SPL12,WRMSSE12)

########################################## LEVEL 11 #######################



st11out<-list()
Level11out<-list()
Level11<-list()
ITEMS=names(table(who$item_id))
STATE=names(table(who$state_id ))

h=1
for(i in STATE){
  for(j in ITEMS){
    Level11[[h]]=as.numeric(colSums(sales[who$state_id==i&who$item_id==j,]))
    st11out[[h]]=as.numeric(colSums(stBU[who$state_id==i&who$item_id==j,]))
    Level11out[[h]]=as.numeric(colSums(outsample[who$state_id==i&who$item_id==j,]))
    
    h=h+1}
  
}


WRMSSE11td<-WRMSSE11bu<-SPL11<-0


for(good in 1:length(Level11)){
  tims=Level11[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)
  err<-Level11out[[good]]-st11out[[good]]
  errtd<-Level11out[[good]]-st$mean
  
  SPL11<-SPL11+spl(y,Level11out[[good]],st,steps)*W$weight[30490+good]
  WRMSSE11bu<-WRMSSE11bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+good]
  WRMSSE11td<-WRMSSE11td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+good]
  
}

c(SPL11,WRMSSE11bu,WRMSSE11td)


######################################### LEVEL 10 #######################

st10out<-list()
Level10out<-list()
Level10<-list()
ITEMS=names(table(who$item_id))

h=1
for(j in ITEMS){
  Level10[[h]]=as.numeric(colSums(sales[who$item_id==j,]))
  st10out[[h]]=as.numeric(colSums(stBU[who$item_id==j,]))
  Level10out[[h]]=as.numeric(colSums(outsample[who$item_id==j,]))
  h=h+1}


WRMSSE10td<-WRMSSE10bu<-SPL10<-0

for(good in 1:length(Level10)){
  tims=Level10[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)
  err<-Level10out[[good]]-st10out[[good]]#st11out[[good]]
  errtd<-Level10out[[good]]-st$mean
  SPL10<-SPL10+spl(y,Level10out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+good]
  WRMSSE10bu<-WRMSSE10bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+good]
  WRMSSE10td<-WRMSSE10td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+good]
}

c(WRMSSE10td,WRMSSE10bu,SPL10)

########################################## LEVEL 9 #######################

st9out<-list()
Level9out<-list()
Level9<-list()
ITEMS=names(table(who$dept_id))
STATE=names(table(who$store_id))

h=1
for(i in STATE){
  for(j in ITEMS){
    Level9[[h]]=as.numeric(colSums(sales[who$store_id==i&who$dept_id==j,]))
    Level9out[[h]]=as.numeric(colSums(outsample[who$store_id==i&who$dept_id==j,]))
    st9out[[h]]=as.numeric(colSums(stBU[who$store_id==i&who$dept_id==j,]))
    h=h+1}
}

WRMSSE9td<-WRMSSE9bu<-SPL9<-0

for(good in 1:length(Level9)){
  tims=Level9[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)
  
  err<-Level9out[[good]]-st9out[[good]]
  errtd<-Level9out[[good]]-st$mean
  SPL9<-SPL9+spl(y,Level9out[[good]],st,steps)*W$weight[30490+9147+3049+good]
  
  WRMSSE9bu<-WRMSSE9bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+good]
  WRMSSE9td<-WRMSSE9td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+good]
  
}

c(WRMSSE9td,WRMSSE9bu,SPL9)

########### LEVEL 8 #################################

st8out<-list()
Level8<-list()
Level8out<-list()
ITEMS=names(table(who$cat_id))
STATE=names(table(who$store_id))

h=1
for(i in STATE){
  for(j in ITEMS){
    Level8[[h]]=as.numeric(colSums(sales[who$store_id==i&who$cat_id==j,]))
    Level8out[[h]]=as.numeric(colSums(outsample[who$store_id==i&who$cat_id==j,]))
    st8out[[h]]=as.numeric(colSums(stBU[who$store_id==i&who$cat_id==j,]))
    h=h+1}
  
}

WRMSSE8td<-WRMSSE8bu<-SPL8<-0


for(good in 1:length(Level8)){
  
  tims=Level8[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)
  err<-Level8out[[good]]-st8out[[good]]
  errtd<-Level8out[[good]]-st$mean
  WRMSSE8bu<-WRMSSE8bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+good]
  WRMSSE8td<-WRMSSE8td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+good]
  SPL8<-SPL8+spl(y,Level8out[[good]],st,steps)*W$weight[30490+9147+3049+70+good]
}


######################################## Level 7 ######################

st7out<-list()
Level7<-list()
Level7out<-list()
ITEMS=names(table(who$dept_id ))
STATE=names(table(who$state_id))

h=1
for(i in STATE){
  for(j in ITEMS){
    Level7[[h]]=as.numeric(colSums(sales[who$state_id==i&who$dept_id==j,]))
    Level7out[[h]]=as.numeric(colSums(outsample[who$state_id==i&who$dept_id==j,]))
    st7out[[h]]=as.numeric(colSums(stBU[who$state_id==i&who$dept_id==j,]))
    h=h+1}
  
}


WRMSSE7td<-WRMSSE7bu<-SPL7<-0

for(good in 1:length(Level7)){
  tims=Level7[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)
  
  err<-Level7out[[good]]-st7out[[good]]
  errtd<-Level7out[[good]]-st$mean
  
  WRMSSE7bu<-WRMSSE7bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+good]
  WRMSSE7td<-WRMSSE7td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+good]
  SPL7<-SPL7+spl(y,Level7out[[good]],st,steps)*W$weight[30490+9147+3049+70+30+good]
}

c(WRMSSE7td,WRMSSE7bu,SPL7)

###################################  LEvel 6 ######################

st6out<-list()
Level6<-list()
Level6out<-list()
ITEMS=names(table(who$cat_id ))
STATE=names(table(who$state_id))

h=1
for(i in STATE){
  for(j in ITEMS){
    Level6[[h]]=as.numeric(colSums(sales[who$state_id==i&who$cat_id==j,]))
    Level6out[[h]]=as.numeric(colSums(outsample[who$state_id==i&who$cat_id==j,]))
    st6out[[h]]=as.numeric(colSums(stBU[who$state_id==i&who$cat_id==j,]))
    h=h+1}
  
}

WRMSSE6td<-WRMSSE6bu<-SPL6<-0

for(good in 1:length(Level6)){
  tims=Level6[[good]]
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)  
  err<-Level6out[[good]]-st6out[[good]]
  errtd<-Level6out[[good]]-st$mean
  
  SPL6<-SPL6+spl(y,Level6out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+good]  
  WRMSSE6td<-WRMSSE6td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+good]
  WRMSSE6bu<-WRMSSE6bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+good]
  
}

c(WRMSSE6bu,WRMSSE6td,SPL6)

#######################################   Level 5


Level5<-st5out<-list()
Level5out<-list()
ITEMS=names(table(who$cat_id ))
STATE=names(table(who$dept_id))

h=1
for(i in STATE){
  Level5[[h]]=as.numeric(colSums(sales[who$dept_id==i,]))
  Level5out[[h]]=as.numeric(colSums(outsample[who$dept_id==i,]))
  st5out[[h]]=as.numeric(colSums(stBU[who$dept_id==i,]))
  h=h+1
  
}

WRMSSE5td<-WRMSSE5bu<-SPL5<-0

for(good in 1:length(Level5)){
  tims=ts(Level5[[good]],frequency = 7)
  
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)  
  err<-Level5out[[good]]-st5out[[good]]
  errtd<-Level5out[[good]]-st$mean
  
  SPL5<-SPL5+spl(y,Level5out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+9+good]
  WRMSSE5bu<-WRMSSE5bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+good]
  WRMSSE5td<-WRMSSE5td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+good]
  
}

##################################### LEVEL 4 #################################
st4out<-list()
Level4<-list()
Level4out<-list()
ITEMS=names(table(who$cat_id))

h=1
for(j in ITEMS){
  Level4[[h]]=as.numeric(colSums(sales[who$cat_id==j,]))
  Level4out[[h]]=as.numeric(colSums(outsample[who$cat_id==j,]))
  st4out[[h]]=as.numeric(colSums(stBU[who$cat_id==j,]))
  h=h+1;
}

WRMSSE4td<-WRMSSE4bu<-SPL4<-0

for(good in 1:length(Level4)){
  tims=Level4[[good]]
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)  
  err<-Level4out[[good]]-st4out[[good]]
  errtd<-Level4out[[good]]-st$mean
  
  WRMSSE4bu<-WRMSSE4bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+good]
  WRMSSE4td<-WRMSSE4td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+good]
  SPL4<-SPL4+spl(y,Level4out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+9+7+good]
}

##################################### LEVEL 3 #################################
st3out<-list()
Level3<-list()
Level3out<-list()
STORE=names(table(who$store_id))

h=1
for(j in STORE){
  Level3[[h]]=as.numeric(colSums(sales[who$store_id==j,]))
  Level3out[[h]]=as.numeric(colSums(outsample[who$store_id==j,]))
  st3out[[h]]=as.numeric(colSums(stBU[who$store_id==j,]))
  h=h+1
  
}

WRMSSE3bu<-WRMSSE3td<-SPL3<-0

for(good in 1:length(Level3)){
  
  tims=Level3[[good]]
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  
  st=StrucThetaM5(y,steps)  
  err<-Level3out[[good]]-st3out[[good]]
  errtd<-Level3out[[good]]-st$mean
  
  WRMSSE3bu<-WRMSSE3bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+good]
  WRMSSE3td<-WRMSSE3td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+good]
  SPL3<-SPL3+spl(y,Level3out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+9+7+3+good]
}


##################################### LEVEL 2 #################################
st2out<-list()
Level2<-list()
Level2out<-list()
STATE=names(table(who$state_id))

h=1
for(j in STATE){
  Level2[[h]]=as.numeric(colSums(sales[who$state_id==j,]))
  Level2out[[h]]=as.numeric(colSums(outsample[who$state_id==j,]))
  st2out[[h]]=as.numeric(colSums(stBU[who$state_id==j,]))
  h=h+1
  
}

WRMSSE2td<-WRMSSE2bu<-SPL2<-0

for(good in 1:length(Level2)){
  tims=Level2[[good]]
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  
  st=StrucThetaM5(y,steps)  
  err<-Level2out[[good]]-st2out[[good]]
  errtd<-Level2out[[good]]-st$mean
  
  SPL2<-SPL2+spl(y,Level2out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+9+7+3+10+good]
  WRMSSE2bu<-WRMSSE2bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+10+good]
  WRMSSE2td<-WRMSSE2td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+10+good]
  
}

##################################### LEVEL 1 #################################
st1out<-list()
Level1<-list()
Level1out<-list()


h=1

Level1[[h]]=as.numeric(colSums(sales[,]))
Level1out[[h]]=as.numeric(colSums(outsample[,]))
st1out[[h]]=as.numeric(colSums(stBU[,]))
h=h+1;

WRMSSE1bu<-WRMSSE1td<-SPL1<-0


for(good in 1:length(Level1)){
  tims=ts(Level1[[good]],frequency = 7)
  for(t in 1:length(tims)){if(tims[t]!=0){y=tims[t:length(tims)];break}}
  st=StrucThetaM5(y,steps)  
  err<-Level1out[[good]]-st1out[[good]]
  errtd<-Level1out[[good]]-st$mean
  
  SPL1<-SPL1+spl(y,Level1out[[good]],StrucThetaM5(y,steps),steps)*W$weight[30490+9147+3049+70+30+21+9+7+3+10+3+good]
  WRMSSE1bu<-WRMSSE1bu+sqrt(mean(err^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+10+3+good]
  WRMSSE1td<-WRMSSE1td+sqrt(mean(errtd^2)/mean(diff(y)^2))*W$weight[30490+9147+3049+70+30+21+9+7+3+10+3+good]
  
}

resultsbu=c(WRMSSE1bu,WRMSSE2bu,WRMSSE3bu,WRMSSE4bu,WRMSSE5bu,WRMSSE6bu,WRMSSE7bu,WRMSSE8bu,WRMSSE9bu,WRMSSE10bu,WRMSSE11bu,WRMSSE12)
resultstd=c(WRMSSE1td,WRMSSE2td,WRMSSE3td,WRMSSE4td,WRMSSE5td,WRMSSE6td,WRMSSE7td,WRMSSE8td,WRMSSE9td,WRMSSE10td,WRMSSE11td,WRMSSE12)


rbind(round(c(resultstd,mean(resultstd)),3),round(c(resultsbu,mean(resultsbu)),3))

round(c(SPL1,SPL2,SPL3,SPL4,SPL5,SPL6,SPL7,SPL8,SPL9,SPL10,SPL11,SPL12,
        mean(c(SPL1,SPL2,SPL3,SPL4,SPL5,SPL6,SPL7,SPL8,SPL9,SPL10,SPL11,SPL12))),3)
