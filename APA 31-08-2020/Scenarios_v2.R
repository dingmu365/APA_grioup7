
# p = 11; # number of total covariates
# asym = .5 ;# whether treatment effects are distributed asymmetrically across treated and control
# n = 1000 ;# total size of the dataset
# propens = .5; #treatment probability
# sig = .01;
# treatsize = .5; # treatment effect size
# levsize = 1;
# mu_num;
# tau_num;
# option_num=1


Data_Gen<-function(p = 11, # number of total covariates
                   asym = .5 ,# whether treatment effects are distributed asymmetrically across treated and control
                   n = 1000 ,# total size of the dataset
                   propens = .5, #treatment probability
                   sig = .01,
                   treatsize = .5, # treatment effect size
                   levsize = 1,
                   mu_num,
                   tau_num,
                   option_num,
                   method="non_tr" # since the honesty causal tree need three samples
                   )
{

# draw W
w = rbinom(n, 1, propens)
# draw X
set.seed(123)
X <- as.data.frame( matrix(rnorm(n * p), nrow = n, ncol = p))
mu<-Data_gen(X,mu_num)

tau<-Data_gen(X,tau_num)
# generate outcomes as function of treatment status, mu, tau, and noise
y<-y_gen(mu,asym,tau,X,option_num=option_num,w,sig)
opt_num=0-option_num
y_<-y_gen(mu,asym,tau,X,option_num=opt_num,w,sig)
# create formulas for estimation
f <- ""
nextx <- ""
if (p>1) {
  for (ii in 1:(p-1)) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
    f <- paste(f, nextx, "+", sep="")
  }
  f <- paste(f, "x", ii+1, sep="")
} else if (p==1) {
  f <- "x1"
}

for (ii in 1:p) {
  nextx <- paste("x",ii, sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
}
nameall <- c( name,  "y", "w", "tau_true")

tau_true <- (1-2*w)*(y_ - y)  


data_gen<-as.list(NULL)
if(method=="honesty")
{
  
  ntr <- round(.333*n)
  nest <- round(.333*n)
  ntest <- n - ntr - nest
  dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
  dataEst <- data.frame(X[(ntr+1):(ntr+nest),], y[(ntr+1):(ntr+nest)], w[(ntr+1):(ntr+nest)], tau_true[(ntr+1):(ntr+nest)])
  dataTest <- data.frame(X[(ntr+nest+1):n,], y[(ntr+nest+1):n], w[(ntr+nest+1):n], tau_true[(ntr+nest+1):n])
  
  
  names(dataTrain)=nameall
  names(dataEst)=nameall
  names(dataTest)=nameall
  data_gen[["dataEst"]]<-dataEst
  
}

else if(method=="non_tr")
{
  #using 90% generated data as the training data, and the rest 10% as the test data
  
  ntr <- round(.9*n)
  ntest <- n - ntr 
dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], tau_true[1:ntr])
dataTest <- data.frame(X[(ntr+1):n,], y[(ntr+1):n], w[(ntr+1):n], tau_true[(ntr+1):n])
names(dataTrain)=nameall
names(dataTest)=nameall
}
data_gen[["dataTrain"]]<-dataTrain
data_gen[["dataTest"]]<-dataTest
data_gen[["w"]]<-w 
data_gen[["p"]]<-p
data_gen[["f"]]<-f
return (data_gen)
}




Data_gen<- function(gen_Data,function_num)
{
  rp_gen<-gen_Data
  rp_gen$transition<-rep(0,nrow(rp_gen))
  # 1-8 are from [1]
  if(function_num==1)
  {
    
  }else if(function_num==2){
    
    rp_gen$transition[rp_gen$V1>1]<-rp_gen$V1[rp_gen$V1>1]-5
   
  }else if(function_num==3){
    
    rp_gen$transition<-2 * rp_gen$V1 -4
    
    
  }else if(function_num==4){
    
    rp_gen$transition<-rp_gen$V2*rp_gen$V4*rp_gen$V6 + 2 *rp_gen$V2*rp_gen$V4*(1-rp_gen$V6)+
      4*rp_gen$V2*(1-rp_gen$V4)*(1-rp_gen$V6)+5*(1-rp_gen$V2)*rp_gen$V4*rp_gen$V6+6*(1-rp_gen$V2)*rp_gen$V4*(1-rp_gen$V6)+
      7*(1-rp_gen$V2)*(1-rp_gen$V4)*rp_gen$V6+7*(1-rp_gen$V2)*(1-rp_gen$V4)*rp_gen$V6+8*(1-rp_gen$V2)*(1-rp_gen$V4)*(1-rp_gen$V6)
    
    
  }else if(function_num==5){
    
    rp_gen$transition<-rp_gen$V1+rp_gen$V3+rp_gen$V5+rp_gen$V7+rp_gen$V8+rp_gen$V9-2
    
    
  }else if(function_num==6){
    
    rp_gen$transition<-2*rp_gen$V8*rp_gen$V9
    rp_gen$transition[rp_gen$V1>1&rp_gen$V3>0]<- rp_gen$transition[rp_gen$V1>1&rp_gen$V3>0]+4*rp_gen$V1[rp_gen$V1>1&rp_gen$V3>0]*rp_gen$V3[rp_gen$V1>1&rp_gen$V3>0]
    rp_gen$transition[rp_gen$V5>1&rp_gen$V7>0]<- rp_gen$transition[rp_gen$V5>1&rp_gen$V7>0]+4*rp_gen$V5[rp_gen$V5>1&rp_gen$V7>0]*rp_gen$V7[rp_gen$V5>1&rp_gen$V7>0]
    
    
  }else if(function_num==7){
    
    rp_gen$transition<-0.5*(rp_gen$V1^2+rp_gen$V2+rp_gen$V3^2+rp_gen$V4+rp_gen$V5^2+rp_gen$V6+rp_gen$V7^2+rp_gen$V8+rp_gen$V9^2-11)
    
  }else if(function_num==8){
    
    #rp_gen$transition<-(1/2^0.5)*(Data_gen(4,gen_Data)+Data_gen(5,gen_Data))
    rp_gen$transition<-(1/2^0.5)*(rp_gen$V2*rp_gen$V4*rp_gen$V6 + 2 *rp_gen$V2*rp_gen$V4*(1-rp_gen$V6)+
                                    4*rp_gen$V2*(1-rp_gen$V4)*(1-rp_gen$V6)+5*(1-rp_gen$V2)*rp_gen$V4*rp_gen$V6+6*(1-rp_gen$V2)*rp_gen$V4*(1-rp_gen$V6)+
                                    7*(1-rp_gen$V2)*(1-rp_gen$V4)*rp_gen$V6+7*(1-rp_gen$V2)*(1-rp_gen$V4)*rp_gen$V6+8*(1-rp_gen$V2)*(1-rp_gen$V4)*(1-rp_gen$V6)
                                  +rp_gen$V1+rp_gen$V3+rp_gen$V5+rp_gen$V7+rp_gen$V8+rp_gen$V9-2)
    
  }else {
    print("only selection 1-x")
  }
  return(rp_gen$transition)
}

y_gen<-function(mu,asym,tau,X,option_num,w,sig)
{
  w_<-1-w
  # option 1 is the y equals to the treatment effect
  if(option_num==1)
  {y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + rnorm(n,0,sig)}
  else if(option_num==-1)
  {
  y <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + rnorm(n,0,sig)
  }
  # option 2 is the y equals to the treatment effect adds a linear relation with other variables
  else if(option_num==2)
  {
    c<- seq(1,ncol(X),4)
    set.seed(123)
    w_x<- runif(length(c), 1, 10)
    y <- mu + asym*w*tau + (asym-1)*(1-w)*tau + w_x*X[,c] + rnorm(n,0,sig)
    
  }
  else if(option_num==-2)
  {
    c<- seq(1,ncol(X),4)
    set.seed(123)
    w_x<- runif(length(c), 1, 10)
    y <- mu + asym*w_*tau + (asym-1)*(1-w_)*tau + w_x*X[,c] + rnorm(n,0,sig)
    
  }
  # option 3 is the y equals to the treatment effect adds more complicated relation of other variables
  else if(option_num==3)
  {
    
  }
  else if(option_num==-3)
  {
    
  }
  # option 4 is the y presented the interaction among the treatmnet effect and other variables 
  else if(option_num==4)
  {
    
  }
  else if(option_num==-4)
  {
    
  }
  
  return(y)
}