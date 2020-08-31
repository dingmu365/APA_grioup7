
# if(!require("rlang")) install.packages("rlang")
 library("rlang")
# 
# if(!require("dplyr")) install.packages("dplyr")
 library(dplyr)
# 
# if(!require("data.table")) install.packages("data.table")
 library(data.table)
# 
#  install.packages('lazyeval')
#  if(!require("devtools")) install.packages("devtools");
# if(!require("rpart")) install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
# if(!require("rpart.plot")) install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
# if(!require("reshape2")) install.packages("reshape2", dependencies=TRUE, repos='http://cran.us.r-project.org')
# if(!require("plyr")) install.packages("plyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
# install.packages("processx")
#install.packages("usethis")
#install.packages("glue")
library(glue)
library(usethis)
library(processx)
library(devtools)



library(rpart)
library(rpart.plot)
library(reshape2)
#library(plyr)
install_github("susanathey/causalTree",force=TRUE)
library(causalTree)

source("causalMARS.R")
source("truncpow.R")
source("myridge.R")
source("predict.causalMARS.R")
source("makebx.newmars.R")
source("predict.bagged.causalMARS.R")
source("bagged.causalMARS.R")

library(Rmisc)

library("grf")
source("Scenarios.R")


c_tree<-function(data,dataTest,split.Bucket.temp=F,split.Rule.temp="policy"
                 ,cv.option.temp = "CT")
{
  # split.Bucket.temp=F
  # split.Rule.temp="policy"
  # cv.option.temp = "CT"
  
  dataTrain<-data[["data"]]
  dataTest<-dataTest[["data"]]
  ntr <-nrow(dataTrain)- round(.333*nrow(dataTrain))
  dataEst<-dataTrain[(ntr+1):nrow(dataTrain),]
  dataTrain<-dataTrain[1:ntr,]
  rownames(dataEst)<-1:nrow(dataEst)
  p <- data[["p"]]# number of total covariates
  f <- data[["f"]]# formulas for estimation
  ncov_sample<-floor(p/3) #number of covariates (randomly sampled) to use to build tree
  ncolx<-p #total number of covariates
  xvalvec = sample(5, nrow(dataTrain), replace=TRUE)
  
  tree<-honest.causalTree(as.formula(paste("y~",f)), data = dataTrain, treatment = dataTrain$w,est_data = dataEst, est_treatment = dataEst$w, split.Rule = split.Rule.temp, cv.option = cv.option.temp, split.Honest = T, cv.Honest = T, split.Bucket = split.Bucket.temp, xval = xvalvec, cp = 0, minsize = 25, propensity = 0.5)
  
  opcpid <- which.min(tree$cp[,4])
  opcp <- tree$cp[opcpid,1]
  tree_prune <- prune(tree, cp = opcp)
  
  # rpart.plot(tree_prune)
  #can use tree_prune or tree as trained causal tree model
  predicthonest <-predict(tree_prune,newdata=dataTest,type="vector")
  rm(tree_prune)
  tree<-NULL
  return(predicthonest)
}

c_forest<-function(data,dataTest,num_tree=100)
{
  #str_idx,end_idx
  dataTrain<-data[["data"]]
  dataTest<-dataTest[["data"]]
  w <- data[["w"]]
  p <- data[["p"]]# number of total covariates
  f <- data[["f"]]# formulas for estimation
  ncov_sample<-floor(p/3) #number of covariates (randomly sampled) to use to build tree
  ncolx<-p #total number of covariates
  #xvalvec = sample(5, nrow(dataTrain), replace=TRUE)
  cf <- causalForest(as.formula(paste("y~",f)), data=dataTrain, treatment=w,
                     split.Rule="CT", split.Honest=T,  split.Bucket=F, bucketNum = 5,bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 2L, split.alpha = 0.5, cv.alpha = 0.5,sample.size.total = floor(nrow(dataTrain) / 2), sample.size.train.frac = .5, mtry = ceiling(ncol(dataTrain)/3),nodesize = 3, num.trees=num_tree,ncolx=ncolx,ncov_sample=ncov_sample) 
  
  cfpredtest <- predict(cf, newdata=dataTest,predict.all = TRUE, type="vector")
  cf_result<-list(cf=cf, dataTrain=dataTrain, dataTest=dataTest, prediction=cfpredtest)
  rm(cf)
  
  return(cf_result)
}


c_cm<-function(data,dataTest, nbag =100)
{
  dataTrain <- data[["data"]]
  dataTest<- dataTest[["data"]]
  w <- data[["w"]]
  p <- data[["p"]]# number of total covariates
  f <- data[["f"]]# formulas for estimation
  fit_bcm = bagged.causalMARS(subset(dataTrain, select = -c(tau_true,w,y)), dataTrain$w, dataTrain$y, nbag)
  
  pred_bcm = predict(fit_bcm, newx =subset(dataTest, select = -c(tau_true,w,y))) 
  
  cm_result<-list(fit_bcm=fit_bcm, prediction=pred_bcm)
  return(cm_result)
}

coverage<- function(true_tau, CI)
{
  counting = sum(true_tau>=CI["lower"] & true_tau<= CI["upper"])
  return((counting+0.0)/length(true_tau))
}

CI_general<-function(data,accuracy=0.95,method)
{
  
  CI<-as.data.frame(NULL)
  for(i in 1:nrow(data))
  {
    n<-CI(as.numeric(data[i,]),accuracy)
    CI[i,"upper"]<-as.numeric(n["upper"]) 
    CI[i,"lower"]<-as.numeric(n["lower"])
  }
  
  return(CI)
}


library(knitr)
library(DescTools)


#wrapping all the models together and store all the results into a data.frame
wrapping_multiply_model<-function (mu_num ,
                                   tau_num ,
                                   option_num,
                                   p = 11, 
                                   n = 1000,
                                   times=500,
                                   num_tree=10,
                                   asym = .5,
                                   propens = .5,
                                   sig = .01,
                                   treatsize = .5,
                                   levsize = 1)
{
  dataTest<-Data_Gen(p,asym,n=10000,propens, sig, treatsize, levsize, mu_num, tau_num, option_num,seed_num=365)
  hat_tau_ct<-list()
  hat_tau_cf<-list()
  hat_tau_grf<-list()
  hat_tau_grf_sd<-list()
  hat_tau_cm<-list()
  data_stee<-list()
  
  # 
  # nrow(data[["data"]])
  for(i in 1:times)
  {
    data=NULL
    ct=NULL
    data<-Data_Gen(p,asym,n,propens, sig, treatsize, levsize, mu_num, tau_num,
                   option_num,seed_num=i*100)
    
    ct<-c_tree(data,dataTest)
    hat_tau_ct[[i]]<-ct
    
    cf<-c_forest(data,dataTest,num_tree)
    hat_tau_cf[[i]]<-cf[["prediction"]]
    
    general_rf<-c_grf(data,dataTest,num.tree=num_tree)
    hat_tau_grf[[i]]<-general_rf[["prediction"]]
    
    causal_cm<-c_cm (data,dataTest)
    hat_tau_cm[[i]]<-causal_cm[["prediction"]]
    print("i")
  } 
  
  true_tau<-dataTest[["data"]] $tau_true
  result_tb<-data.frame(NULL)
  ##################################################################  
  data_pannel_ct<-data.frame(NULL)
  
  hat_tau_ct<-hat_tau_ct
  for (i in 1:n)
  {
    
    trans<-vector()
    for(j in 1:length(hat_tau_ct))
    {
      trans[j]<-hat_tau_ct[[j]][i]
    }
    data_pannel_ct<-rbind(data_pannel_ct,trans)
  }
  nrow(data_pannel_ct)
  colnames(data_pannel_ct)<-paste0("time_",1:length(hat_tau_ct))
  
  index=1
  #calculate the MSE of causal tree
  result_tb[index,1]=MSE(rowMeans(data_pannel_ct),true_tau)
  #calculate the CI of ct
  CI_ct<-CI_general(data_pannel_ct,0.95)    
  #calculate the CI mean_width of 
  result_tb[index,2]=as.numeric(colMeans(CI_ct["upper"]-CI_ct["lower"]))
  #calculate the CI covertage percentage of causal tree
  result_tb[index,3]=coverage(true_tau ,CI_ct)
  #calculate the ration of CI covertage/mean_width of causal tree
  result_tb[index,4]=( result_tb[index,3]+0.0)/ result_tb[index,2]
  #calculate the CI_q of causal tree
  CI_ct_q<- CI_quantile(data_pannel_ct,0.95,method = "no_tree")
  #calculate the CI_q mean_width of causal tree
  result_tb[index,5]=as.numeric(colMeans(CI_ct_q["upper"]-CI_ct_q["lower"]))
  #calculate the CI_q covertage percentage of causal tree
  result_tb[index,6]=coverage(true_tau ,CI_ct_q)
  #calculate the ration of CI_q covertage/mean_width of causal tree
  result_tb[index,7]=( result_tb[index,6]+0.0)/ result_tb[index,5]
  
  ##################################################################
  data_pannel_cf<-data.frame(NULL)
  
  hat_tau_cf<-hat_tau_cf
  for (i in 1:n)
  {
    
    trans<-vector()
    for(j in 1:length(hat_tau_cf))
    {
      trans[j]<-hat_tau_cf[[j]][["aggregate"]][i]
    }
    data_pannel_cf<-rbind(data_pannel_cf,trans)
  }
  nrow(data_pannel_cf)
  colnames(data_pannel_cf)<-paste0("time_",1:length(hat_tau_cf))
  
  index=2
  #calculate the MSE of causal forest
  result_tb[index,1]<-MSE(rowMeans(data_pannel_cf),true_tau)
  #calculate the CI of causal forest
  CI_cf<-CI_general(data_pannel_cf, 0.95)
  #calculate the CI mean_width of causal forest
  result_tb[index,2]=as.numeric(colMeans(CI_cf["upper"]-CI_cf["lower"]))
  #calculate the CI covertage percentage of causal forest
  result_tb[index,3]=coverage(true_tau ,CI_cf)
  #calculate the ration of CI covertage/mean_width of  causal forest
  result_tb[index,4]=( result_tb[index,3]+0.0)/ result_tb[index,2]
  #calculate the CI_q of causal forest
  CI_cf_q<-CI_quantile(data_pannel_cf,0.95,method = "no_tree")
  #calculate the CI_q mean_width of causal forest
  result_tb[index,5]=as.numeric(colMeans(CI_cf_q["upper"]-CI_cf_q["lower"]))
  #calculate the CI_q covertage percentage of causal forest
  result_tb[index,6]=coverage(true_tau ,CI_cf_q)
  #calculate the ration of CI_q covertage/mean_width of causal forest 
  result_tb[index,7]= (result_tb[index,6]+0.0)/ result_tb[index,5]
  
  ##################################################################
  data_pannel_grf<-data.frame(NULL)
  
  hat_tau_grf<-hat_tau_grf
  for (i in 1:n)
  {
    
    trans<-vector()
    for(j in 1:length(hat_tau_grf))
    {
      trans[j]<-hat_tau_grf[[j]][i,"predictions"]
    }
    data_pannel_grf<-rbind(data_pannel_grf,trans)
  }
  nrow(data_pannel_grf)
  colnames(data_pannel_grf)<-paste0("time_",1:length(hat_tau_grf))
  
  index=3
  #calculate the MSE of grf  
  result_tb[index,1]=MSE(rowMeans(data_pannel_grf),true_tau)
  #calculate the CI of grf
  CI_grf<-CI_general(data_pannel_grf,0.95)
  #calculate the CI mean_width of grf
  result_tb[index,2]<-as.numeric(colMeans(CI_grf["upper"]-CI_grf["lower"]))
  #calculate the CI covertage percentage of grf
  result_tb[index,3]= coverage(true_tau ,CI_grf)
  #calculate the ration of CI covertage/mean_width of grf
  result_tb[index,4]=(result_tb[index,3]+0.0)/result_tb[index,2]
  #calculate the CI_q of grf
  CI_grf_q<-CI_quantile(data_pannel_grf,0.95,method = "no_tree")
  #calculate the CI_q mean_width of grf
  result_tb[index,5]=as.numeric(colMeans(CI_grf_q["upper"]-CI_grf_q["lower"]))
  #calculate the CI_q covertage percentage of grf
  result_tb[index,6]=coverage(true_tau ,CI_grf_q)
  #calculate the ration of CI_q covertage/mean_width of grf
  result_tb[index,7]=(result_tb[index,6]+0.0)/result_tb[index,5]
  
  ##################################################################
  data_pannel_cm<-data.frame(NULL)
  hat_tau_cm<-hat_tau_cm
  for (i in 1:n)
  {
    
    trans<-vector()
    for(j in 1:length(hat_tau_cm))
    {
      trans[j]<-hat_tau_cm[[j]][i]
    }
    data_pannel_cm<-rbind(data_pannel_cm,trans)
  }
  nrow(data_pannel_cm)
  colnames(data_pannel_cm)<-paste0("time_",1:length(hat_tau_cm))
  index=4
  #calculate the MSE of cm
  result_tb[index,1]=MSE(rowMeans(data_pannel_cm),true_tau)
  #calculate the CI of cm
  CI_cm<-CI_general(data_pannel_cm,0.95)
  #calculate the CI mean_width of cm
  result_tb[index,2]<-as.numeric(colMeans(CI_cm["upper"]-CI_cm["lower"]))
  #calculate the CI covertage percentage of cm
  result_tb[index,3]<-coverage(true_tau ,CI_cm)
  #calculate the ration of CI covertage/mean_width of cm
  result_tb[index,4]=( result_tb[index,3]+0.0)/ result_tb[index,2]
  
  #calculate the CI_q of cm
  CI_cm_q<-CI_quantile(data_pannel_grf,0.95,method = "no_tree")
  #calculate the CI_q mean_width of cm 
  result_tb[index,5]<-as.numeric(colMeans(CI_cm_q["upper"]-CI_cm_q["lower"]))
  #calculate the CI_q covertage percentage of cm
  result_tb[index,6]<-coverage(true_tau ,CI_cm_q)
  #calculate the ration of CI_q covertage/mean_width of cm
  result_tb[index,7]<-(  result_tb[index,6]+0.0)/  result_tb[index,5]
  
  colnames(result_tb)<-c("MSE","mean_width","coverage","ratio","mean_width_q","coverage_q","ratio_q");
  rownames(result_tb)<-c("ct","cf","grf","cm")
  
  return(result_tb)
}


mu_num = 2
tau_num = 7
option_num = 1
p=50
n=500
tb_7<-wrapping_multiply_model (mu_num=2,tau_num = 7, option_num= 1,times= 500,p=p,n=n,num_tree=50)

# 
# tb_8[1,]=c(5.999975,0.14290587,0.0004,0.002799045,3.145833,0.0047,0.001494040)
# tb_8[2,]=c(6.021385,0.08767860,0.0001,0.001140529,1.930145,0.0032,0.001657907)
# tb_8[3,]=c(6.010697,0.08066666,0.0002,0.002479339,1.782521,0.0027,0.001514708)
# tb_8[4,]=c(9.739055,0.30477674,0.0016,0.005249744,1.782521,0.0027,0.001514708)



