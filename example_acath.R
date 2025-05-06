#Code for Medical Data Examples: Duke Cardiac Catheterization Diagnostic Dataset

#Here we apply four methods: logistic, single index, nonparametric minmax and WYIM to acath_male60 data, the data has been preprocessed
#For each method we record: computing time, linear combination coefficient vector, cutoff point and according WYI
#Here we set w=n1/n, another situation in paper is w=0.5. The w can be set between 0 and 1 and remenber to change the title and names of plots with w

source("WYIM_code.R")
library(np) #for single index method
library(SLModels) #for nonparametric minmax method

data=read.csv("acath_male60.csv")

w=472/523 #n1/n
#w=0.5 #set another weight

y<-as.matrix(data[,5])
x<-as.matrix(data[,2:4])


######------------Method 1: logistic regression----------------------------------- 
start_time <- Sys.time()
out<-glm(y~x-1,family=binomial(link=logit))$coef
end_time <- Sys.time()

bhat_logistic=out/sqrt(sum(out^2)) 

result_logistic=list(time_logistic=as.numeric(end_time - start_time),
                     bhat_logistic=bhat_logistic,
                     cutoff_logistic=youden_w(y,x%*%bhat_logistic,w)$cutoff,
                     WYI_logistic=youden_w(y,x%*%bhat_logistic,w)$WYI)
                     

######----------- Method 2: single index------------------------------------------
start_time <- Sys.time()
out=npindexbw(x,as.vector(y))
end_time <- Sys.time()

bhat_si=out$beta
bhat_si=bhat_si/sqrt(sum(bhat_si^2))

result_si=list(time_si=as.numeric(end_time - start_time),
               bhat_si=bhat_si,
               cutoff_si=youden_w(y,x%*%bhat_si,w)$cutoff,
               WYI_si=youden_w(y,x%*%bhat_si,w)$WYI)



######------------Method 3: non-parametric minimax----------------------------------
start_time <- Sys.time()
out=SLModels(cbind(data.frame(x),y),"minmax")
end_time <- Sys.time()

bhat_mm=as.numeric(out[1:2])
bhat_mm=bhat_mm/sqrt(sum(bhat_mm^2))

result_mm=list(time_mm=as.numeric(end_time - start_time),
               bhat_mm=bhat_mm,
               cutoff_mm=youden_w(y,cbind(apply(x,1,max),apply(x,1,min))%*%bhat_mm,w)$cutoff,
               WYI_mm=youden_w(y,cbind(apply(x,1,max),apply(x,1,min))%*%bhat_mm,w)$WYI)
#score of mm is bhat_mm[1]*max(x)+bhat_mm[2]*min(x)


######------------Method 4: DM youden start with logistic---------------------------
start_time <- Sys.time()
out=WYIM(x,y,w,beta_start=rep(1/sqrt(3),3)) #beta start from equal coefficients for each covariate
end_time <- Sys.time()

bhat_wyim=out$bhat

result_wyim=list(time_dmks=as.numeric(end_time - start_time),
                 bhat_wyim=bhat_wyim, 
                 cutoff_wyim=out$cutoff,
                 WYI_wyim=out$WYI)

###---------------present the results--------------------------------------------
result_wyim
result_logistic
result_si
result_mm


