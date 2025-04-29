#Code for Medical Data Examples: Diagnosis Alzheimer With Handwriting Dataset

#Here we apply three forward methods: logistic, nonparametric minmax and WYIM to DARWIN_30 data, the data has been preprocessed
#For each method we record: linear combination coefficient vector, cutoff point and according WYI in final model
#The plots show conditional density curves of linear combination x%*%bhat given y=0 and y=1
#Here we set w=0.5. The w can be set between 0 and 1 

source("WYIM_code.R")
library(SLModels) #for nonparametric minmax method
library(ggplot2) #for density curve
library(patchwork)

data=read.csv("DARWIN_30.csv")

x=as.matrix(data[,2:31])
y=data[,32]
w=0.5
d=dim(x)[2]


###-------------------method 1:logistic forward-----------------------------------
bhat_logistic=matrix(0,d,d)
youden_logistic=rep(0,d)
mu_logistic=rep(0,d)

uni_youden=rep(0,d)
for(h in 1:d){   #evaluate the WYI with x[,h]
  uni_youden[h]=youden_w(y,x[,h],w)$WYI
}
max_ind=which.max(uni_youden)
x=x[,c(max_ind,setdiff(1:d, max_ind))]  #bring covariate with max WYI to column 1

#forward selection start with one covariate
mu_logistic[1]=youden_w(y,x[,1],w)$cutoff 
youden_logistic[1]=youden_w(y,x[,1],w)$WYI 
bhat_logistic[1,]=c(1,rep(0,d-1)) #beta, WYI and cutoff with one covariate
zero_ind=which(bhat_logistic[1,]==0)
nonzero_ind=which(bhat_logistic[1,]!=0) 
best_youden=1

for(i in 2:d){ #with i-1 covariates selected, search for i_th covariate
  youden_path=rep(0,d-i+1)
  b_path=matrix(0,d-i+1,d)
  mu_path=rep(0,d-i+1)
  for(j in 1:(d-i+1)){ #evaluate the WYI when adding the remaining d-i+1 cadidate covariates as i_th covariate
    out<-glm(y~x[,c(nonzero_ind,zero_ind[j])]-1,family=binomial(link=logit))$coef    
    b_path[j,c(nonzero_ind,zero_ind[j])]<-out/sqrt(sum(out^2))
    youden_path[j]=youden_w(y,x[,c(nonzero_ind,zero_ind[j])]%*%b_path[j,c(nonzero_ind,zero_ind[j])],w)$WYI
    mu_path[j]=youden_w(y,x[,c(nonzero_ind,zero_ind[j])]%*%b_path[j,c(nonzero_ind,zero_ind[j])],w)$cutoff
    
  }
  youden_logistic[i]=max(youden_path) #select i_th covariate according to max WYI
  bhat_logistic[i,]=b_path[which.max(youden_path),]
  mu_logistic[i]=mu_path[which.max(youden_path)]
  zero_ind=which(bhat_logistic[i,]==0) #index for remaining covariates 
  nonzero_ind=which(bhat_logistic[i,]!=0) #index for i covariates selected
  
  
  if(youden_logistic[i]>youden_logistic[i-1]){
    best_youden=i
  }
  if(i-best_youden>=2){ #terminate when two successive variable additions failed to enhance WYI
    break
  }
}
#result for logistic forward:number of selected covariates, bhat, cutoff, WYI in final model
result_logistic=list(d=best_youden,bhat_logistic=bhat_logistic[best_youden,],
                     cutoff_logistic=mu_logistic[best_youden],
                     WYI_logistic=youden_logistic[best_youden])
score_logistic=x%*%bhat_logistic[best_youden,]
curve_logistic_Az=data.frame(score_logistic,y=as.factor(y))



####-------------------------method 2: nonparametric minmax forward selection---------------------------------------------
bhat_mm=matrix(0,d,2)
youden_mm=rep(0,d)
mu_mm=rep(0,d)

uni_youden=rep(0,d)
for(h in 1:d){   #evaluate the WYI with x[,h]
  uni_youden[h]=youden_w(y,x[,h],w)$WYI
}
max_ind=which.max(uni_youden)
x=x[,c(max_ind,setdiff(1:d, max_ind))]  #bring covariate with max WYI to column 1

#forward selection start with one covariate
mu_mm[1]=youden_w(y,x[,1],w)$cutoff
youden_mm[1]=youden_w(y,x[,1],w)$WYI
bhat_mm[1,]=c(1,0) #beta, WYI and cutoff with one covariate
index=rep(0,d) #index of covariate selected in each round
zero_ind=seq(2,30)
nonzero_ind=1
best_youden=1

for(i in 2:d){ #with i-1 covariates selected, search for i_th covariate
  ks_path=rep(0,d-i+1)
  b_path=matrix(0,d-i+1,2) #coefficient of mm is 2-dimensional
  t_path=rep(0,d-i+1)
  for(j in 1:(d-i+1)){ #evaluate the WYI when adding the remaining d-i+1 cadidate covariates as i_th covariate
    out=SLModels(cbind(data.frame(x[,c(nonzero_ind,zero_ind[j])]),y),"minmax")
    bhat_minmax<-as.numeric(out[1:2])
    b_path[j,]<-bhat_minmax/sqrt(sum(bhat_minmax^2))
    ks_path[j]=youden_w(y,cbind(apply(x[,c(nonzero_ind,zero_ind[j])],1,max),apply(x[,c(nonzero_ind,zero_ind[j])],1,min))%*%b_path[j,],w)$WYI
    t_path[j]=youden_w(y,cbind(apply(x[,c(nonzero_ind,zero_ind[j])],1,max),apply(x[,c(nonzero_ind,zero_ind[j])],1,min))%*%b_path[j,],w)$cutoff
  }
  youden_mm[i]=max(ks_path) #select i_th covariate according to max WYI
  index[i]=zero_ind[which.max(ks_path)] #index of i_th covariate selected
  bhat_mm[i,]=b_path[which.max(ks_path),]
  nonzero_ind=c(nonzero_ind,index[i]) #index for i covariates selected
  zero_ind=zero_ind[-which.max(ks_path)] #index for remaining covariates 
  mu_mm[i]=t_path[which.max(ks_path)]
  
  if(youden_mm[i]>youden_mm[i-1]){
    best_youden=i
    x_ind=nonzero_ind
  }
  if(i-best_youden>=2){ #terminate when two successive variable additions failed to enhance WYI
    break
  }
}
#result for mm forward: bhat, cutoff, WYI and index of slected covariates in final model
result_mm=list(bhat_mm=bhat_mm[best_youden,],
                     cutoff_mm=mu_mm[best_youden],
                     WYI_mm=youden_mm[best_youden],ind_selected=x_ind)
score_mm=cbind(apply(x[,x_ind],1,max),apply(x[,x_ind],1,min))%*%bhat_mm[best_youden,]
curve_mm_Az=data.frame(score_mm,y=as.factor(y))


###--------------------method 3:WYIM forward--------------------------------------
out=WYIM_forward(x,y,w)
bhat_wyim=out$bhat

#result for mm forward: number of covariates, bhat, cutoff, WYI in final model
result_wyim=list(d=out$d,
                 bhat_wyim=out$bhat,
                 cutoff_wyim=out$cutoff,
                 WYI_wyim=out$WYI)
score_wyim=x%*%bhat_wyim
curve_wyim_Az=data.frame(score_wyim,y=as.factor(y))

#bhat, cutoff and WYI during WYIM forward selection
forward_wyim=list(bhat_forward=out$bhat_forward, 
                  cutoff_forward=out$cutoff_forward,
                  WYI_forward=out$WYI_forward)
###--------------------------present result----------------------------------------
result_wyim
result_logistic
result_mm

##conditional density curves for different methods----------------------------------
density_logistic=
  ggplot(curve_logistic_Az, aes(x =score_logistic, color = y)) +
  geom_density(size=0.7,bw=0.7) +  
  geom_vline(xintercept =result_logistic$cutoff_logistic,linetype = "dashed",size = 0.7,color='navy') +
  ylab("")+xlab("")+
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC2"))+ # 0 red
  scale_x_continuous(limits = c(-5,5))+scale_y_continuous(limits = c(0,0.56))+
  theme_minimal()+
  theme(axis.text = element_text(color="black",size = 15),title = element_text(color="black",size = 15),legend.position = "none")+
  annotate("text", x = 2.5, y = 0.5, label = "logistic",size=6)

density_mm=
  ggplot(curve_mm_Az, aes(x =score_mm, color = y)) +
  geom_density(size=0.7,bw=0.7) +  
  geom_vline(xintercept =result_mm$cutoff_mm,linetype = "dashed",size = 0.7,color='navy') +
  ylab("")+xlab("")+
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC2"))+ # 0 red
  scale_x_continuous(limits = c(-5,5))+scale_y_continuous(limits = c(0,0.56))+
  theme_minimal()+
  theme(axis.text = element_text(color="black",size = 15),title = element_text(color="black",size = 15),legend.position = "none")+
  annotate("text", x = 2.5, y = 0.5, label = "MM",size=6)

density_wyim=
  ggplot(curve_wyim_Az, aes(x =score_wyim, color = y)) +
  geom_density(size=0.7,bw=0.7) +  
  geom_vline(xintercept =result_wyim$cutoff_wyim,linetype = "dashed",size = 0.7,color='navy') +
  ylab("")+xlab("")+
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC2"))+ # 0 red
  scale_x_continuous(limits = c(-5,5))+scale_y_continuous(limits = c(0,0.56))+
  theme_minimal()+
  theme(axis.text = element_text(color="black",size = 15),title = element_text(color="black",size = 15),legend.position = "none")+
  annotate("text", x = 2.5, y = 0.5, label = "WYIM",size=6)

density_Az=density_wyim+density_logistic+density_mm+plot_layout(ncol = 3)
ggsave("density_Az.pdf",width =15,height =4.5,family="Helvetica")
