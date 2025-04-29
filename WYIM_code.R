#Here we provide three functions used in our paper:youden_w, WYIM and WYIM_forward



##---------------------------------------------------------------------------------------------------
##calculate WYI based on binary Y, the score and weight 0<w<1----------------------------------------
##return WYI and cutoff point------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------
youden_w=function(y,score,w=0.5){
  n0=sum(y==0)
  n1=sum(y==1)
  n=n0+n1
  
  n0_w=n0/(2-2*w) 
  n1_w=n1/(2*w)       
  muL<-sort(score)         #generate ascending-sorted scores as cutoff value candidates  
  y<-y[order(score)]      #sort y in ascending order of score
  youden<-numeric(n) 
  
  if(y[1]==0){             #calculate WYI with cutoff muL[1]
    youden[1]=1/n0_w
  }else{
    youden[1]=-1/n1_w
  }
  
  for(i in 2:n){               #calculate WYI with cutoff from muL[2] to muL[n]
    if(y[i]==0){
      youden[i]<-youden[i-1]+1/n0_w
    }else{
      youden[i]<-youden[i-1]-1/n1_w
    }
  }
  
  d<-length(which(youden==max(youden)))  #count occurrences of max WYI value
  muhat<-muL[which(youden==max(youden))[d]]  # select muL vector as cutoff according to last occurrence of max WYI
  
  return(list(WYI=max(youden)+2*w-1,cutoff=muhat))  
}

##---------------------------------------------------------------------------------------------------
##Code of proposed WYIM method-----------------------------------------------------------------------
##calculate coefficients of linear combination, cutoff based on binary Y, non-zero covariates x, 0<w<1 and initial value beta_start---------------
##return coefficient vector, iteration round, according WYI and cutoff-------------------------------
##---------------------------------------------------------------------------------------------------
WYIM=function(x,y,w=0.5,beta_start,maxiter=500,tol=10^(-6)){
  n0=sum(y==0)
  n1=sum(y==1)
  n=n0+n1
  n0_w=n0/(2-2*w)
  n1_w=n1/(2*w)
  
  p<-dim(x)[2]  # number of covariates
  
  if(p==1){ 
    beta2<-1
    i<-0
    youden_max<-youden_w(y,x,w)$WYI	#bhat=1 when p=1
  }
  
  if(p>1){
    beta1<-beta_start
    out<-youden_w(y,x%*%beta1,w)
    mu1<-out$cutoff   # use (beta1,mu1) as the starting point
    youden_max<-out$WYI
    
    beta2<-beta1    # each iteration update (beta1,mu1) to (beta2,mu2)
    mu2<-mu1
    
    i<-1   # iteration indicator
    while(i<maxiter){         
      for(k in 1:p){	  # at the ith iteration, update beta[k] one by one
        if(abs(beta2[k])>0){	 # if one beta[k] shrink to 0 at one iteration, then leave it as 0 afterwards	
          cutoff<-sort((mu1-as.matrix(x[,-k])%*%beta2[-k])/x[,k])  #generate ascending-sorted break points as beta2[k] value candidates  
          ord<-order((mu1-as.matrix(x[,-k])%*%beta2[-k])/x[,k])
          xk<-x[ord,k]
          yy<-y[ord]	#sort xk and y in ascending order of break point
          
          youden_0<-rep(0,n)
          youden_1<-rep(0,n)
          
          
          youden_0[1]<-sum(cutoff[1]*xk[which(yy==0)]<=cutoff[which(yy==0)]*xk[which(yy==0)])/n0_w
          youden_1[1]<-sum(cutoff[1]*xk[which(yy==1)]<=cutoff[which(yy==1)]*xk[which(yy==1)])/n1_w
          
          for(ii in 2:n){ #evaluate youden_0 and youden_1 at each break point by following rules
            if(yy[ii-1]==0 & yy[ii]==0){
              if(xk[ii-1]>0 & xk[ii]>0){
                youden_0[ii]<-youden_0[ii-1]-1/n0_w
              }
              if(xk[ii-1]<0 & xk[ii]<0){
                youden_0[ii]<-youden_0[ii-1]+1/n0_w
              }
              if(xk[ii-1]*xk[ii]<0){
                youden_0[ii]<-youden_0[ii-1]
              }							
              youden_1[ii]<-youden_1[ii-1]						
            }
            if(yy[ii-1]==0 & yy[ii]==1){
              if(xk[ii-1]>0){
                youden_0[ii]<-youden_0[ii-1]-1/n0_w
              }else{youden_0[ii]<-youden_0[ii-1]}
              if(xk[ii]<0){
                youden_1[ii]<-youden_1[ii-1]+1/n1_w
              }else{youden_1[ii]<-youden_1[ii-1]}
            }
            if(yy[ii-1]==1 & yy[ii]==0){
              if(xk[ii]<0){
                youden_0[ii]<-youden_0[ii-1]+1/n0_w
              }else{youden_0[ii]<-youden_0[ii-1]}
              if(xk[ii-1]>0){
                youden_1[ii]<-youden_1[ii-1]-1/n1_w
              }else{youden_1[ii]<-youden_1[ii-1]}
            }
            if(yy[ii-1]==1 & yy[ii]==1){
              if(xk[ii-1]>0 & xk[ii]>0){
                youden_1[ii]<-youden_1[ii-1]-1/n1_w
              }
              if(xk[ii-1]<0 & xk[ii]<0){
                youden_1[ii]<-youden_1[ii-1]+1/n1_w
              }
              if(xk[ii-1]*xk[ii]<0){
                youden_1[ii]<-youden_1[ii-1]
              }							
              youden_0[ii]<-youden_0[ii-1]	
            }
          }						
          
          youden_<-youden_0-youden_1   #obtain WYI at each break point
          d<-length(which(youden_==max(youden_)))  #count occurrences of max WYI value
          beta2[k]<-cutoff[which(youden_==max(youden_))[d]]	# select break points as beta2[k] according to last occurrence of max WYI
        }			
      }
      
      if(sum(beta2^2)==0){
        youden_max<-youden_w(y,rep(0,n),w)$WYI
        break
      }
      
      beta2<-beta2/sqrt(sum(beta2^2)) #normalize beta2 to unit length (L2-normalization)
      out<-youden_w(y,x%*%beta2,w) 
      mu2<-out$cutoff #cutoff with score x%*%beta2
      youden_max<-out$WYI # WYI with score x%*%beta2
      
      if(abs(1-sum(beta1*beta2))<tol){ # check convergence of beta
        break # terminate iteration when beta align
      }
      
      beta1<-beta2 #update beta vector
      mu1<-mu2 #update cutoff
      i<-i+1		
    }	
  }
  return(list(bhat=beta2,iter=i,WYI=youden_max,cutoff=mu2)) 
}


##---------------------------------------------------------------------------------------------------
##code for forward variable selection of WYIM method-------------------------------------------------
##calculate the coefficient vector and cutoff point after forward selection--------------------------
##return number of covariates selected, coefficient vector, cutoff point and WYI after forward selection and those values during forward selection
##---------------------------------------------------------------------------------------------------
WYIM_forward<-function(x,y,w=0.5,maxiter=500,tol=10^(-6)){
  
  d=dim(x)[2] #dimension of covariates
  beta=matrix(0,d,d)
  youden=rep(0,d)
  mu=rep(0,d)
  
  uni_youden=rep(0,d)
  for(h in 1:d){   #evaluate the WYI with x[,h]
    uni_youden[h]=youden_w(y,x[,h],w)$WYI
  }
  nonzero_ind=which.max(uni_youden)   
  beta[1,nonzero_ind]=1              #select first covariate into the model according to max WYI
  zero_ind=which(beta[1,]==0)              #index of remaining covariates
  youden[1]=uni_youden[nonzero_ind]              # WYI with one covariate
  mu[1]=youden_w(y,x%*%beta[1,],w)$cutoff              # cutoff with covariate
  best_youden=1
  
  #forward selection
  for(i in 2:d){      #with i-1 covariates selected, search for i_th covariate 
    youden_path=rep(0,d-i+1) 
    mu_path=rep(0,d-i+1) 
    b_path=matrix(0,d-i+1,d) 
    
    for(j in 1:(d-i+1)){   #evaluate the WYI when adding the remaining d-i+1 cadidate covariates as i_th covariate
      start=beta[i-1,] 
      start[zero_ind[j]]=1/sqrt(i)   #enable x[,zero_ind[j]] into update process since in WYIM method if beta[zero_ind[j]] shrink to 0 then leave it as 0 afterwards 
      beta_start=start/sqrt(sum(start))   #normalize beta_start to unit length
      
      out=WYIM(x,y,w,beta_start=beta_start,maxiter=maxiter,tol=tol)
      youden_path[j]=out$WYI   #WYI while adding beta[zero_ind[j]]
      b_path[j,]=out$bhat    #coefficient vector while adding beta[zero_ind[j]]
      mu_path[j]=out$cutoff
    }
    youden[i]=max(youden_path) #select i_th covariate according to max WYI
    mu[i]=mu_path[which.max(youden_path)]
    beta[i,]=b_path[which.max(youden_path),] #coefficient vector and cutoff with i covariates selected
    zero_ind=which(beta[i,]==0)  #index for remaining covariates 
    nonzero_ind=which(beta[i,]!=0)  #index for i covariates selected
    
    
    if(youden[i]>youden[i-1]){
      best_youden=i
    }
    if(i-best_youden>=2){ #terminate when two successive variable additions failed to enhance WYI
      break
    }
  }
  
  return(list(d=best_youden,bhat=beta[best_youden,],WYI=youden[best_youden],cutoff=mu[best_youden],
              bhat_forward=beta[1:best_youden,],WYI_forward=youden[1:best_youden],cutoff_forward=mu[1:best_youden]))
 #p,bhat,WYI,cutoff note for final results and bhat_forward,WYI_forward,cutoff_forward show those values during forward selection
  }