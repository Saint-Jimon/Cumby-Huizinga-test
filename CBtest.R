#######################################
########## Cumby-Huizinga #############
#######################################

#Cumby-Huizinga autocorrelation test
#Author: Stijn Jansen
#Based on: Cumby, R. and Huizinga J., (1992) Econometrica; 
#ivactest Baum; Schaffer

CBtest<-function(x,s){
  model <- plm:::describe(x, "model")
  effect <- plm:::describe(x, "effect")
  theta <- x$ercomp$theta

  X <- model.matrix(x,model=model,effect=effect,theta=theta,rhs=1)
  y <- pmodel.response(model.frame(x), model = model, effect = effect,theta=theta)
  W <- model.matrix(x,model=model,effect=effect, theta=theta ,rhs = 2)
  res<-residuals(x)
  
  
  T<-length(y)
  k<-length(X[1,])
  h<-length(W[1,])
  sigma_sq<-(crossprod(res,res)/T)[1]
  U<-matrix(0,nrow=T,ncol=s)
#   for(i in 1:T){
#     for(j in 1:s){
#       if((i-j-s)>0){ U[i,j]<-x$res[i-j-s]}
#       else U[i,j]<-0
#      
#     }d
#   }
  
  for(j in 1:s){
    for(i in 1:T){
      if((i-j)>0) U[i,j]<-res[(i-j)] #q=0
    }
  }
  
  B<- -(crossprod(U,X)/T)/sigma_sq
  A_t<- crossprod(W,W)/T
  D<- T*solve(crossprod(X,W)%*%solve(A_t)%*%crossprod(W,X))%*% crossprod(X,W)%*%solve(A_t)
  phi<-matrix(0,nrow=2*s,ncol=h+s)
  phi[1,1:h]<-B%*%D  
  phi[2,h+s]<- diag(s)/sigma_sq
  
  eta1<-res*W
  eta2<-res*U
  
  
  
  eta<-t(cbind(eta1,eta2))

#create psi--------------
  add <- function(x) Reduce("+", x)
  N<-ceiling(T^(0.25))
  
#(26) Gaussian weight = 1 for n = 0
  l0<-vector("list",length=(T))
  for(t in 1:(T)) {
    l0[[t]] = eta[,t]%*%t(eta[,t])
  }
  R0<-add(l0)/T
	#psi<- R0

	pvecp = vector("list",length=N)
	pvecm = vector("list",length=N)
	
	
	for(n in 1:N) {
	  
	#positive part N
	  lp<-vector("list",length=(T-n-1))
	  for(t in (n+1):(T)) {
	    lp[[t-n]] = eta[,t]%*%t(eta[,(t-n)])
	  }
	  pvecp[[n]]<-add(lp)/length(lp) #ipv /T
	 
	  #negative part N 
	  
	  lm<-vector("list",length=(T-n-1))
	    for(t in 1:(T-n)) {
	     lm[[t]]<-eta[,t]%*%t(eta[,(t+n)])
	    }
	  pvecm[[n]] =add(lm)/length(lm)  #ipv /T
	}
	#crossprod(eta[,1:(T-n)],eta[,])
	
	tnsq = 2*(N^2)
	

	w_n[1:N]<-exp((-(1:N)^2)/tnsq)
	
	for(i in 1:N){
	  pvecp[[i]]<- w_n[i]*pvecp[[i]]
	  pvecm[[i]]<- w_n[i]*pvecm[[i]]
	}
	
	
	pvesum<-vector("list",length=N)
	
	for(n in 1:N){
	  pvesum[[n]]<-pvecp[[n]]+pvecm[[n]]
	}
	
	psi<-add(pvesum)+R0
	
# 	psilist<-vector("list",length=N)
# 
# 	for(i in 1:N) {
# 	  psilist[[i]] = exp(-(i^2)/tnsq) * (pvesum[[i]])
# 	}
# 	psi<-add(psilist)+R0

	
	
	#Rn en volgende opnieuw
	
	
	
	
	


# Vr<-crossprod(eta2,eta2)*sqrt(sigma_sq)^(-4)
# C<-crossprod(eta2,eta1)*sqrt(sigma_sq)^(-2)
# omega<-crossprod(W,x$res)%*%crossprod(x$res,W)/T
# 
# psi<-matrix(nr=h+s,nc=h+s)
# psi[1:h,1:h]<-omega
# psi[(h+1):(h+s),1:h]<- sigma_sq*C
# psi[1:h,(h+1):(h+s)]<-sigma_sq*t(C)
# psi[(h+1):(h+s),(h+1):(h+s)]<-Vr*sigma_sq^2
#   
# #   
# #   eta_temp<-matrix(nrow=T,ncol=h+s) 
# #   for (t in (s+1):T){
# #     for(j in 1:h){
# #       eta_temp[t,j]<-x$res[t]*W[t,j]
# #     }
# #     for(j in (h+1):(h+s)){
# #       #print(t-q-(j-h+1))
# #       eta_temp[t,j]<-x$res[t]*x$res[(t-j+h+1)]
# #     }
# #   }
# #   
# #   #eta<-t(eta_temp[(s+1):T,]) #eerste s kolommen zijn verwijderd
# #   eta<-t(eta_temp)
# #   
# #   #code based on stata ivactest
#   N<-ceiling(T^(0.25))
#   #n=0
#   l0<-vector("list",length=T-s)
#    for(t in (1):(T-1)){
#     l0[[t]]<-tcrossprod(eta[,t],eta[,t])
#    }
#   R0<-add(l0)
#   psi<-R0
#   
#   
# 
#   
#   lp<-vector("list",length=N)
#   lm<-vector("list",length=N)
#   
#   
#   for(i in 1:N) {
#     for(t in (i+1):(T-1)) {
#       lp[[i]]<- tcrossprod(eta[,t]*eta[,(t-i)])/N
#     }
#       for(t in 1:(T-i)) {
#       lm[i] <- tcrossprod(eta[,t]*eta[,(t+i)])/N
#       }
#   }
#   RP<-add(lp)
#   RM<-add(lm)
#   
#   R_n<-RP+RM
#   psi<-(R0+R_n)/T
#   
# # #in progress-----------------------
#   tnsq = 2*N^2
#   
#   for(i in 1:N) {
#     psi = psi + exp(-i^2/tnsq)* (RP + RM)
#   }
#   
# #   R_n<-tcrossprod(eta,eta)/T
# #   psi<-R_n+R0
#   
#   for(i in 1:N) {
#     psi = psi + exp(-i^2/tnsq)* (lp[[i]]+lm[[i]])
#   }
#   #oude code
# #   totlist<-vector("list",length=T-3-s)
# #   l1 <-  vector("list", length=828-1-s)
# #   
# #   lag<-1
# #   for( lag in 1:825){
# #     for(i in (1+s+lag):828){
# #       l1[[(i-s-lag)]] <-tcrossprod(eta[,i],eta[,(i-lag)])
# #     }
# #     totlist[[lag]]<-(add(l1)/T)
# #   }
# # 
# #   #add <- function(x) Reduce("+", x)
# #   #R_1<- tcrossprod(eta[,(T-s)],eta[,(T-s-1)]) #kolom T is verwijderd
# #   
# #   psi<-add(totlist)
#   
#   
  
#  
	
	
  testmatrix<-phi%*%psi%*%t(phi)
  #r<-crossprod(res[2:T],res[1:(T-1)])/crossprod(res,res)
  r<-matrix(acf(res,lag.max=s,plot=FALSE,na.action=na.pass)$acf[2:(s+1)])

  
  I<-T*t(r)%*%solve(testmatrix[2,2]+testmatrix[1,1]+testmatrix[2,1]+testmatrix[1,2])%*%r

  
  result <- list()
  result[1]<-I
  result[2]<-1-pchisq(I,s)
  result
}

add <- function(x) Reduce("+", x)

