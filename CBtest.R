#######################################
########## Cumby-Huizinga #############
#######################################

#Cumby-Huizinga autocorrelation test
#Based on: Cumby, R. and Huizinga J., (1992) Econometrica; 
#ivactest Baum; Schaffer

CBtest<-function(x,s,data){
  require(ibdreg)
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
  
  for(j in 1:s){
    for(i in 1:T){
      if((i-j)>0) U[i,j]<-res[(i-j)] #q=0
    }
  }
  
  B<- -(crossprod(U,X)/T)/sigma_sq
  A_t<- crossprod(W,W)/T
  D<- T*Ginv(crossprod(X,W)%*%Ginv(A_t,eps=1e-25)$Ginv%*%crossprod(W,X),eps=1e-25)$Ginv%*% crossprod(X,W)%*%Ginv(A_t,eps=1e-25)$Ginv #bekijken want soms solv(A_t) niet mogelijk
  
#   D<- T*solve(crossprod(X,W)%*%solve(A_t)%*%crossprod(W,X))%*% crossprod(X,W)%*%solve(A_t) #bekijken want soms solv(A_t) niet mogelijk
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
	
  w_n<-matrix(nr=N)
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

