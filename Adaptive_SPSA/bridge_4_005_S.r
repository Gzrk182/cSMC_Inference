#!/usr/bin/env Rscript

library(MASS)
library(mvtnorm)
library(mvnfast)
#library(retry,lib.loc = "/storage/maths/strbdh/rlib")

######################stuff needed to run cSMC######################
log_sum_exp<-function(w){
  c<-max(w)
  return(c+log(sum(exp(w-c))))
}
extract_likelihood<-function(w,is.resamp,if.print=F){
  l<-log_sum_exp(w[1,])-log(dim(w)[2])
  for (i in 2:dim(w)[1]){
    if(is.resamp[i]==1){
      l_now<-log_sum_exp(w[i,])-log(dim(w)[2])
    }else{
      l_now<-log_sum_exp(w[i,])-log_sum_exp(w[(i-1),])
      
    }
    if(if.print==T){
      print(l_now)
    }
    l<-l+l_now
  }
  
  return(l)
}
#computes stably the log of sd of vector of tiny magnitudes
log_sd<-function(x){
  log_m_2<-log_sum_exp(2*x)-log(length(x))
  log_m_1<-log_sum_exp(x)-log(length(x))
  interm<-log_m_2-2*log_m_1
  log_var<-log(exp(interm)-1)+2*log_m_1
  return(0.5*log_var)
}
rel_sd<-function(L_vec,k){
  l<-length(L_vec)
  log_sd(L_vec[max((l-k+1),1):l])-(log_sum_exp(L_vec[max((l-k+1),1):l])-log(min(k,l)))
}

if_double_N<-function(L_vec,N_vec,k){
  l<-length(L_vec)
  if(l<k){
    return(F)
  }else{
    cond_1<-l>=k
    cond_2<-N_vec[l]==N_vec[max((l-k+1),1)]
    cond_3<-sum(diff(L_vec[max((l-k+1),1):l])>0)!=(k-1)
    return(cond_1&cond_2&cond_3)
  }
}
det_small<-function(x){
  
  if (length(x)==1) {
    return(x)
  }else{
    return(x[1,1]*x[2,2]-x[1,2]^2)}
}

inv_small<-function(x){
  if (length(x)==1) {
    return(1/x)
  }else{
    return((1 / (x[1,1]*x[2,2]-x[1,2]^2)) * matrix(c(x[2,2], -x[1,2], -x[1,2], x[1,1]), nrow = 2))
    
  }
}




scale_const_compute<-function(x,x_mean,x_cov,psi_mean,psi_cov){
  p<-length(psi_mean)
  mu_1<-x_mean(x)
  mu_2<-psi_mean
  mu_1m2<-mu_1-mu_2
  V_1<-x_cov(x)
  part_1<-log((2*pi)^p*(det_small(V_1+psi_cov)))*(-0.5)
  part_2<--0.5*t(mu_1m2)%*%inv_small(V_1+psi_cov)%*%mu_1m2
  #returns log of scaling constant
  return(part_1+part_2)
}


twist_f_specif<-function(x,x_mean,x_cov,psi_mean,psi_cov){
  
  p<-length(psi_mean)
  mu_1<-x_mean(x)
  mu_2<-psi_mean
  L_1<-inv_small(x_cov(x))
  L_2<-inv_small(psi_cov)
  V_n<-inv_small(L_1+L_2)
  mu_n<-V_n%*%(L_1%*%mu_1+L_2%*%mu_2)
  return(list(mu_n=mu_n,V_n=V_n))
}


twisted_g<-function(x_t,y_t,potential_function,psi_mean,psi_cov,scale_const){
  result<-potential_function(x_t,y_t)+scale_const-mvnfast::dmvn(x_t,psi_mean,psi_cov,log=T)
  #returns log of value of potential function
  return(result)
}

bpf<-function(N,y,tk_list,alpha){
  #browser()
  y<-as.matrix(y)
  tk_char<-c("x_mean","x_cov","potential_function")
  tk_now<-tk_list[[1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  T_max<-dim(y)[1]-1
  p_max<-p<-dim(y)[2]
  mvn_array<-rnorm(T_max*N*p_max)
  if(!is.null(attributes(potential_function))){
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }}
  x_0<-MASS::mvrnorm(n=N,x_mean(0),x_cov(rep(0)))
  y_0<-y[1,]
  w_0<-numeric(N)
  w<-matrix(nrow = T_max+1,ncol = N)
  W<-matrix(nrow = T_max+1,ncol = N)
  x<-vector(mode='list', length=T_max+1)
  x[[1]]<-x_0
  for (n in 1:N) {
    w_0[n]<-potential_function(x_0[n,],y_0)
  }
  is.resamp<-numeric(T_max+1)
  W_0<-exp(w_0-log_sum_exp(w_0))
  W[1,]<-W_0
  w[1,]<-w_0
  is.resamp[1]<-0
  
  for (t in 1:T_max) {
    if(1/sum(W[t,]^2)<alpha*N){
      A_t<-sample(1:N,N,replace = T, prob = W[t,])
      w_t_hat<-rep(0,N)
      is.resamp[t+1]<-1
    }else{
      A_t<-1:N
      w_t_hat<-w[t,]
      is.resamp[t+1]<-0
    }
    tk_now<-tk_list[[t+1]]
    for (i in 1:length(tk_char)) {
      assign(tk_char[i],tk_now[[i]])
    }
    if(!is.null(attributes(potential_function))){
      if(names(attributes(potential_function))=="p_now"){
        p<-attributes(potential_function)$p_now
      }}
    y_t<-y[t+1,]
    x_t<-matrix(nrow = N,ncol = p)
    w_t<-numeric(N)
    mvn_now<-matrix(mvn_array[((t-1)*N*p_max+1):(t*N*p_max)][1:(p*N)],nrow = N)
    for (n in 1:N) {
      mean_x<-x_mean(x[[t]][A_t[n],])
      cov_x<-x_cov(x[[t]][A_t[n],])
      x_t[n,]<-mean_x+t(chol(cov_x))%*%mvn_now[n,]
      w_t[n]<-w_t_hat[n]+potential_function(x_t[n,],y_t)
    }
    W_t<-exp(w_t-log_sum_exp(w_t))
    W[t+1,]<-W_t
    w[t+1,]<-w_t
    x[[t+1]]<-x_t
    #print(t)
  }
  return(list(x=x,w=w,is.resamp=is.resamp,idd=mean(1/rowSums(W^2))))
  
}

prod_cross<-function(x){
  p<-length(x)
  result<-vector()
  for (i in 1:p) {
    x_now<-x[1]
    result<-c(result,x_now*x)
    x<-x[-1]
  }
  return(result)
}

matrix_place<-function(x){
  n<-floor(sqrt(length(x)*2))
  result<-matrix(rep(0,n^2),nrow = n,ncol = n)
  result[lower.tri(result,diag = T)]<-x
  result<-result+t(result)
  return(result/2)
}

matrix_make<-function(x){
  N<-dim(x)[1]
  p<-dim(x)[2]
  quadr_block<-matrix(nrow = N,ncol = p*(p+1)/2)
  for (n in 1:N) {
    quadr_block[n,]<-prod_cross(x[n,])
  }
  return(cbind(quadr_block,x,rep(1,n)))
}

recu_f_approx_1<-function(x,y,tk_list){
  err_sum<-0
  y<-as.matrix(y)
  N<-dim(x[[1]])[1]
  T_max<-dim(y)[1]-1
  p<-dim(x[[1]])[2]
  psi_t<-numeric(N)
  
  tk_char<-c("x_mean","x_cov","potential_function")
  tk_now<-tk_list[[T_max+1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  if(!is.null(attributes(potential_function))){
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }}
  
  for (n in 1:N) {
    psi_t[n]<-potential_function(x[[T_max+1]][n,],y[T_max+1,])
  }
  Z<-matrix_make(as.matrix(x[[T_max+1]]))
  beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)
  
  beta_t<-beta_t[-length(beta_t)]
  Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)
  
  beta_t<-beta_t[-(1:(p*(p+1)/2))]
  mu_t<-Sigma_t%*%beta_t
  sig_eigen<-eigen(Sigma_t,symmetric = T)
  if(!all(sig_eigen$values>0)){
    a<-sig_eigen$values
    a[which(sig_eigen$value<=0)]<-1e-2
    Sigma_t<-sig_eigen$vectors%*%diag(a,nrow=p)%*%t(sig_eigen$vectors)
    err_sum<-err_sum+1
  }
  #Sigma_t<-diag(1e10,p)}
  opt_par_T_max<-list(mu_t,Sigma_t)
  psi_spec<-vector(mode = "list",1+T_max)
  psi_spec[[T_max+1]]<-opt_par_T_max
  
  #plot(x[[T_max+1]],psi_t)
  for(t in T_max:1){
    potential_function<-tk_list[[t]][[3]]
    x_mean<-tk_list[[t+1]][[1]]
    x_cov<-tk_list[[t+1]][[2]]
    if(!is.null(attributes(potential_function))){
      if(names(attributes(potential_function))=="p_now"){
        p<-attributes(potential_function)$p_now
      }}
    psi_t<-numeric(N)
    scale_const_vec<-numeric(N)
    den_vec<-numeric(N)
    
    for (n in 1:N) {
      scale_const_vec[n]<-scale_const_compute(x[[t]][n,],x_mean,x_cov,psi_spec[[t+1]][[1]],psi_spec[[t+1]][[2]])
      den_vec[n]<-potential_function(x[[t]][n,],y[t,])
    }
    
    psi_t<-den_vec+scale_const_vec
    
    
    Z<-matrix_make(x[[t]])
    beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)
    
    
    beta_t<-beta_t[-length(beta_t)]
    Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)
    
    beta_t<-beta_t[-(1:(p*(p+1)/2))]
    mu_t<-Sigma_t%*%beta_t
    sig_eigen<-eigen(Sigma_t,symmetric = T)
    if(!all(sig_eigen$values>0)){
      a<-sig_eigen$values
      a[which(sig_eigen$value<=0)]<-1e-2
      Sigma_t<-sig_eigen$vectors%*%diag(a,nrow=p)%*%t(sig_eigen$vectors)
      err_sum<-err_sum+1
    }
    opt_par_t<-list(mu_t,Sigma_t)
    psi_spec[[t]]<-opt_par_t
  }
  print(err_sum)
  return(psi_spec)
  
  #return(psi_T_max)
}

recu_f_approx_2<-function(x,y,tk_list){
  err_sum<-0
  y<-as.matrix(y)
  N<-dim(x[[1]])[1]
  T_max<-dim(y)[1]-1
  p<-dim(x[[1]])[2]
  psi_t<-numeric(N)
  
  tk_char<-c("x_mean","x_cov","potential_function")
  tk_now<-tk_list[[T_max+1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  if(names(attributes(potential_function))=="p_now"){
    p<-attributes(potential_function)$p_now
  }
  for (n in 1:N) {
    psi_t[n]<-potential_function(x[[T_max+1]][n,],y[T_max+1,])
  }
  Z<-matrix_make(as.matrix(x[[T_max+1]]))
  beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)
  
  beta_t<-beta_t[-length(beta_t)]
  Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)
  
  beta_t<-beta_t[-(1:(p*(p+1)/2))]
  mu_t<-Sigma_t%*%beta_t
  sig_eigen<-eigen(Sigma_t)
  if(!all(sig_eigen$values>0)){
    Sigma_t<-diag(1e5,p)
    err_sum<-err_sum+1
  }
  #Sigma_t<-diag(1e10,p)}
  opt_par_T_max<-list(mu_t,Sigma_t)
  psi_spec<-vector(mode = "list",1+T_max)
  psi_spec[[T_max+1]]<-opt_par_T_max
  
  #plot(x[[T_max+1]],psi_t)
  for(t in T_max:1){
    potential_function<-tk_list[[t]][[3]]
    x_mean<-tk_list[[t+1]][[1]]
    x_cov<-tk_list[[t+1]][[2]]
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }
    psi_t<-numeric(N)
    scale_const_vec<-numeric(N)
    den_vec<-numeric(N)
    
    for (n in 1:N) {
      scale_const_vec[n]<-scale_const_compute(x[[t]][n,],x_mean,x_cov,psi_spec[[t+1]][[1]],psi_spec[[t+1]][[2]])
      den_vec[n]<-potential_function(x[[t]][n,],y[t,])
    }
    
    psi_t<-den_vec+scale_const_vec
    
    
    Z<-matrix_make(x[[t]])
    beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)
    
    beta_t<-beta_t[-length(beta_t)]
    Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)
    
    beta_t<-beta_t[-(1:(p*(p+1)/2))]
    mu_t<-Sigma_t%*%beta_t
    sig_eigen<-eigen(Sigma_t)
    if(!all(sig_eigen$values>0)){
      Sigma_t<-diag(1e5,p)
      err_sum<-err_sum+1
    }
    opt_par_t<-list(mu_t,Sigma_t)
    psi_spec[[t]]<-opt_par_t
    #print(t)
  }
  print(err_sum)
  return(psi_spec)
  
  #return(psi_T_max)
}



twist_bpf<-function(N,y,tk_list,alpha,psi_spec,show_ESS=F){
  y<-as.matrix(y)
  T_max<-dim(y)[1]-1
  p_max<-p<-dim(y)[2]
  tk_char<-c("x_mean","x_cov","potential_function")
  tk_now<-tk_list[[1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  mvn_array<-rnorm(T_max*N*p_max)
  if(!is.null(attributes(potential_function))){
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }}
  f_0_tst<-twist_f_specif(rep(0,p),x_mean,x_cov,psi_spec[[1]][[1]],psi_spec[[1]][[2]])
  x_0_tst<-matrix(nrow = N,ncol = p)
  mvn_0<-MASS::mvrnorm(n=N,rep(0,p),diag(rep(1,p)))
  for (n in 1:N) {
    x_0_tst[n,]<-f_0_tst[[1]]+t(chol(f_0_tst[[2]]))%*%mvn_0[n,]
  }
  y_0<-y[1,]
  w_0_tst<-numeric(N)
  w_tst<-matrix(nrow = T_max+1,ncol = N)
  W_tst<-w_tst
  x_tst<-vector(mode='list', length=T_max+1)
  x_tst[[1]]<-x_0_tst
  psi_0_tilde<-scale_const_compute(rep(0,p),x_mean,x_cov,psi_spec[[1]][[1]],psi_spec[[1]][[2]])
  tk_char<-c("x_mean","x_cov")
  tk_now<-tk_list[[2]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  for (n in 1:N) {
    w_0_tst[n]<-twisted_g(x_0_tst[n,],y_0,potential_function,psi_spec[[1]][[1]],psi_spec[[1]][[2]],scale_const_compute(x_0_tst[n,],x_mean,x_cov,psi_spec[[2]][[1]],psi_spec[[2]][[2]]))+psi_0_tilde
  }
  W_0_tst<-exp(w_0_tst-log_sum_exp(w_0_tst))
  W_tst[1,]<-W_0_tst
  w_tst[1,]<-w_0_tst
  is.resamp<-numeric(T_max+1)
  is.resamp[1]<-0
  
  for (t in 1:(T_max-1)) {
    if(1/sum(W_tst[t,]^2)<alpha*N){
      A_t<-sample(1:N,N,replace = T, prob = W_tst[t,])
      w_t_tst_hat<-rep(0,N)
      is.resamp[t+1]<-1
    }else{
      A_t<-1:N
      w_t_tst_hat<-w_tst[t,]
      is.resamp[t+1]<-0
    }
    tk_char<-c("x_mean","x_cov","potential_function")
    tk_now<-tk_list[[t+1]]
    for (i in 1:length(tk_char)) {
      assign(tk_char[i],tk_now[[i]])
    }
    if(!is.null(attributes(potential_function))){
      if(names(attributes(potential_function))=="p_now"){
        p<-attributes(potential_function)$p_now
      }}
    y_t<-y[t+1,]
    x_t_tst<-matrix(nrow = N,ncol = p)
    w_t_tst<-numeric(N)
    mvn_now<-matrix(mvn_array[((t-1)*N*p_max+1):(t*N*p_max)][1:(p*N)],nrow = N)
    psi_t_tilde<-numeric(N)
    for (n in 1:N) {
      f_t_tst<-twist_f_specif(x_tst[[t]][A_t[n],],x_mean,x_cov,psi_spec[[t+1]][[1]],psi_spec[[t+1]][[2]])
      x_t_tst[n,]<-f_t_tst[[1]]+t(chol(f_t_tst[[2]]))%*%mvn_now[n,]
    }
    tk_char<-c("x_mean","x_cov")
    tk_now<-tk_list[[t+2]]
    for (i in 1:length(tk_char)) {
      assign(tk_char[i],tk_now[[i]])
    }
    for (n in 1:N) {
      psi_t_tilde[n]<-scale_const_compute(x_t_tst[n,],x_mean,x_cov,psi_spec[[t+2]][[1]],psi_spec[[t+2]][[2]])
      w_t_tst[n]<-w_t_tst_hat[n]+twisted_g(x_t_tst[n,],y_t,potential_function,psi_spec[[t+1]][[1]],
                                           psi_spec[[t+1]][[2]],psi_t_tilde[n])
    }
    W_t_tst<-exp(w_t_tst-log_sum_exp(w_t_tst))
    W_tst[t+1,]<-W_t_tst
    w_tst[t+1,]<-w_t_tst
    x_tst[[t+1]]<-x_t_tst
  }
  
  if(1/sum(W_tst[T_max,]^2)<alpha*N){
    A_t<-sample(1:N,N,replace = T, prob = W_tst[T_max,])
    w_t_tst_hat<-rep(0,N)
    is.resamp[T_max+1]<-1
  }else{
    A_t<-1:N
    w_t_tst_hat<-w_tst[T_max,]
    is.resamp[T_max+1]<-0
  }
  tk_char<-c("x_mean","x_cov","potential_function")
  tk_now<-tk_list[[T_max+1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  if(!is.null(attributes(potential_function))){
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }}
  y_t<-y[T_max+1,]
  x_t_tst<-matrix(nrow = N,ncol = p)
  w_t_tst<-numeric(N)
  mvn_now<-matrix(mvn_array[((T_max-1)*N*p_max+1):(T_max*N*p_max)][1:(p*N)],nrow = N)
  for (n in 1:N) {
    f_t_tst<-twist_f_specif(x_tst[[T_max]][A_t[n],],x_mean,x_cov,psi_spec[[T_max+1]][[1]],psi_spec[[T_max+1]][[2]])
    x_t_tst[n,]<-f_t_tst[[1]]+t(chol(f_t_tst[[2]]))%*%mvn_now[n,]
    w_t_tst[n]<-w_t_tst_hat[n]+twisted_g(x_t_tst[n,],y_t,potential_function,psi_spec[[T_max+1]][[1]],psi_spec[[T_max+1]][[2]],0)
  }
  W_t_tst<-exp(w_t_tst-log_sum_exp(w_t_tst))
  W_tst[T_max+1,]<-W_t_tst
  w_tst[T_max+1,]<-w_t_tst
  x_tst[[T_max+1]]<-x_t_tst
  if(show_ESS==T){
    print(1/rowSums(W_tst^2))
    return(list(x_tst=x_tst,w_tst=w_tst,is.resamp=is.resamp,ESS=1/rowSums(W_tst^2)))
  }else{
    return(list(x_tst=x_tst,w_tst=w_tst,is.resamp=is.resamp,idd=mean(1/rowSums(W_tst^2))))
  }
  
}



cSMC<-function(N,y,tk_list,k=3,alpha=0.5,if_print=F,keep_psi_spec=F){
  #initialize with a bpf
  pf_result<-bpf(N,y,tk_list,alpha)
  L_bpf<-extract_likelihood(pf_result[[2]],pf_result[[3]])
  l<-0
  L_vec<-vector()
  n_1<-c("Iteration","N","L","#Resamp")
  n_2<-c(0,N,L_bpf,sum(pf_result[[3]]))
  if(if_print==T){
    print(paste(n_1,n_2,sep = "="))
  }
  I<-0
  keep_vec<-vector()
  while ((l<k)){
    I<-I+1
    
    psi_spec<-recu_f_approx_2(pf_result[[1]],y,tk_list)
    pf_result<-twist_bpf(N,y,tk_list,alpha,psi_spec)
    L_vec<-c(L_vec,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    if(sum(pf_result[[3]])==0){
      l<-l+1
      keep_vec<-c(keep_vec,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    }
    if(I>1000){
      l<-l+1
      keep_vec<-c(keep_vec,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    }
    n_1<-c("Iteration","N","L","#Resamp")
    n_1<-c("Iteration","#Resamp","L","#MeanESS")
    n_2<-c(l,N,L_vec[length(L_vec)],sum(pf_result[[3]]))
    n_2<-c(I,sum(pf_result[[3]]),L_vec[length(L_vec)],pf_result[[4]])
    if(if_print==T){
      print(paste(n_1,n_2,sep = "="))
    }
  }
  if(keep_psi_spec==F){return(mean(keep_vec))
  }else
  {
    return(psi_spec)
  }
  
}
######################END##########################


#################FHN simulator using Strang#####################
f_t<-function(t,e,g,b,s_1,s_2,X_0){
  X_0<-as.numeric(X_0)
  a<--2*t/e
  comp_1<-X_0[1]/sqrt(exp(a)+X_0[1]^2*(1-exp(a)))
  comp_2<-b*t+X_0[2]
  return(c(comp_1,comp_2))
}
f_t_1<-function(t,e,g,b,s_1,s_2,v_0){
  v_0<-as.numeric(v_0)
  a<--2*t/e
  comp_1<-v_0/sqrt(exp(a)+v_0^2*(1-exp(a)))
  return(comp_1)
}
f_t_2<-function(t,e,g,b,s_1,s_2,u_0){
  u_0<-as.numeric(u_0)
  comp_2<-b*t+u_0
  return(comp_2)
}
f_inv<-function(X,t,e){
  a<-exp(-2*t/e)
  numtr<-a*X^2
  denom<-1-(1-a)*X^2
  if(sum(denom<0)!=0){
    browser()
  }
  return(sign(X)*sqrt(numtr/denom))
}
matrix_eAt<-function(t,e,g,b,s_1,s_2){
  k<-4*g/e-1
  if(k==0){
    return(exp(-t/2)*matrix(c(1+t/2,e*t/4,-t/e,1-t/2),nrow = 2))
  }else if(k>0){
    c_t<-cos(0.5*sqrt(k)*t)
    s_t<-sin(0.5*sqrt(k)*t)/sqrt(k)
    return(exp(-t/2)*matrix(c(c_t+s_t,2*g*s_t,-2*s_t/(e),c_t-s_t),nrow = 2))
  }else if(k<0){
    c_t<-cosh(0.5*sqrt(-k)*t)
    s_t<-sinh(0.5*sqrt(-k)*t)
    return(exp(-t/2)*matrix(c(c_t+s_t,2*g*s_t,-2*s_t/(e),c_t-s_t),nrow = 2))
  }
}


matrix_Ct<-function(t,e,g,b,s_1,s_2){
  k<-4*g/e-1
  if(k==0){
    c_11<-exp(-t)/(4*e^2)*(4*s_2^2*(-2+2*exp(t)-t*(2+t))+e^2*s_1^2*(-10+10*exp(t)-t*(6+t)))
    c_12<-exp(-t)/(8*e)*(-4*s_2^2*t^2+e^2*s_1^2*(4*exp(t)-(2+t)^2))
    c_21<-c_12
    c_22<-exp(-t)/16*(4*s_2^2*(-2+2*exp(t)-(t-2)*t)+e^2*s_1^2*(-2+2*exp(t)-t*(2+t)))
    return(matrix(c(c_11,c_21,c_12,c_22),nrow = 2))
  }else if(k>0){
    c_t<-cos(sqrt(k)*t)
    s_t<-sin(sqrt(k)*t)/sqrt(k)
    c_11<-e*exp(-t)/(2*g*k)*(-4*g/e^2*(s_1^2*g+s_2^2/e)+k*exp(t)*(s_1^2*(1+g/e)+s_2^2/e^2)
                             +c_t*(s_1^2*(1-3*g/e)+s_2^2/e^2)-k*s_t*(s_1^2*(1-g/e)+s_2^2/e^2))
    c_12<-e*exp(-t)/(2*k)*(s_1^2*k*exp(t)-2/e*(s_1^2*g+s_2^2/e)
                           +c_t*(s_1^2*(1-2*g/e)+2*s_2^2/e^2)-k*s_t*s_1^2)
    c_21<-c_12
    c_22<-e*exp(-t)/(2*k)*((s_2^2/e+s_1^2*g)*(c_t-4*g/e+k*exp(t))+k*s_t*(s_2^2/e-s_1^2*g))
    return(matrix(c(c_11,c_21,c_12,c_22),nrow = 2))
  }else if(k<0){
    c_t<-cosh(sqrt(-k)*t)
    s_t<-sinh(sqrt(-k)*t)
    c_11<-e*exp(-t)/(2*g*k)*(-4*g/e^2*(s_1^2*g+s_2^2/e)+k*exp(t)*(s_1^2*(1+g/e)+s_2^2/e^2)
                             +c_t*(s_1^2*(1-3*g/e)+s_2^2/e^2)-k*s_t*(s_1^2*(1-g/e)+s_2^2/e^2))
    c_12<-e*exp(-t)/(2*k)*(s_1^2*k*exp(t)-2/e*(s_1^2*g+s_2^2/e)
                           +c_t*(s_1^2*(1-2*g/e)+2*s_2^2/e^2)-k*s_t*s_1^2)
    c_21<-c_12
    c_22<-e*exp(-t)/(2*k)*((s_2^2/e+s_1^2*g)*(c_t-4*g/e+k*exp(t))+k*s_t*(s_2^2/e-s_1^2*g))
    return(matrix(c(c_11,c_21,c_12,c_22),nrow = 2))
  }
}

FHN_S_Gnr<-function(T_max=10,T_step=10^-4,e=1,g=20,b=0.1,s_1=0.1,s_2=0.2,X_0=c(0,0)){
  T_total<-ceiling(T_max/T_step)
  X_matrix<-matrix(nrow = 2,ncol = 1+T_total)
  X_matrix[,1]<-X_0
  Z_matrix<-matrix(nrow = 2,ncol = 1+T_total)
  Z_matrix[,1]<-X_0
  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  mvn_array<-MASS::mvrnorm(T_total,rep(0,2),Cdelta)
  for (i in 1:T_total) {
    z<-eAdelta%*%f_t(T_step/2,e,g,b,s_1,s_2,X_matrix[,i])+mvn_array[i,]
    Z_matrix[,i+1]<-z
    X_matrix[,i+1]<-f_t(T_step/2,e,g,b,s_1,s_2,z)
    
  }
  return(list(X=data.frame(t(X_matrix)),Z=data.frame(t(Z_matrix))))
}


############################COPYEVERYTIMEEND#################################

######################stuff needed to run the hypoelliptic Lie-Trotter bridge k=10 result######################

list_extrop<-function(l_1,l_2){return(list((l_1[[1]]+l_2[[1]])/2,(l_1[[2]]+l_2[[2]])/2))}

FHN_S_bridge_2_list_gen<-function(V, par_vec, T_step){
  T_step_0<-T_step
  #initiallize
  par_name<-c("e","g","b","s_1","s_2")
  for (i in 1:5) {
    assign(par_name[i],par_vec[i])
  }
  
  
  T_max<-length(V)-1
  
  V_0<-f_inv(V,T_step_0/2,e)[1]
  a_0<-exp(-T_step_0/e)
  f_0<-expression(z_t/sqrt(a_0+z_t^2*(1-a_0)))
  Df_0<-D(f_0,name="z_t")
  
  eAdelta_0<-matrix_eAt(T_step_0,e,g,b,s_1,s_2)
  Cdelta_0<-matrix_Ct(T_step_0,e,g,b,s_1,s_2)
  mu_a_0<-eAdelta_0[1,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
    eAdelta_0[1,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0)
  mu_b_0<-eAdelta_0[2,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
    eAdelta_0[2,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0)
  tk_list<-vector(mode = "list",1+2*T_max)
  x_ini_mean<-function(x){
    return(mu_b_0+Cdelta_0[2,1]*(V[1]-mu_a_0)/Cdelta_0[1,1])
  }
  x_ini_cov<-function(x){
    return(Cdelta_0[2,2]-Cdelta_0[2,1]^2/Cdelta_0[1,1])
  }
  potential_0<-function(x,y){
    z_t<-V_0
    return(mvnfast::dmvn(V_0,eAdelta_0[1,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
                           eAdelta_0[1,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0),
                         Cdelta_0[1,1],log = T)-log(abs(eval(expr = Df_0))))
  }
  attributes(potential_0)<-list(p_now=1)
  tk_list[[1]]<-list(x_ini_mean,x_ini_cov,potential_0)
  #shrink the time scale
  T_step<-T_step/2
  V<-f_inv(V,T_step/2,e)
  a<-exp(-T_step/e)
  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  
  ########
  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  x_mean_ux<-function(V){
    force(V)
    return(
      function(x){
        return(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,c(V,x)))
      }
    )
  }
  x_cov_ux<-function(x){
    return(Cdelta)
  }
  potential_ux<-function(x,y){
    return(0)
  }
  attributes(potential_ux)<-list(p_now=2)
  x_mean_xx<-function(x){
    return(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x))
  }
  x_cov_xx<-function(x){
    return(Cdelta)
  }
  potential_xx<-function(x,y){
    return(0)
  }
  attributes(potential_xx)<-list(p_now=2)
  potential_xx_last<-function(V){
    force(V)
    potential<-function(x,y){
      z_t<-V
      mu<-(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x))[1]
      return(mvnfast::dmvn(V,mu,Cdelta[1,1],log = T)-log(abs(eval(expr = Df))))
    }
    attributes(potential)<-list(p_now=2)
    return(potential)
  }
  x_mean_xu<-function(V){
    force(V)
    return(
      function(x){
        meanx<-eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x)
        mu_b<-meanx[2]
        mu_a<-meanx[1]
        return(mu_b+Cdelta[2,1]*(V-mu_a)/Cdelta[1,1])
      }
    )
  }
  x_cov_xu<-function(x){
    return(Cdelta[2,2]-Cdelta[2,1]^2/Cdelta[1,1])
  }
  potential_xu<-function(x,y){
    return(0)
  }
  attributes(potential_xu)<-list(p_now=1)
  for (t in 2:(1+T_max)) {
    tk_list[[(t-2)*2+2]]<-list(x_mean_ux(V[t-1]),x_cov_ux,potential_xx_last(V[t]))
    tk_list[[(t-1)*2+1]]<-list(x_mean_xu(V[t]),x_cov_xu,potential_xu)
  }
  return(tk_list[-length(tk_list)])
}


FHN_S_bridge_list_gen<-function(V, par_vec, T_step,K){
  T_step_0<-T_step
  #initiallize
  par_name<-c("e","g","b","s_1","s_2")
  for (i in 1:5) {
    assign(par_name[i],par_vec[i])
  }
  
  T_max<-length(V)-1
  V_0<-f_inv(V,T_step_0/2,e)[1]
  a_0<-exp(-T_step_0/e)
  f_0<-expression(z_t/sqrt(a_0+z_t^2*(1-a_0)))
  Df_0<-D(f_0,name="z_t")
  eAdelta_0<-matrix_eAt(T_step_0,e,g,b,s_1,s_2)
  Cdelta_0<-matrix_Ct(T_step_0,e,g,b,s_1,s_2)
  mu_a_0<-eAdelta_0[1,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
    eAdelta_0[1,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0)
  mu_b_0<-eAdelta_0[2,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
    eAdelta_0[2,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0)
  tk_list<-vector(mode = "list",1+K*T_max)
  x_ini_mean<-function(x){
    return(mu_b_0+Cdelta_0[2,1]*(V[1]-mu_a_0)/Cdelta_0[1,1])
  }
  x_ini_cov<-function(x){
    return(Cdelta_0[2,2]-Cdelta_0[2,1]^2/Cdelta_0[1,1])
  }
  potential_0<-function(x,y){
    z_t<-V_0
    return(mvnfast::dmvn(V_0,eAdelta_0[1,1]*f_t_1(T_step_0,e,g,b,s_1,s_2,0)+
                           eAdelta_0[1,2]*f_t_2(T_step_0,e,g,b,s_1,s_2,0),
                         Cdelta_0[1,1],log = T)-log(abs(eval(expr = Df_0))))
  }
  attributes(potential_0)<-list(p_now=1)
  tk_list[[1]]<-list(x_ini_mean,x_ini_cov,potential_0)
  #shrink the time scale
  T_step<-T_step/(K)
  V<-f_inv(V,T_step/2,e)
  a<-exp(-T_step/e)
  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  ########
  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  x_mean_ux<-function(V){
    force(V)
    return(
      function(x){
        return(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,c(V,x)))
      }
    )
  }
  x_cov_ux<-function(x){
    return(Cdelta)
  }
  potential_ux<-function(x,y){
    return(0)
  }
  attributes(potential_ux)<-list(p_now=2)
  x_mean_xx<-function(x){
    return(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x))
  }
  x_cov_xx<-function(x){
    return(Cdelta)
  }
  potential_xx<-function(x,y){
    return(0)
  }
  attributes(potential_xx)<-list(p_now=2)
  potential_xx_last<-function(V){
    force(V)
    potential<-function(x,y){
      z_t<-V
      mu<-(eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x))[1]
      return(mvnfast::dmvn(V,mu,Cdelta[1,1],log = T)-log(abs(eval(expr = Df))))
    }
    attributes(potential)<-list(p_now=2)
    return(potential)
  }
  x_mean_xu<-function(V){
    force(V)
    return(
      function(x){
        meanx<-eAdelta%*%f_t(T_step,e,g,b,s_1,s_2,x)
        mu_b<-meanx[2]
        mu_a<-meanx[1]
        return(mu_b+Cdelta[2,1]*(V-mu_a)/Cdelta[1,1])
      }
    )
  }
  x_cov_xu<-function(x){
    return(Cdelta[2,2]-Cdelta[2,1]^2/Cdelta[1,1])
  }
  potential_xu<-function(x,y){
    return(0)
  }
  attributes(potential_xu)<-list(p_now=1)
  for (t in 2:(1+T_max)) {
    tk_list[[(t-2)*K+2]]<-list(x_mean_ux(V[t-1]),x_cov_ux,potential_ux)
    for (k in 1:(K-3)) {
      tk_list[[(t-2)*K+2+k]]<-list(x_mean_xx,x_cov_xx,potential_xx)
    }
    tk_list[[(t-1)*K]]<-list(x_mean_xx,x_cov_xx,potential_xx_last(V[t]))
    tk_list[[(t-1)*K+1]]<-list(x_mean_xu(V[t]),x_cov_xu,potential_xu)
  }
  return(tk_list[-length(tk_list)])
}



cSMC_4_slow<-function(tk_list_2,tk_list_4,y_2,y_4){
  psi_spec_0<-cSMC(10,y_2,tk_list_2,0.5,if_print = T,k=3,keep_psi_spec = T)
  #step 2: create initialization psi_spec for n_bridge=4
  psi_spec_1<-vector(mode = "list", 4*999)
  psi_spec_1[[1]]<-psi_spec_0[[1]]
  for (t in 2:(999)) {
    psi_spec_1[[(t-2)*4+2]]<-psi_spec_0[[(t-2)*2+2]]
    
    #psi_spec_0[[(t-2)*2+2]][[1]]+c(psi_spec_0[[(t-2)*2+1]][[1]],V[-1][t-1])
    
    psi_spec_1[[(t-2)*4+3]]<-list(rep(0,2),diag(1e10,2))
    
    psi_spec_1[[(t-1)*4]]<-list(rep(0,2),diag(1e10,2))
    psi_spec_1[[(t-1)*4+1]]<-psi_spec_0[[(t-1)*2+1]]
  }
  t<-1000
  psi_spec_1[[(t-2)*4+2]]<-psi_spec_0[[(t-2)*2+2]]
  psi_spec_1[[(t-2)*4+3]]<-list(rep(0,2),diag(1e10,2))
  psi_spec_1[[(t-1)*4]]<-list(rep(0,2),diag(1e10,2))
  
  #run with n_bridge=4
  
  #cook the new policies
  pf_result_1<-twist_bpf(100,y_4,tk_list_4,0.5,psi_spec_1)
  #pf_result_1<-bpf(1000,y_4,tk_list_4,0.5)
  #better now
  extract_likelihood(pf_result_1[[2]],pf_result_1[[3]])
  sum(pf_result_1[[3]])
  psi_spec_2<-recu_f_approx_1(pf_result_1[[1]],y_4,tk_list_4)
  #visible improvements
  pf_result_2<-twist_bpf(50,y_4,tk_list_4,0.5,psi_spec_2)
  extract_likelihood(pf_result_2[[2]],pf_result_2[[3]])
  sum(pf_result_2[[3]])
  
  #keep doing this
  pf_result<-pf_result_2
  ggg<-0
  hhh<-0
  ess<-pf_result[[4]]/50
  l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
  while (ggg<3) {
    if(ess>0.975&l_est<(-1e5)){
      psi_spec<-recu_f_approx_2(pf_result[[1]],y_4,tk_list_4)
    }else{
      psi_spec<-recu_f_approx_1(pf_result[[1]],y_4,tk_list_4)
    }
    
    pf_result<-twist_bpf(50,y_4,tk_list_4,0.5,psi_spec)
    if(sum(pf_result[[3]])==0){
      ggg<-ggg+1
      hhh<-hhh+extract_likelihood(pf_result[[2]],pf_result[[3]])
    }
    ess<-pf_result[[4]]/50
    l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
    print(c(sum(pf_result[[3]]),extract_likelihood(pf_result[[2]],pf_result[[3]]),pf_result[[4]]))
  }
  return(list(hhh/3,psi_spec))
}

cSMC_4<-function(psi_spec_ini,tk_list_4,y_4,k){
  
  #run with n_bridge=4
  
  #cook the new policies
  pf_result_1<-twist_bpf(20,y_4,tk_list_4,0.5,psi_spec_ini)
  #pf_result_1<-bpf(1000,y_4,tk_list_4,0.5)
  #better now
  #keep doing this
  pf_result<-pf_result_1
  ggg<-0
  hhh<-vector()
  ess<-pf_result[[4]]/20
  l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
  while (ggg<k) {
    if(ess>0.975&l_est<(-1e5)){
      psi_spec<-recu_f_approx_2(pf_result[[1]],y_4,tk_list_4)
    }else{
      psi_spec<-recu_f_approx_1(pf_result[[1]],y_4,tk_list_4)
    }
    pf_result<-twist_bpf(20,y_4,tk_list_4,0.5,psi_spec)
    if(sum(pf_result[[3]])==0){
      ggg<-ggg+1
      hhh<-c(hhh,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    }
    ess<-pf_result[[4]]/20
    l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
    print(c(sum(pf_result[[3]]),extract_likelihood(pf_result[[2]],pf_result[[3]]),pf_result[[4]]))
  }
  #pf_100<-twist_bpf(100,y_4,tk_list_4,0.5,psi_spec)
  return(list(median(hhh),psi_spec))
}
cSMC_8<-function(psi_spec_ini,tk_list_8,y_8,k){
  
  #run with n_bridge=4
  
  #cook the new policies
  pf_result_1<-twist_bpf(20,y_8,tk_list_8,0.5,psi_spec_ini)
  #keep doing this
  pf_result<-pf_result_1
  ggg<-0
  hhh<-vector()
  ess<-pf_result[[4]]/20
  l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
  I<-0
  while (ggg<k&I<20) {
    I<-I+1
    if(ess>0.975&l_est<(-1e5)){
      psi_spec<-recu_f_approx_2(pf_result[[1]],y_8,tk_list_8)
    }else{
      psi_spec<-recu_f_approx_1(pf_result[[1]],y_8,tk_list_8)
    }
    pf_result<-twist_bpf(20,y_8,tk_list_8,0.5,psi_spec)
    if(sum(pf_result[[3]])==0){
      ggg<-ggg+1
      hhh<-c(hhh,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    }
    ess<-pf_result[[4]]/20
    l_est<-extract_likelihood(pf_result[[2]],pf_result[[3]])
    print(c(sum(pf_result[[3]]),extract_likelihood(pf_result[[2]],pf_result[[3]]),pf_result[[4]]))
  }
  pf_out<-twist_bpf(20,y_8,tk_list_8,0.5,psi_spec)
  l_out<-extract_likelihood(pf_out[[2]],pf_out[[3]])
  return(list(l_out,psi_spec,I))
}
vec_div<-function(x_1,x_2){
  result<-matrix(1:(length(x_1)^2),nrow = length(x_1))
  for (i in 1:length(x_1)) {
    result[,i]<-x_1[i]/x_2
  }
  return(result)
}
force_spd<-function(x){
  p<-dim(x)[1]
  x_sq<-x%*%x
  Sigma_t<-x_sq
  sig_eigen<-eigen(Sigma_t)
  a<-sig_eigen$values
  a[which(sig_eigen$value<=1e-20)]<-1e-2
  Sigma_t<-sig_eigen$vectors%*%diag(sqrt(a),nrow=p)%*%t(sig_eigen$vectors)
  return(Sigma_t)
}
seed<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
a_vec<-1*(51:1e4)^(-0.602)
c_vec<-5e-3*(1:1e4)^(-0.101)
d_vec<-1e-2*(1:1e4)^(-0.101)
par_0<-(c(0.09,1.4,0.7,0.35))

K<-length(a_vec)
n_bridge<-4
T_step<-5e-2
T_max<-50;
T_step_sim<-5e-3
set.seed(seed)
FHN_S_obs<-FHN_S_Gnr(T_max,T_step_sim,0.1,1.5,0.8,0,0.3)
set.seed(NULL)
V_20_e<-FHN_S_obs$X$X1[1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
par_vec<-matrix(1:(5*(1+K)),ncol=5)
scale_par<-c(0.1,1.5,0.8,0.3)*c(1,2,2,1)
#scale_par<-c(1,1,1,1)
ini_par<-par_0/scale_par
par_vec[1,]<-c(ini_par*scale_par,0)
y_2<-matrix(1:((2*(999))*2),ncol=2)
y_4<-matrix(1:((4*(999))*2),ncol=2)

tk_list_2<-FHN_S_bridge_2_list_gen(V_20_e[-1],append(par_vec[1,1:4],
                                                      0,after = 3),T_step)
tk_list_4<-FHN_S_bridge_list_gen(V_20_e[-1],append(par_vec[1,1:4],
                                                    0,after = 3),T_step,4)
y_m_f<-cSMC_4_slow(tk_list_2,tk_list_4,y_2,y_4)

h_k_old<-0
i<-1
while (i<200) {
  
  theta_old<-par_vec[i,1:4]
  #create initialization psi_spec for n_bridge=4
  psi_spec_ini<-y_m_f[[2]]
  
  
  Delta_now<-sample(c(1,-1),4,replace = T,prob = rep(0.5,2))
  
  y_p_f<-cSMC_4(psi_spec_ini,FHN_S_bridge_list_gen(V_20_e[-1],append((theta_old/scale_par+c_vec[i]*Delta_now)*scale_par,
                                                                      0,after = 3),T_step,n_bridge),y_4,2)
  psi_spec_ini<-y_p_f[[2]]
  y_m_f<-cSMC_4(psi_spec_ini,FHN_S_bridge_list_gen(V_20_e[-1],append((theta_old/scale_par-c_vec[i]*Delta_now)*scale_par,
                                                                      0,after = 3),T_step,n_bridge),y_4,2)
  
  delta_now<-sample(c(1,-1),4,replace = T,prob = rep(0.5,2))
  
  psi_spec_ini<-y_p_f[[2]]
  y_p_f_d<-cSMC_4(psi_spec_ini,FHN_S_bridge_list_gen(V_20_e[-1],append((theta_old/scale_par+d_vec[i]*delta_now+c_vec[i]*Delta_now)*scale_par,
                                                                        0,after = 3),T_step,n_bridge),y_4,2)
  psi_spec_ini<-y_p_f_d[[2]]
  y_m_f_d<-cSMC_4(psi_spec_ini,FHN_S_bridge_list_gen(V_20_e[-1],append((theta_old/scale_par+d_vec[i]*delta_now-c_vec[i]*Delta_now)*scale_par,
                                                                        0,after = 3),T_step,n_bridge),y_4,2)
  
  y_p<-y_p_f[[1]]
  y_m<-y_m_f[[1]]
  y_p_d<-y_p_f_d[[1]]
  y_m_d<-y_m_f_d[[1]]
  g_1<-(y_p_d-y_p)/(d_vec[i]*delta_now)
  g_2<-(y_m_d-y_m)/(d_vec[i]*delta_now)
  g_d<-g_1-g_2
  g_m<-vec_div(g_d,2*c_vec[i]*Delta_now)
  h_k_now<-((i-1)/(i))*h_k_old+(1/(i))*0.5*(g_m+t(g_m))
  g_k<-(y_p-y_m)/(2*c_vec[i]*Delta_now)
  h_k_final<-force_spd(h_k_now)
  gain<-(a_vec[i]*solve(h_k_final,tol = 1e-50)%*%g_k)
  
  if(max(abs(gain))>0.1){
    gain<-gain/(10*max(abs(gain)))
  }
  print(gain)
  par_new<-((theta_old/scale_par+gain)*scale_par)
  l_new<-cSMC_4(psi_spec_ini,FHN_S_bridge_list_gen(V_20_e[-1],append(par_new,
                                                                      0,after = 3),T_step,n_bridge),y_4,2)[[1]]
  if(l_new>((y_p+y_m)/2-0.1)){
    par_vec[i+1,1:4]<-par_new
    par_vec[i+1,5]<-l_new
    
  }else{
    par_vec[i+1,1:4]<-par_vec[i,1:4]
    par_vec[i+1,5]<-par_vec[i,5]
    print("blocked")
  }
  print(c((par_vec[i+1,1:4]),l_new,i))
  h_k_old<-h_k_now
  i<-i+1
}
par_out<-par_vec[1:i,]
save(par_out,file = paste("par_",seed,".rda",sep = ""))






