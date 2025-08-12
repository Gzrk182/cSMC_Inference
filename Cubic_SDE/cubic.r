library(doParallel)
library(MASS)
library(mvtnorm)
library(mvnfast)
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

scale_const_compute<-function(x,x_mean,x_cov,psi_mean,psi_cov){
  p<-length(psi_mean)
  mu_1<-x_mean(x)
  mu_2<-psi_mean
  mu_1m2<-mu_1-mu_2
  V_1<-x_cov(x)
  part_1<-log((2*pi)^p*(det(V_1+psi_cov)))*(-0.5)
  part_2<--0.5*t(mu_1m2)%*%chol2inv(chol(V_1+psi_cov))%*%mu_1m2
  #returns log of scaling constant
  return(as.matrix(part_1+part_2))
}

twist_f_specif<-function(x,x_mean,x_cov,psi_mean,psi_cov){
  p<-length(psi_mean)
  mu_1<-x_mean(x)
  mu_2<-psi_mean
  L_1<-chol2inv(chol(x_cov(x)))
  L_2<-chol2inv(chol(psi_cov))
  V_n<-chol2inv(chol(L_1+L_2))
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

  beta_matrix<-beta_t
  beta_t<-beta_t[-length(beta_t)]
  Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)

  beta_t<-beta_t[-(1:(p*(p+1)/2))]
  mu_t<-Sigma_t%*%beta_t
  sig_eigen<-eigen(Sigma_t)
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
  psi_matrix<-psi_t
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
    psi_matrix<-rbind(psi_matrix,psi_t)

    Z<-matrix_make(x[[t]])
    beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)

    beta_matrix<-rbind(beta_matrix,beta_t)
    beta_t<-beta_t[-length(beta_t)]
    Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]),tol = 1e-100)

    beta_t<-beta_t[-(1:(p*(p+1)/2))]
    mu_t<-Sigma_t%*%beta_t
    sig_eigen<-eigen(Sigma_t)
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
  beta_matrix<-beta_t
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
  psi_matrix<-psi_t
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
    psi_matrix<-rbind(psi_matrix,psi_t)

    Z<-matrix_make(x[[t]])
    beta_t<-qr.solve(t(Z)%*%Z,t(Z)%*%psi_t,tol=1e-100)
    beta_matrix<-rbind(beta_matrix,beta_t)
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
  W_tst<-matrix(nrow = T_max+1,ncol = N)
  x_tst<-vector(mode='list', length=T_max+1)
  x_tst[[1]]<-x_0_tst
  psi_0_tilde<-scale_const_compute(rep(0,p),x_mean,x_cov,psi_spec[[1]][[1]],psi_spec[[1]][[2]])
  tk_char<-c("x_mean","x_cov")
  tk_now<-tk_list[[2]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  for (n in 1:N) {

    psi_1_tilde<-scale_const_compute(x_0_tst[n,],x_mean,x_cov,psi_spec[[2]][[1]],psi_spec[[2]][[2]])
    w_0_tst[n]<-twisted_g(x_0_tst[n,],y_0,potential_function,psi_spec[[1]][[1]],psi_spec[[1]][[2]],psi_1_tilde)+psi_0_tilde
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
  if(keep_psi_spec==F){pf_result<-twist_bpf(N,y,tk_list,alpha,psi_spec)

  return(extract_likelihood(pf_result[[2]],pf_result[[3]]))
  }else
  {
    return(psi_spec)
  }

}
cubic_sampler_strang<-function(T_step,T_max,sigma){
  N<-ceiling(T_max/T_step)
  x_vec<-numeric(N+1)
  x_vec[1]<-0
  normal_array<-rnorm(N)
  z_vec<-numeric(N+1)
  z_vec[1]<-0
  xi_sd<-sqrt((sigma^2/2)*(1-exp(-2*T_step)))
  for (t in 1:N) {
    x_hat<-(x_vec[t])/sqrt(exp(-T_step)+x_vec[t]^2*(1-exp(-T_step)))
    x_tilde<-exp(-T_step)*x_hat+xi_sd*normal_array[t]
    x_vec[t+1]<-(x_tilde)/sqrt(exp(-T_step)+x_tilde^2*(1-exp(-T_step)))
    z_vec[t+1]<-x_tilde
  }
  return(list(x=x_vec,z=z_vec))
}

cubic_sampler_strang_ini<-function(T_step,T_max,sigma,ini){
  N<-ceiling(T_max/T_step)
  x_vec<-numeric(N+1)
  x_vec[1]<-ini
  normal_array<-rnorm(N)
  z_vec<-numeric(N+1)
  z_vec[1]<-ini
  xi_sd<-sqrt((sigma^2/2)*(1-exp(-2*T_step)))
  for (t in 1:N) {
    x_hat<-(x_vec[t])/sqrt(exp(-T_step)+x_vec[t]^2*(1-exp(-T_step)))
    x_tilde<-exp(-T_step)*x_hat+xi_sd*normal_array[t]
    x_vec[t+1]<-(x_tilde)/sqrt(exp(-T_step)+x_tilde^2*(1-exp(-T_step)))
    z_vec[t+1]<-x_tilde
  }
  return(list(x=x_vec,z=z_vec))
}


lh_LT<-function(data_obs,T_step,sigma){
  xi_var<-(sigma^2/2)*(1-exp(-2*T_step))
  l_t<-dnorm(data_obs[1],0,sqrt(xi_var),log = T)
  for (t in 2:length(data_obs)) {
    x<-data_obs[t-1]
    x_hat<-(x)/sqrt(exp(-2*T_step)+x^2*(1-exp(-2*T_step)))
    l_t<-c(l_t,dnorm(data_obs[t],exp(-T_step)*x_hat,sqrt(xi_var),log = T))
  }
  return(l_t)
}

lh_S<-function(z,T_step,sigma){
  f<-expression(x/sqrt(exp(-T_step)+x^2*(1-exp(-T_step))))
  Df<-D(f,name="x")
  xi_var<-(sigma^2/2)*(1-exp(-2*T_step))
  x<-z[1]
  l_t<-dnorm(z[1],0,sqrt(xi_var),log = T)-log(abs(eval(expr = Df)))
  for (t in 2:length(z)) {
    x_hat<-(z[t-1])/sqrt(exp(-2*T_step)+z[t-1]^2*(1-exp(-2*T_step)))
    x<-z[t]
    l_t<-c(l_t,dnorm(z[t],exp(-T_step)*x_hat,sqrt(xi_var),log = T)-log(abs(eval(expr = Df))))
  }
  return(l_t)
}

cubic_bridge_EM_list_gen<-function(T_step,x_start,x_end,sigma,k){
  force(x_start);force(x_end)
  tk_list<-vector(mode = "list",k-1)
  T_step<-T_step/k
  x_ini_mean<-function(x){
    return(x_start-x_start^3*T_step)
  }
  x_ini_cov<-function(x){
    return(sigma^2*T_step)
  }
  x_mean<-function(x){
    return(x-x^3*T_step)
  }
  x_cov<-function(x){
    return(sigma^2*T_step)
  }
  potential_function<-function(x,y){
    return(0)
  }
  tk_list[[1]]<-list(x_ini_mean,x_ini_cov,potential_function)
  for (i in 2:(k-2)) {
    tk_list[[i]]<-list(x_mean,x_cov,potential_function)
  }
  potential_T<-function(x,y){
    return(dnorm(x_end,x_mean(x),sqrt(x_cov(x)),log = T))
  }
  tk_list[[k-1]]<-list(x_mean,x_cov,potential_T)
  return(tk_list)
}


cubic_bridge_LT_list_gen<-function(T_step,x_start,x_end,sigma,k){
  force(x_start);force(x_end)
  tk_list<-vector(mode = "list",k-1)
  T_step<-T_step/k
  xi_var<-(sigma^2/2)*(1-exp(-2*T_step))
  x_ini_mean<-function(x){

    return(exp(-T_step)*(x_start)/sqrt(exp(-2*T_step)+x_start^2*(1-exp(-2*T_step))))
  }
  x_ini_cov<-function(x){
    return(xi_var)
  }
  x_mean<-function(x){
    x_hat<-(x)/sqrt(exp(-2*T_step)+x^2*(1-exp(-2*T_step)))
    return(exp(-T_step)*x_hat)
  }
  x_cov<-function(x){
    return(xi_var)
  }
  potential_function<-function(x,y){
    return(0)
  }
  tk_list[[1]]<-list(x_ini_mean,x_ini_cov,potential_function)
  for (i in 2:(k-2)) {
    tk_list[[i]]<-list(x_mean,x_cov,potential_function)
  }
  potential_T<-function(x,y){
    return(dnorm(x_end,x_mean(x),sqrt(x_cov(x)),log = T))
  }
  tk_list[[k-1]]<-list(x_mean,x_cov,potential_T)
  return(tk_list)
}

f_inv<-function(X,t,e){
  a<-exp(-2*t/e)
  numtr<-a*X^2
  denom<-1-(1-a)*X^2
  if(sum(denom<0)!=0){
    browser()
  }else{

    return(sign(X)*sqrt(numtr/denom))}
}
cubic_bridge_S_list_gen<-function(T_step,x_start,x_end,sigma,k){


  T_step<-T_step/k
  z_start<-f_inv(x_start,T_step/2,1)
  z_end<-f_inv(x_end,T_step/2,1)
  is_both_num<-is.nan(c(z_start,z_end))
  while (any(is_both_num)) {
    z_start<-f_inv(x_start,T_step/2,1)
    z_end<-f_inv(x_end,T_step/2,1)
    is_both_num<-is.nan(c(z_start,z_end))
    if(any(is_both_num)){
      T_step<-T_step/2
      k<-2*k
    }
  }
  tk_list<-vector(mode = "list",k-1)
  xi_var<-(sigma^2/2)*(1-exp(-2*T_step))
  force(z_start);force(z_end)
  f<-expression(z_end/sqrt(exp(-T_step)+z_end^2*(1-exp(-T_step))))
  Df<-D(f,name="z_end")
  x_ini_mean<-function(x){
    return(exp(-T_step)*(z_start)/sqrt(exp(-2*T_step)+z_start^2*(1-exp(-2*T_step))))
  }
  x_ini_cov<-function(x){
    return(xi_var)
  }
  x_mean<-function(x){
    x_hat<-(x)/sqrt(exp(-2*T_step)+x^2*(1-exp(-2*T_step)))
    return(exp(-T_step)*x_hat)
  }
  x_cov<-function(x){
    return(xi_var)
  }
  potential_function<-function(x,y){
    return(0)
  }
  tk_list[[1]]<-list(x_ini_mean,x_ini_cov,potential_function)
  for (i in 2:(k-2)) {
    tk_list[[i]]<-list(x_mean,x_cov,potential_function)
  }
  potential_T<-function(x,y){
    return(dnorm(z_end,x_mean(x),sqrt(x_cov(x)),log = T)-log(abs(eval(expr = Df))))
  }
  tk_list[[k-1]]<-list(x_mean,x_cov,potential_T)
  return(tk_list)
}


###########generate a grid around sigma=20
seed<-as.integer(231012)
result<-vector(mode="list",length = 5)
par_grid<-seq(10,35,by=0.5)
result_vec_LT<-vector(mode = "numeric",length = length(par_grid))
result_vec_S<-vector(mode = "numeric",length = length(par_grid))

T_step<-1e-1
T_max<-100;
T_step_sim<-1e-4
sigma<-20
set.seed(231012)
data_obs<-cubic_sampler_strang_ini(T_step_sim,T_max,sigma,0)
data_obs_1000<-data_obs[[1]][1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
set.seed(NULL)
k_vec<-c(5,10,20,40,100)
cl<-makeCluster(12)
registerDoParallel(cl)

for (i in 1:5) {


k<-k_vec[i]

for (g in 1:length(par_grid)) {

  result_vec_LT[g]<-sum(unlist(foreach(t=1:1000) %dopar% {tk_now<-cubic_bridge_LT_list_gen(
    T_step,data_obs_1000[t],data_obs_1000[t+1],par_grid[g],k);cSMC(100,1:length(tk_now),tk_now)}))
  result_vec_S[g]<-sum(unlist(foreach(t=1:1000) %dopar% {tk_now<-cubic_bridge_LT_list_gen(
    T_step,data_obs_1000[t],data_obs_1000[t+1],par_grid[g],k);cSMC(100,1:length(tk_now),tk_now)}))
  print(c(result_vec_LT[g],result_vec_S[g],par_grid[g]))
}

result[[i]]<-cbind(result_vec_LT,result_vec_S)
print(i)

}

save(result,file = "ps_likelihood.rda")

sig_vec<-c(25,30,35,40)
explosion_result<-matrix(1:28,nrow = 7)
for (j in 1:4) {
  err_k<-vector()
  T_step<-1e-1
  T_max<-100;
  T_step_sim<-1e-4
  sigma<-sig_vec[j]
  data_obs<-cubic_sampler_strang_ini(T_step_sim,T_max,sigma,ini)
  data_obs_1000<-data_obs[[1]][1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
  for(k in 4:10){
    error_count <- 0
    abs_gt_20_count <- 0
    for (i in 1:1000) {
      tryCatch({
        result <- cSMC(10,1:(k-1),cubic_bridge_EM_list_gen(T_step,data_obs_1000[i],data_obs_1000[i+1],sigma,k))

        if (abs(result) > 10) {
          abs_gt_20_count <<- abs_gt_20_count + 1
        }
      }, error = function(e) {
        error_count <<- error_count + 1
      })
    }

    err_k<-c(err_k,error_count+abs_gt_20_count)
  }
  explosion_result[,j]<-err_k
}

explosion_list<-vector(mode = "list",length = 9)
explosion_list[[1]]<-explosion_result

for (m in 2:9) {
  sig_vec<-c(25,30,35,40)
  explosion_result<-matrix(1:28,nrow = 7)
  for (j in 1:4) {
    err_k<-vector()
    T_step<-1e-1
    T_max<-100;
    T_step_sim<-1e-4
    sigma<-sig_vec[j]
    data_obs<-cubic_sampler_strang_ini(T_step_sim,T_max,sigma,ini)
    data_obs_1000<-data_obs[[1]][1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
    for(k in 4:10){
      error_count <- 0
      abs_gt_20_count <- 0
      for (i in 1:1000) {
        tryCatch({
          result <- cSMC(10,1:(k-1),cubic_bridge_EM_list_gen(T_step,data_obs_1000[i],data_obs_1000[i+1],sigma,k))

          if (abs(result) > 10) {
            abs_gt_20_count <<- abs_gt_20_count + 1
          }
        }, error = function(e) {
          error_count <<- error_count + 1
        })
      }

      err_k<-c(err_k,error_count+abs_gt_20_count)
    }
    explosion_result[,j]<-err_k
  }
  explosion_list[[m]]<-explosion_result
}
