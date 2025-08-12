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
  if(names(attributes(potential_function))=="p_now"){
    p<-attributes(potential_function)$p_now
  }
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
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }
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
  if(names(attributes(potential_function))=="p_now"){
    p<-attributes(potential_function)$p_now
  }
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
