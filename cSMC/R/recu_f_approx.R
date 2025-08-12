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
prod_cross(1:3)

matrix_place<-function(x){
  n<-floor(sqrt(length(x)*2))
  result<-matrix(rep(0,n^2),nrow = n,ncol = n)
  result[lower.tri(result,diag = T)]<-x
  result<-result+t(result)
  return(result/2)
}

matrix_place(prod_cross(1:4))

qr.solve(matrix_place(prod_cross(1:4)))

matrix_make<-function(x){
  N<-dim(x)[1]
  p<-dim(x)[2]
  quadr_block<-matrix(nrow = N,ncol = p*(p+1)/2)
  for (n in 1:N) {
    quadr_block[n,]<-prod_cross(x[n,])
  }
  return(cbind(quadr_block,x,rep(1,n)))
}

#' Title
#'
#' @param x
#' @param y
#' @param tk_list
#'
#' @return
#' @export
recu_f_approx<-function(x,y,tk_list){
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
  beta_t<-solve.qr(qr(Z),psi_t)
  beta_matrix<-beta_t
  beta_t<-beta_t[-length(beta_t)]
  Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]))

  beta_t<-beta_t[-(1:(p*(p+1)/2))]
  mu_t<-Sigma_t%*%beta_t
  if(det(as.matrix(Sigma_t))<0){Sigma_t<-diag(1e3,p)}
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
    beta_t<-solve.qr(qr(Z),psi_t)
    beta_matrix<-rbind(beta_matrix,beta_t)
    beta_t<-beta_t[-length(beta_t)]
    Sigma_t<--0.5*qr.solve(matrix_place(beta_t[1:(p*(p+1)/2)]))

    beta_t<-beta_t[-(1:(p*(p+1)/2))]
    mu_t<-Sigma_t%*%beta_t
    if(det(as.matrix(Sigma_t))<0){Sigma_t<-diag(1e3,p)}
    opt_par_t<-list(mu_t,Sigma_t)
    psi_spec[[t]]<-opt_par_t
  }

  return(psi_spec)

  #return(psi_T_max)
}

