

log_sum_exp<-function(w){
  c<-max(w)
  return(c+log(sum(exp(w-c))))
}






#' Title
#'
#' @param w
#' @param is.resamp
#' @param if.print
#'
#' @return
#' @export
#'
#' @examples
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






#computes the density of a bunch of mvn variable with
#the same mean and diagonal cov matrix


#' @export
LGSSM<-function(A,T_max,y_cov){
  p<-dim(A)[1]
  X<-matrix(nrow = T_max+1,ncol = p)
  Y<-matrix(nrow = T_max+1,ncol = p)
  X[1,]<-mvnfast::rmvn(n=1,rep(0,p),diag(rep(1,p)))
  Y[1,]<-X[1,]+mvnfast::rmvn(n=1,rep(0,p),diag(y_cov,nrow = p))
  for (i in 1:T_max) {
    X[i+1,]<-A%*%X[i,]+as.numeric(mvnfast::rmvn(n=1,rep(0,p),diag(rep(1,p))))
    Y[i+1,]<-X[i+1,]+as.numeric(mvnfast::rmvn(n=1,rep(0,p),diag(y_cov,nrow = p)))
  }
  return(list(X=X,Y=Y))
}

#' @export
Kal_Simple<-function(A,y,y_cov){
  p<-dim(A)[1]
  T_max<-dim(y)[1]-1
  S<-B<-diag(rep(1,p))
  R<-diag(y_cov,nrow = p)
  m<-matrix(nrow = 1+T_max,ncol=p)
  Q<-vector(mode = "list",1+T_max)
  l<-mvnfast::dmvn(y[1,],rep(0,p),S+R,log = T)
  #initialize i.e. x_0|y_0
  Q_1<-solve(diag(rep(1,p))+solve(R))
  m_0<-Q_1%*%solve(R)%*%y[1,]
  m[1,]<-m_0
  Q[[1]]<-Q_1
  for (t in 1:T_max) {
    E_t<-A%*%Q[[t]]%*%t(A)+S
    Q_t<-E_t%*%(diag(rep(1,p))-solve(E_t+R)%*%E_t)
    m_t<-(diag(rep(1,p))-E_t%*%solve(E_t+R))%*%A%*%m[t,]+E_t%*%solve(E_t+R)%*%y[1+t,]
    Q[[1+t]]<-Q_t
    m[1+t,]<-m_t
    l_now<-mvnfast::dmvn(y[1+t,],A%*%m[t,],E_t+R,log = T)
    l<-l+l_now
  }
  return(l)
}

#' @export
grid_sim<-function(M,N,tk_list,y){
  y<-as.matrix(y)
  T_max<-dim(y)[1]-1
  #specify the grid
  grid_N<-seq(from=-M,to=M,length.out=N)
  tk_char<-c("x_mean","x_cov","y_mean","y_cov")
  tk_now<-tk_list[[1]]
  for (i in 1:length(tk_char)) {
    assign(tk_char[i],tk_now[[i]])
  }
  y_mean_vec<-unlist(lapply(grid_N, y_mean))
  y_sd_vec<-sqrt(unlist(lapply(grid_N, y_cov)))
  int_vec_y<-dnorm(y[1],y_mean_vec,y_sd_vec,log = T)
  x_mean_vec<-unlist(lapply(grid_N,x_mean))
  x_sd_vec<-sqrt(unlist(lapply(grid_N,x_cov)))
  int_vec_x<-dnorm(grid_N,x_mean_vec,x_sd_vec,log = T)
  #likelihood factor
  L_factor_now<-log(2*M/N)+log_sum_exp(int_vec_x+int_vec_y)
  L_vec<-L_factor_now
  #update
  p_filt_now<-int_vec_x+int_vec_y-L_factor_now
  for (t in 1:T_max) {
    tk_now<-tk_list[[t]]
    for (i in 1:length(tk_char)) {
      assign(tk_char[i],tk_now[[i]])
    }
    #prediction
    x_mean_vec<-unlist(lapply(grid_N,x_mean))
    x_sd_vec<-sqrt(unlist(lapply(grid_N,x_cov)))
    p_pred_now<-numeric(N)
    for (n in 1:N) {
      int_vec_x<-dnorm(grid_N[n],x_mean_vec,x_sd_vec,log = T)
      p_pred_now[n]<-log(2*M/N)+log_sum_exp(int_vec_x+p_filt_now)
    }
    #likelihood factor
    y_mean_vec<-unlist(lapply(grid_N, y_mean))
    y_sd_vec<-sqrt(unlist(lapply(grid_N, y_cov)))
    int_vec_y<-dnorm(y[1+t],y_mean_vec,y_sd_vec,log = T)
    L_factor_now<-log(2*M/N)+log_sum_exp(int_vec_y+p_pred_now)
    L_vec<-c(L_vec,L_factor_now)
    #update
    p_filt_now<-p_pred_now+int_vec_y-L_factor_now
  }
  return(sum(L_vec))
}

