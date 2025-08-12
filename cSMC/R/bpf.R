
#' The Bootstrap Particle Filter
#'
#' @param N Number of particles
#' @param y Measurement process as T by N matrix
#' @param tk_list A length T list. The t-th element is a list containing the mean and covariance function of the Gaussian transition
#' kernel, and the potential function, at time t. See example.
#' @param alpha A number controls adaptive resampling.
#'
#' @return A list containing the particles (N by p matrix), weights (N by p matrix) and a vector indicating at which times
#' the algorithm has resampled.
#' @export
#'
#' @examples
#' A<-diag(c(0.9,0.7,0.6,0.3,0))
#'A[lower.tri(A)]<-c(0.3,0.1,0.4,0.1,0.2,0.1,0.2,0.1,0.5,0.2)
#'LG_5_x_mean<-function(x){
#'  return(A%*%x)
#'}
#'LG_5_x_cov<-function(x){
#'return(diag(rep(1,5)))
#'}
#'LG_5_potential_log<-function(x,y){
#'  return(mvnfast::dmvn(y,x,diag(rep(0.25,5)),log = T))
#'}
#'tk_list_LG_5<-vector(mode = "list",101)
#'for (i in 1:101) {
#'  tk_list_LG_5[[i]]<-list(LG_5_x_mean,LG_5_x_cov,LG_5_potential_log)
#'}
#'T_data <- 100
#'data_obs<-LGSSM(A,T_data,0.25)
#'y<-data_obs$Y
#'bpf(500,y,tk_list_LG_5,0.5)

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
  if(names(attributes(potential_function))=="p_now"){
    p<-attributes(potential_function)$p_now
  }
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
  set.seed(1)
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
    if(names(attributes(potential_function))=="p_now"){
      p<-attributes(potential_function)$p_now
    }
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
