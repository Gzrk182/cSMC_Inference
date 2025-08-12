

#' Controlled Sequantial Monte Carlo
#'
#' @param N Number of particles
#' @param y Measurement process as T by N matrix
#' @param tk_list A length T list. The t-th element is a list containing the mean and covariance function of the Gaussian transition
#' kernel, and the potential function, at time t. See example. Same as for the `bpf`.
#' @param k Integer that controls stopping. Stops at the k-th iteration where no resampling happened.
#' @param alpha Controls adaptive resampling
#'
#' @return The normalizing constant estimate in log-scale.
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
#'cSMC(500,y,tk_list_LG_5,3,0.5)
cSMC<-function(N,y,tk_list,k=3,alpha=0.5,if_print=F){
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

    psi_spec<-recu_f_approx(pf_result[[1]],y,tk_list)
    pf_result<-twist_bpf(N,y,tk_list,alpha,psi_spec)
    L_vec<-c(L_vec,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    if(sum(pf_result[[3]])==0){
      l<-l+1
      keep_vec<-c(keep_vec,extract_likelihood(pf_result[[2]],pf_result[[3]]))
    }
    if(I>10){
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
  return(mean(keep_vec))
}
