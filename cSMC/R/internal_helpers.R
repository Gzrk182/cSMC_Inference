
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

