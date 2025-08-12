library(cSMC)
library(mvtnorm)
library(mvnfast)
library(MASS)
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

set.seed(231012)
T_max <- 20;T_step <- 2e-5;length.out<-1000
FHN_S_obs_mcmc<-FHN_S_Gnr(T_max,T_step,0.1,1.5,0.8,0,0.3)
set.seed(NULL)
V_mcmc<-FHN_S_obs_mcmc$X$X1[1:(floor(T_max/T_step)+1)%%(round(T_max/T_step)/length.out)==1]
T_step_mcmc<-T_max/length.out
PMMH_S_h_bpf<-function(V,T_step,N_mcmc,N_smc,prop_step,pre_cond,theta_0){
  p<-length(theta_0)
  result<-matrix(nrow = N_mcmc+1, ncol=p)
  result[1,]<-theta_0
  MVN_array<-MASS::mvrnorm(N_mcmc,rep(0,p),prop_step*pre_cond)
  u_array<-runif(N_mcmc)
  eAdelta<-matrix_eAt(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])
  Cdelta<-matrix_Ct(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])
  FHN_S_tk_list<-FHN_S_tk_list_Gnr(V[-1],T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])

  a<-exp(-T_step/theta_0[1])
  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  z_t<-f_inv(V[2],T_step/2,theta_0[1])
  pf_result<-cSMC::bpf(N_smc,V[-(1:2)],FHN_S_tk_list,0.5)
  L_0<-cSMC::extract_likelihood(pf_result[[2]],pf_result[[3]])+
    mvnfast::dmvn(f_inv(V[2],T_step/2,theta_0[1]),eAdelta[1,1]*f_t_1(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0)+
                    eAdelta[1,2]*f_t_2(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0),
                  Cdelta[1,1],log = T)-log(abs(eval(Df)))
  theta_prior<-theta_0
  for (n in 1:N_mcmc) {
    theta_1<-exp(log(theta_0)+MVN_array[n,])
    eAdelta<-matrix_eAt(T_step,theta_1[1],theta_1[2],theta_1[3],0,theta_1[4])
    Cdelta<-matrix_Ct(T_step,theta_1[1],theta_1[2],theta_1[3],0,theta_1[4])
    FHN_S_tk_list<-FHN_S_tk_list_Gnr(V[-1],T_step,theta_1[1],theta_1[2],theta_1[3],0,theta_1[4])
    a<-exp(-T_step/theta_1[1])
    z_t<-f_inv(V[2],T_step/2,theta_1[1])
    pf_result<-cSMC::bpf(N_smc,V[-(1:2)],FHN_S_tk_list,0.5)
    L_1<-cSMC::extract_likelihood(pf_result[[2]],pf_result[[3]])+
      mvnfast::dmvn(f_inv(V[2],T_step/2,theta_1[1]),eAdelta[1,1]*f_t_1(T_step,theta_1[1],theta_1[2],theta_1[3],0,theta_1[4],0)+
                      eAdelta[1,2]*f_t_2(T_step,theta_1[1],theta_1[2],theta_1[3],0,theta_1[4],0),
                    Cdelta[1,1],log = T)-log(abs(eval(Df)))
    accept_prob<-L_1+mvnfast::dmvn(log(theta_1),log(theta_prior),diag(rep(1,p)),log = T)-
      L_0-mvnfast::dmvn(log(theta_0),log(theta_prior),diag(rep(1,p)),log = T)
    if(log(u_array[n])<=accept_prob){
      result[(n+1),]<-theta_1
      theta_0<-theta_1
      L_0<-L_1
    }else{
      result[(n+1),]<-theta_0
    }
    print(c(result[(n),],n))
  }

  return(result)
}


PMMH_result_bpf<-PMMH_S_h_bpf(V_mcmc,T_step_mcmc,1e5,125,prop_step = 2,pre_cond = pred_cov,c(0.1,1.5,0.8,0.3))
save(PMMH_result_bpf,file = "bpf_hypo_1e5.rda")
