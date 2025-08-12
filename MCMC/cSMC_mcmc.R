
library(MASS)
library(mvtnorm)
library(cSMC)
log_sum_exp<-function(w){
  c<-max(w)
  return(c+log(sum(exp(w-c))))
}

FHN_EM_Gnr<-function(T_max=100,T_step=1e-3,e=1,g=5,b=0.1,s_1=0.1,s_2=0.2,X_0=c(0,0)){
  T_total<-floor(T_max/T_step)
  X_matrix<-matrix(nrow = 2,ncol = 1+T_total)
  X_matrix[,1]<-X_0
  for (i in 1:T_total) {
    X_tm1<-X_matrix[,i]
    V_now<-(X_tm1[1]-X_tm1[1]^3-X_tm1[2])/e
    U_now<-g*X_tm1[1]-X_tm1[2]+b
    X_now<-X_tm1+T_step*c(V_now,U_now)+c(rnorm(1,0,sqrt(s_1^2*T_step)),rnorm(1,0,sqrt(s_2^2*T_step)))
    X_matrix[,i+1]<-X_now
    print(i)
  }
  return(data.frame(t(X_matrix)))
}



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
  T_total<-floor(T_max/T_step)
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
#validation


FHN_S_tk_list_Gnr<-function(V,T_step,e,g,b,s_1,s_2){
  V<-f_inv(V,T_step/2,e)
  T_max<-length(V)-1
  tk_list<-vector(mode = "list",T_max)
  a<-exp(-T_step/e)
  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  mu_a_0<-eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,0)+
    eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,0)
  mu_b_0<-eAdelta[2,1]*f_t_1(T_step,e,g,b,s_1,s_2,0)+
    eAdelta[2,2]*f_t_2(T_step,e,g,b,s_1,s_2,0)
  x_0_mean<-function(x){
    return(mu_b_0+Cdelta[2,1]*(V[1]-mu_a_0)/Cdelta[1,1])
  }
  x_0_cov<-function(x){
    return(Cdelta[2,2]-Cdelta[2,1]^2/Cdelta[1,1])
  }
  g_1_potential<-function(x,y){
    z_t<-V[2]
    mu_a_1<-eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)+
      eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,V[1])
    return(mvnfast::dmvn(V[2],mu_a_1,Cdelta[1,1],log = T)-log(abs(eval(expr = Df))))
  }
  tk_list[[1]]<-list(x_0_mean,x_0_cov,g_1_potential)
  sigma_f<-Cdelta[1,1]
  sigma_p<-Cdelta[2,2]
  x_mean<-function(V_tm1,V_t){
    force(V_t);force(V_tm1)
    return(
      function(x){
        mu_b<-eAdelta[2,1]*f_t_1(T_step,e,g,b,s_1,s_2,V_tm1)+
          eAdelta[2,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)
        mu_a<-eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)+
          eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,V_tm1)
        return(mu_b+Cdelta[2,1]*(V_t-mu_a)/Cdelta[1,1])
      }
    )
  }
  x_cov<-function(x){
    return(Cdelta[2,2]-Cdelta[2,1]^2/Cdelta[1,1])
  }
  potential_function<-function(V_tm1,V_t){
    force(V_tm1);force(V_t)
    return(
      function(x,y){
        z_t<-V_t
        mu<-eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)+
          eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,V_tm1)
        return(mvnfast::dmvn(V_t,mu,Cdelta[1,1],log = T)-log(abs(eval(expr = Df))))
      }
    )
  }
  for (t in 2:T_max) {
    tk_list[[t]]<-list(x_mean(V[t-1],V[t]),x_cov,potential_function(V[t],V[t+1]))
  }
  return(tk_list)
}

set.seed(231012)
T_max <- 20;T_step <- 2e-5;length.out<-1000
FHN_S_obs_mcmc<-FHN_S_Gnr(T_max,T_step,0.1,1.5,0.8,0,0.3)
set.seed(NULL)
V_mcmc<-FHN_S_obs_mcmc$X$X1[1:(floor(T_max/T_step)+1)%%(round(T_max/T_step)/length.out)==1]
T_step_mcmc<-T_max/length.out
PMMH_S_h<-function(V,T_step,N_mcmc,N_smc,prop_step,pre_cond,theta_0){
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
  L_0<-cSMC::cSMC(N_smc,V[-(1:2)],FHN_S_tk_list,k=1)+
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
    L_1<-cSMC::cSMC(N_smc,V[-(1:2)],FHN_S_tk_list,k=1)+
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


PMMH_result<-PMMH_S_h(V_mcmc,T_step_mcmc,100000,10,prop_step = 2,pre_cond = cov(log(exploratory_run)),c(0.1,1.5,0.8,0.3))
save(PMMH_result,file = "hypo_1e5.rda")
#####################
theta_0<-colMeans(PMMH_result[-(1:2000),])
set.seed(231012)
T_max <- 20;T_step <- 2e-5;length.out<-1000
FHN_S_obs_mcmc<-FHN_S_Gnr(T_max,T_step,0.1,1.5,0.8,0,0.3)
set.seed(NULL)
V<-FHN_S_obs_mcmc$X$X1[1:(floor(T_max/T_step)+1)%%(round(T_max/T_step)/length.out)==1]
T_step<-T_max/length.out

eAdelta<-matrix_eAt(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])
Cdelta<-matrix_Ct(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])
FHN_S_tk_list<-FHN_S_tk_list_Gnr(V[-1],T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4])
a<-exp(-T_step/theta_0[1])
f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
Df<-D(f,name="z_t")
z_t<-f_inv(V[2],T_step/2,theta_0[1])

bpfVcsmc<-matrix(1:2000,nrow = 1000)
for (m in 1:1000) {
  L_cSMC<-cSMC::cSMC(10,V[-(1:2)],FHN_S_tk_list,k=1)+
    mvnfast::dmvn(f_inv(V[2],T_step/2,theta_0[1]),eAdelta[1,1]*f_t_1(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0)+
                    eAdelta[1,2]*f_t_2(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0),
                  Cdelta[1,1],log = T)-log(abs(eval(Df)))

  result_BPF<-cSMC::bpf(125,V[-(1:2)],FHN_S_tk_list,0.5)
  L_BPF<-extract_likelihood(result_BPF[[2]],result_BPF[[3]])+
    mvnfast::dmvn(f_inv(V[2],T_step/2,theta_0[1]),eAdelta[1,1]*f_t_1(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0)+
                    eAdelta[1,2]*f_t_2(T_step,theta_0[1],theta_0[2],theta_0[3],0,theta_0[4],0),
                  Cdelta[1,1],log = T)-log(abs(eval(Df)))
  bpfVcsmc[m,]<-c(L_cSMC,L_BPF);print(m)
}

save(bpfVcsmc,file = "bpfVcsmc.rda")



