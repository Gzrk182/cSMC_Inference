#library(R.matlab)
#read 1609
#dta<-readMat("~/ybhwtngzmm-1/Cutaneous Stimulation/1609.mat")[[1]][,2]
#subsample so the correct time step is 2e-2
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




FHN_LT_tk_list_Gnr<-function(V,T_step,e,g,b,s_1,s_2,ini=c(0,0)){
  T_max<-length(V)-1
  tk_list<-vector(mode = "list",T_max)

  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  mu_a_0<-eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,ini[1])+
    eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,ini[2])
  mu_b_0<-eAdelta[2,1]*f_t_1(T_step,e,g,b,s_1,s_2,ini[1])+
    eAdelta[2,2]*f_t_2(T_step,e,g,b,s_1,s_2,ini[2])

  x_0_mean<-function(x){
    return(mu_b_0+Cdelta[2,1]*(V[1]-mu_a_0)/Cdelta[1,1])
  }
  x_0_cov<-function(x){
    return(Cdelta[2,2]-Cdelta[2,1]^2/Cdelta[1,1])
  }
  g_1_potential<-function(x,y){
    mu_a_1<-eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)+
      eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,V[1])
    return(mvnfast::dmvn(V[2],mu_a_1,Cdelta[1,1],log = T))
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
        mu<-eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,x)+
          eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,V_tm1)
        return(mvnfast::dmvn(V_t,mu,Cdelta[1,1],log = T))
      }
    )
  }
  for (t in 2:T_max) {
    tk_list[[t]]<-list(x_mean(V[t-1],V[t]),x_cov,potential_function(V[t],V[t+1]))
  }
  return(tk_list)
}


obj_LT_h<-function(T_step,par,V){
  par<-exp(par)
  T_max<-2000
  T_step<-2e-2
  T_step_sim<-2e-2
  FHN_S_obs<-FHN_S_Gnr(T_max,T_step_sim,par[1], par[2], par[3], 0,par[4])
  V_20_h<-FHN_S_obs$X$X1[1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
  v_mean<-mean(FHN_S_obs$X$X1)
  u_mean<-mean(FHN_S_obs$X$X2)
  print(v_mean)
  V<-c(V+v_mean)
  eAdelta<-matrix_eAt(T_step,par[1],par[2],par[3],0,par[4])
  Cdelta<-matrix_Ct(T_step,par[1],par[2],par[3],0,par[4])
  FHN_LT_tk_list<-FHN_LT_tk_list_Gnr(V[-1],T_step,par[1],par[2],par[3],0,par[4],c(v_mean,u_mean))
  result<-cSMC::cSMC(10,V[-(1:2)],FHN_LT_tk_list,k=1)+
    mvnfast::dmvn(V[2],eAdelta[1,1]*f_t_1(T_step,par[1],par[2],par[3],0,par[4],v_mean)+
                    eAdelta[1,2]*f_t_2(T_step,par[1],par[2],par[3],0,par[4],u_mean),
                  Cdelta[1,1],log = T)
  #print(c(par,result))
  return(result)
}




FHN_S_tk_list_Gnr<-function(V,T_step,e,g,b,s_1,s_2,ini=c(0,0)){
  V<-f_inv(V,T_step/2,e)
  T_max<-length(V)-1
  tk_list<-vector(mode = "list",T_max)
  a<-exp(-T_step/e)
  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  eAdelta<-matrix_eAt(T_step,e,g,b,s_1,s_2)
  Cdelta<-matrix_Ct(T_step,e,g,b,s_1,s_2)
  mu_a_0<-eAdelta[1,1]*f_t_1(T_step,e,g,b,s_1,s_2,ini[1])+
    eAdelta[1,2]*f_t_2(T_step,e,g,b,s_1,s_2,ini[2])
  mu_b_0<-eAdelta[2,1]*f_t_1(T_step,e,g,b,s_1,s_2,ini[1])+
    eAdelta[2,2]*f_t_2(T_step,e,g,b,s_1,s_2,ini[2])
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

obj_S_h<-function(T_step,par,V){
  par<-exp(par)
  T_max<-2000
  T_step<-2e-2
  T_step_sim<-2e-2
  FHN_S_obs<-FHN_S_Gnr(T_max,T_step_sim,par[1], par[2], par[3], 0,par[4])
  V_20_h<-FHN_S_obs$X$X1[1:(1+ceiling(T_max/T_step_sim))%%round(T_step/T_step_sim)==1]
  v_mean<-mean(FHN_S_obs$X$X1)
  u_mean<-mean(FHN_S_obs$X$X2)
  print(v_mean)
  V<-c(V+v_mean)
  eAdelta<-matrix_eAt(T_step,par[1],par[2],par[3],0,par[4])
  Cdelta<-matrix_Ct(T_step,par[1],par[2],par[3],0,par[4])
  FHN_S_tk_list<-FHN_S_tk_list_Gnr(V[-1],T_step,par[1],par[2],par[3],0,par[4],c(v_mean,u_mean))

  f<-expression(z_t/sqrt(a+z_t^2*(1-a)))
  Df<-D(f,name="z_t")
  z_t<-f_inv(V[2],T_step/2,par[1])
  result<-cSMC::cSMC(10,V[-(1:2)],FHN_S_tk_list,k=1)+
    dmvn(f_inv(V[2],T_step/2,par[1]),eAdelta[1,1]*f_t_1(T_step,par[1],par[2],par[3],0,par[4],v_mean)+
           eAdelta[1,2]*f_t_2(T_step,par[1],par[2],par[3],0,par[4],u_mean),
         Cdelta[1,1],log = T)-log(abs(eval(Df)))
  #print(c(par,result))
  return(result)
}


a_vec<-1*(51:1e4)^(-0.602)
c_vec<-5e-2*(1:1e4)^(-0.101)
d_vec<-1e-1*(1:1e4)^(-0.101)

par_0<-log(c(0.2, 50,50,1))

K<-length(a_vec)
par_vec<-matrix(1:(5*(1+K)),ncol=5)
scale_par<-c(1,1,1,1)
ini_par<-par_0/scale_par
par_vec[1,]<-c(ini_par*scale_par,0)
T_step<-2e-2
h_k_old<-0
i<-1
while (i<5000) {
  I_vec<-vector()
  theta_old<-par_vec[i,1:4]
  #create initialization psi_spec for n_bridge=4


  Delta_now<-sample(c(1,-1),4,replace = T,prob = rep(0.5,2))

  y_p_f<-obj_S_h(2e-2,(theta_old/scale_par+c_vec[i]*Delta_now)*scale_par,dta)

  y_m_f<-obj_S_h(2e-2,(theta_old/scale_par-c_vec[i]*Delta_now)*scale_par,dta)
  delta_now<-sample(c(1,-1),4,replace = T,prob = rep(0.5,2))
  y_p_f_d<-obj_S_h(2e-2,(theta_old/scale_par+d_vec[i]*delta_now+c_vec[i]*Delta_now)*scale_par,dta)
  y_m_f_d<-obj_S_h(2e-2,(theta_old/scale_par+d_vec[i]*delta_now-c_vec[i]*Delta_now)*scale_par,dta)
  y_p<-y_p_f
  y_m<-y_m_f
  y_p_d<-y_p_f_d
  y_m_d<-y_m_f_d
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
  l_new<-obj_S_h(2e-2,par_new,dta)
  #if(l_new>((y_p+y_m)/2-0.1)){
  par_vec[i+1,1:4]<-par_new
  par_vec[i+1,5]<-l_new

  #  }else{
  #   par_vec[i+1,1:4]<-par_vec[i,1:4]
  #   par_vec[i+1,5]<-par_vec[i,5]
  #   print("blocked")
  #  }
  print(c(exp(par_vec[i+1,1:4]),l_new,i))
  h_k_old<-h_k_now
  i<-i+1
}

