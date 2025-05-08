library(PFSelectcopy)

T_0 = abs(c( 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         -1.72156244e-02,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         -4.50479857e+04,  7.99035330e-03,  0.00000000e+00,  0.00000000e+00))
T_k = abs(c( 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00))

T_0 <- abs(results$coefs[3:152])
T_k <- abs(results$coefs[303:452])

result=PFSelectcopy:::MK.statistic(T_0, T_k, method='median')

kappa <- result[,1]                                                             
tau <- result[,2]
b<-order(tau,decreasing=T)
c_0<-kappa[b]==0
M = 1
Rej.Bound=10000
#calculate ratios for top Rej.Bound tau values
ratio<-c();temp_0<-0
for(i in 1:length(b)){
  #if(i==1){temp_0=c_0[i]}
  temp_0<-temp_0+c_0[i]
  temp_1<-i-temp_0
  temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
  ratio<-c(ratio,temp_ratio)
  if(i>Rej.Bound){break}
}
#calculate q values for top Rej.Bound values
q<-rep(1,length(tau));
if(length(which(tau[b]>0))!=0){
  index_bound<-max(which(tau[b]>0))
  for(i in 1:length(b)){
    temp.index<-i:min(length(b),Rej.Bound,index_bound)
    if(length(temp.index)==0){next}
    q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  q[q>1]<-1
}
undebug(PFSelectcopy:::MK.q.byStat)
q_2 <- PFSelectcopy:::MK.q.byStat(kappa, tau, 1, 10000)
