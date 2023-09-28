#############################################################
# Generate correlated Gaussian outcomes
# for stepped wedge CRTs in PCORI proposal
# By Ying ZHang
# Jan 2021
# Based on Fan Li 2018

# n: Number of clusters
# m: Cluster size
# t: Number of periods
# delta: Effect size in odds ratio (OR)
# s2: Dispersion parameter or total variance (assume to be 1
#     and suppressed in the input)
# beta: Vector of period effects in OR
# alpha: Vector of correlations 
#        alpha_0=within period correlation
#        alpha_1=inter-period correlation
#        alpha_2=individual auto-correlation
#############################################################

contGEN_PCORI_Tutorial<-function(m,t,s,N,b,c,Corr,s2,delta,beta,alpha,marginal_model, max_intervention_periods, Kt){
  require(mvtnorm)
  require(magic)
  ########################################################
  # Create correlation matrix.
  ########################################################
  
   nxch<-function(alpha,N,t){
    sumN = sum(N)
    blockR = matrix(1, N[1],N[1])
    for (j in 2:t){
      subblock = matrix(1, N[j],N[j])
      blockR = adiag(blockR,subblock)
    }
    R = alpha[2]*matrix(1,sumN,sumN) + (alpha[1]-alpha[2])*blockR + (1-alpha[1])*diag(1,sumN,sumN)
    return(R)
  }
  
  #########################################################
  #Generate incidence matrix by missing matrix by each cluster
  #########################################################
  incidence_matrix_i  =function(Kti , Nj){
    Kti_by_Nj = rep(Kti,times = Nj)
    matrix1 = diag(Kti_by_Nj)
    Ki = matrix1[rowSums(matrix1)==1,]
    return(Ki)
  }
  # Create treatment sequences
  trtSeq<-matrix(0,s,t)
  #max_intervention_periods = t- b
  for (seq in 1:s){
    periods_before_int = b+c*(seq-1)
    if (marginal_model==1){#average intervention effect model
      trtSeq[seq,] = cbind(matrix(0, 1,periods_before_int),matrix(1, 1,t-periods_before_int))}
    else{#incremental intervention effect model
      trtSeq[seq,] = cbind(matrix(0, 1,periods_before_int),matrix((1:(t-periods_before_int))/max_intervention_periods, 1,t-periods_before_int))}
    #print(trtSeq[seq,])
  }
  #print (trtSeq)
  
  ########################################################
  # returns variance matrix of Gaussian variables with dispersion
  # parameter s2 and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,s2){
    return(s2*r)
  }
  
  # Simulate correlated continuous outcomes
  y<-NULL
  Zi =NULL # true covariates matrix in the 
  time_effect = cbind(matrix(1,t,1),c(1:(t))); # for continuous time effect
  theta = rbind(beta,delta) # total parameters

  for(i in 1:s){
    Ni=N[i,]
    if (Corr=="NE"){ri<-nxch(alpha,Ni,t)}else{print("The correlation structure is not available in the function")}
    v<-cor2var(ri,s2)   # v <- cov matrix
    Xi = cbind(time_effect, matrix(trtSeq[i,],t,1)) # Xi, design matrix in the ith sequence
    Zi = Xi[rep(1:nrow(Xi), times = Ni), ] # repeat the row of Xi by N
    Ki  = incidence_matrix_i(Kt[i,] , Ni)
    Zi_star  = Ki%*%Zi
    print(Zi_star)
    v_star = Ki%*%v%*%t(Ki)
    #print(r_star)
    u_star <-c(Zi_star%*%theta) # mu for the ith sequence
    m0=0
    repeat{
      m0 = m0+1
      yi = t(rmvnorm(1,u_star,v_star))
      y<-rbind(y,yi) # simulate data matrix
      if(m0==m){break}
    } 
  }
  
  # Return simulated data matrix
  return(y)
}