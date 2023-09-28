#########################################
# Simulation program for continuous outcomes, pcgs in the tutorial paper
#########################################

set.seed(96338, kind = "Mersenne-Twister", normal.kind = "Inversion")
marginal_model = 2 # incremental intervention effect model
Corr ="NE"
t=22 #period
s=6 # step
m=1 # cluster per step
n=s*m # total cluster
weights = c(0.1,0.2,0.5,0.1,0.1)
values <- c(2,3,4,5,6) # sample the cluster period size
N <- matrix(sample(values, size = 6*22, prob = weights,replace = TRUE),6,22)
b=7 # The last period before intervention in the 1st cluster
c=2 # 2 periods to train the staff between control and intervention
delta=3 # assumed intervention effect after 10 periods
beta = matrix(c(18.5, 0.02),2,1) # baseline mean +period effect
alpha = c(0.1,0.05)# ICC under Nested exchangeable correlation structure
max_intervention_periods= 10
Kt = matrix(c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
              0, 1, 1, 1, 1, 1 ,1, 0, 0 ,1, 1, 1, 1 ,1, 1, 1 ,1,1 ,0, 0, 0, 0,
              0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1 ,1, 0 ,0, 0,
              0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0 ,
              0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0,
              0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1),s,t,byrow=TRUE)


########################################################
# incidence_matrix 
# expand the incidence matrix Kt defined in cluster into subject-specify matrix
# Returns a subject-specify incidence matrix K

########################################################


incidence_matrix_i  =function(Kti , Nj){
  Kti_by_Nj = rep(Kti,times = Nj)
  matrix1 = diag(Kti_by_Nj)
  Ki = matrix1[rowSums(matrix1)==1,]
  return(Ki)
}


folder = "/Users/zying/Downloads/New Folder With Items/R codes to generate pcgs for CH study"
# Source program code
source(paste0(folder,"/Continuous_PCORI/TutorialPaper_CH_GEN_CONT.R"))

# Create X matrix
# Create treatment sequences
# remember to change treatment sequences under the specific designs
trtSeq<-matrix(0,s,t)
time_effect = cbind(matrix(1,t,1),c(1:(t))); # for continuous time effect
theta = rbind(beta,delta) # total parameters
s2 =2 # The variance of the continuous variable is 1 
# The following loop genearte X and Z, the design matrix in the regression model and correlation model
generateX = function(s,m,t,b,c,Kt,N,Corr){
  X  =NULL
  for(i in 1:s){
    Ni = N[i,]
    periods_before_int = b+c*(i-1)
    if (marginal_model==1){
      trtSeq[i,] = cbind(matrix(0, 1,periods_before_int),matrix(1, 1,t-periods_before_int))}else{
        trtSeq[i,] = cbind(matrix(0, 1,periods_before_int),matrix((1:(t-periods_before_int))/max_intervention_periods, 1,t-periods_before_int))}
    Xi = cbind(time_effect, matrix(trtSeq[i,],t,1)) # Xi, design matrix in the ith sequence
    Ki  = incidence_matrix_i(Kt[i,] , Ni)
    Xi_star = Ki%*%Xi[rep(1:nrow(Xi), times = Ni), ]
    m0=0
    cluster_sub_size = sum(N[Kt[i,]==1])# calculate the real cluster size

    repeat{
      m0 = m0+1
      X =rbind(X, Xi_star) # repeat the row of Xi by N
      if(m0==m){break}
    }
   
  }
  return(X)
}

X= generateX(s,m,t,b,c,Kt,N,Corr)

# Create id
cluster_sub_size =  rowSums(N*Kt)
id<-rep(1:n,cluster_sub_size)
clsize<-cluster_sub_size
start_time <- Sys.time()  
nerror = 0
y<-try(contGEN_PCORI_Tutorial(m,t,s,N,b,c,Corr,s2, delta,beta,alpha,marginal_model, max_intervention_periods, Kt),silent=TRUE)
pcgs_integer = round(y,digits=0)
output_data = data.frame(cbind(id,y,X,pcgs_integer))
colnames(output_data)=c("clusters","y", "int","periods", "intervention","pcgs_integer")
output_data$intervention_binary = as.numeric(output_data$intervention>0)
write.csv(output_data, "/Users/zying/Downloads/sample_data_PCGS_tutorial_paper_ICC01.csv")

