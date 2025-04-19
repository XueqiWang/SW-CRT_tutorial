# analysis of CH

library(geeCRT)

##################################################################
# CH
##################################################################

CH_c <- read.csv("../../Data/sample_data_PCGS_tutorial_paper_ICC01.csv", header = TRUE)

# cluster identifier
id <- (CH_c$clusters)

# responses
y <- CH_c$pcgs_integer

# period identifier
period <- as.numeric(factor(CH_c$periods))

# treatment effect (from SAS code)
incre_trt <- ex_incre_trt <- NULL
for (i in 1:nrow(CH_c)){
  if (CH_c[i, 2] == 1){
    if (CH_c[i, 5] <= 5){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-7)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-7)/5,1))
    }
  }
  if (CH_c[i, 2] == 2){
    if (CH_c[i, 5] <= 7){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-9)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-9)/5,1))
    }
  }
  if (CH_c[i, 2] == 3){
    if (CH_c[i, 5] <= 9){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-11)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-11)/5,1))
    }
  }
  if (CH_c[i, 2] == 4){
    if (CH_c[i, 5] <= 11){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-13)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-13)/5,1))
    }
  }
  if (CH_c[i, 2] == 5){
    if (CH_c[i, 5] <= 13){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-15)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-15)/5,1))
    }
  }
  if (CH_c[i, 2] == 6){
    if (CH_c[i, 5] <= 15){
      incre_trt = c(incre_trt, CH_c[i, 8])
      ex_incre_trt = c(ex_incre_trt, CH_c[i, 8])
    } else{
      incre_trt = c(incre_trt, (CH_c[i, 5]-17)/10)
      ex_incre_trt = c(ex_incre_trt, min((CH_c[i, 5]-17)/5,1))
    }
  }
}

CH_c$incre_trt = incre_trt
CH_c$ex_incre_trt = ex_incre_trt

# Z matrix (from R vignettes for geeCRT package)
m <- as.matrix(table(id, period))
n <- dim(m)[1]
t <- dim(m)[2]

createzCrossSec <- function(m) {
  Z <- NULL
  n <- dim(m)[1]
  
  for (i in 1:n) {
    alpha_0 <- 1
    alpha_1 <- 2
    n_i <- c(m[i, ])
    n_length <- length(n_i)
    POS <- matrix(alpha_1, sum(n_i), sum(n_i))
    loc1 <- 0
    loc2 <- 0
    
    for (s in 1:n_length) {
      n_t <- n_i[s]
      if (n_t == 0) next
      loc1 <- loc2 + 1
      loc2 <- loc1 + n_t - 1
      
      for (k in loc1:loc2) {
        for (j in loc1:loc2) {
          if (k != j) {
            POS[k, j] <- alpha_0
          } else {
            POS[k, j] <- 0
          }
        }
      }
    }
    
    zrow <- diag(2)
    z_c <- NULL
    
    for (j in 1:(sum(n_i) - 1)) {
      for (k in (j + 1):sum(n_i)) {
        z_c <- rbind(z_c, zrow[POS[j, k], ])
      }
    }
    
    Z <- rbind(Z, z_c)
  }
  
  return(Z)
}
Z <- createzCrossSec(m)



#################################################
# Constant intervention effects model
#################################################

# marginal mean design matrix
trt <- CH_c$intervention_binary
X <- NULL
for(i in 1:nrow(CH_c)){
  # add one intercept and the continuous periods effect based on former codes
  X <- rbind(X,c(1, CH_c$periods[i]-1, trt[i]))
}

# GEE/MAEE with NE
fit_cie_nex_maee <- geemaee(y=y, X=X, id=id, Z=Z, family="continuous", epsilon=1e-8, alpadj=TRUE, makevone = FALSE)
fit_cie_nex_maee
# GEE for correlated data 
# Number of Iterations: 63 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 17.6773659 0.35121338 0.26599454 0.29793924 0.33470660 0.34440531
# [2,]    1  0.2028370 0.03579409 0.03546168 0.03982535 0.04494164 0.04683231
# [3,]    2 -0.5448703 0.38025187 0.32486723 0.36707560 0.41702513 0.42724114
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.2202014 0.07160415 0.07918255 0.08788288 0.07821985
# [2,]     1 0.1445320 0.05806730 0.06392278 0.07059154 0.06263240



#################################################
# Incremental
#################################################

# marginal mean design matrix
trt <- CH_c$incre_trt
X <- NULL
for(i in 1:nrow(CH_c)){
  # add one intercept and the continuous periods effect based on former codes
  X <- rbind(X,c(1, CH_c$periods[i]-1, trt[i]))
}

# GEE/MAEE with NE
fit_ie_nex_maee <- geemaee(y=y, X=X, id=id, Z=Z, family="continuous", epsilon=1e-8, alpadj=TRUE, makevone = FALSE)
fit_ie_nex_maee
# GEE for correlated data 
# Number of Iterations: 45 
# Results for marginal mean parameters 
# Beta    Estimate  MB-stderr  BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 18.59559468 0.26663460 0.113762261 0.12159928 0.13309462 0.15036717
# [2,]    1  0.02791844 0.02881419 0.008922223 0.01135618 0.01519675 0.02378498
# [3,]    2  2.93257733 0.55561560 0.273283137 0.32422436 0.39204455 0.50988923
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.10229568 0.03615120 0.04152101 0.04856669 0.04643938
# [2,]     1 0.05235939 0.03930308 0.04585742 0.05360908 0.05141365



#################################################
# Extended incremental
#################################################

# marginal mean design matrix
trt <- CH_c$ex_incre_trt
X <- NULL
for(i in 1:nrow(CH_c)){
  # add one intercept and the continuous periods effect based on former codes
  X <- rbind(X,c(1, CH_c$periods[i]-1, trt[i]))
}

# GEE/MAEE with NE
fit_eie_nex_maee <- geemaee(y=y, X=X, id=id, Z=Z, family="continuous", epsilon=1e-8, alpadj=TRUE, makevone = FALSE)
fit_eie_nex_maee
# GEE for correlated data 
# Number of Iterations: 37 
# Results for marginal mean parameters 
# Beta    Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 18.54566582 0.27747849 0.10440933 0.11029040 0.12090856 0.16577497
# [2,]    1  0.03110682 0.03272875 0.01412225 0.01687668 0.02130993 0.02443233
# [3,]    2  1.82453071 0.41762586 0.27437463 0.32422067 0.38896731 0.40870512
# 
# Results for correlation parameters 
# Alpha   Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.11992166 0.04243059 0.04892489 0.05681714 0.04723524
# [2,]     1 0.04809048 0.02122223 0.02280733 0.02467764 0.02213543


