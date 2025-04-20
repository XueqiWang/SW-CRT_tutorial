# analysis of HIV Testing

library(haven)
library(dplyr)
library(data.table)
library(geeCRT)
library(EigenR)
library(MASS)

# Source the binomial_maee function 
# which uses a different inverse function: Eigen_pinverse rather than ginv to aviod the warning from SVD function in the ginv function
source("binomial_maee.R")

##################################################################
# HIV Testing
##################################################################

hivtesting <- read_sas("hivtesting.sas7bdat")

hivtest <- hivtesting %>%
  arrange(clusternum, ID, time)

# cluster identifier
id <- as.numeric(hivtest$clusternum)

# responses
y <- as.numeric(hivtest$hivt)

### Constant intervention effects model

# marginal mean design matrix
period <- as.matrix(hivtest[, 10:13])
stra <- as.numeric(hivtest$Shandong)
trt <- as.numeric(hivtest$intervention)
X <- cbind(period, stra, trt)


#################################################
# Constrained Block Exchangeable Correlation
#################################################

# Create Z matrix
dt <- as.data.table(hivtest)
setorder(dt, clusternum, ID, time)

dt[, personcnt := 0L]
dt[, obsnum := 1:.N, by = clusternum]
dt[, personcnt := rleid(ID), by = clusternum]

dt1 <- copy(dt)
dt2 <- copy(dt)

setnames(dt1, c("ID", "time", "hivt", "obsnum", "personcnt"),
         c("ID1", "time1", "hivt1", "obsnum1", "personcnt1"))

setnames(dt2, c("ID", "time", "hivt", "obsnum", "personcnt"),
         c("ID2", "time2", "hivt2", "obsnum2", "personcnt2"))

pairs <- dt1[dt2, on = .(clusternum), allow.cartesian = TRUE, nomatch = 0]
pairs <- pairs[obsnum1 < obsnum2]

pairs[, z1 := (ID1 != ID2 & time1 == time2)]
pairs[, z2 := (ID1 != ID2 & time1 != time2)]
pairs[, z3 := (ID1 == ID2)]

summary_table <- pairs[, .(
  sumz1 = sum(z1),
  sumz2 = sum(z2),
  sumz3 = sum(z3)
)]
print(summary_table)

setorder(pairs, clusternum, obsnum1, obsnum2)
Z <- as.matrix(pairs[, .(z1, z3)]*1)

# GEE/MAEE with CBE
fit_cie_cbe_maee <- binomial_maee(y=y, X=X, id=id, Z=Z, link="logit", maxiter=100, epsilon=1e-4, alpadj=TRUE, makevone=FALSE)
fit_cie_cbe_maee
# GEE for correlated data 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta    Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -1.45832201 0.1527441 0.09923136 0.11589011  0.1357349 0.11233249
# [2,]    1 -0.96068817 0.1540112 0.11783113 0.13228059  0.1496527 0.14075314
# [3,]    2 -0.88368979 0.1739720 0.08857089 0.09855129  0.1130534 0.09240868
# [4,]    3 -0.67658650 0.1982700 0.13809556 0.15545421  0.1779356 0.15954028
# [5,]    4 -0.00169826 0.1271550 0.12172866 0.14461308  0.1723164 0.14603365
# [6,]    5  0.27323544 0.1526455 0.11372203 0.13488412  0.1631521 0.13334493
# 
# Results for correlation parameters 
# Alpha   Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.01507586 0.005689983 0.006090455 0.006521222 0.006097931
# [2,]     1 0.21726462 0.021007606 0.023070277 0.025359235 0.022489565


#################################################
# Extended 3-Dependence Correlation
#################################################

# Create Z matrix
dt <- as.data.table(hivtest)
setorder(dt, clusternum, ID, time)

dt[, personcnt := 0L]
dt[, obsnum := 1:.N, by = clusternum]
dt[, personcnt := rleid(ID), by = clusternum]

dt1 <- copy(dt)
dt2 <- copy(dt)

setnames(dt1, c("ID", "time", "hivt", "obsnum", "personcnt"),
         c("ID1", "time1", "hivt1", "obsnum1", "personcnt1"))

setnames(dt2, c("ID", "time", "hivt", "obsnum", "personcnt"),
         c("ID2", "time2", "hivt2", "obsnum2", "personcnt2"))

pairs <- dt1[dt2, on = .(clusternum), allow.cartesian = TRUE, nomatch = 0]
pairs <- pairs[obsnum1 < obsnum2]

pairs[, z1 := (ID1 != ID2 & time1 == time2)]
pairs[, z2 := (ID1 == ID2 & (time2 - time1 == 1))]
pairs[, z3 := (ID1 == ID2 & (time2 - time1 == 2))]
pairs[, z4 := (ID1 == ID2 & (time2 - time1 == 3))]

summary_table <- pairs[, .(
  sumz1 = sum(z1),
  sumz2 = sum(z2),
  sumz3 = sum(z3),
  sumz4 = sum(z4)
)]
print(summary_table)

setorder(pairs, clusternum, obsnum1, obsnum2)
Z <- as.matrix(pairs[, .(z1, z2, z3, z4)]*1)

# GEE/MAEE with E3DE
fit_cie_e3d_maee <- binomial_maee(y=y, X=X, id=id, Z=Z, link="logit", maxiter=100, epsilon=1e-4, alpadj=TRUE, makevone=FALSE)
fit_cie_e3d_maee
# GEE for correlated data 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta      Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 -1.4582318133 0.1527954 0.09853341 0.11512802  0.1349268 0.11173264
# [2,]    1 -0.9581537175 0.1540961 0.11882130 0.13340506  0.1509202 0.14236631
# [3,]    2 -0.8853537558 0.1742509 0.08618084 0.09607324  0.1106203 0.09032969
# [4,]    3 -0.6808808319 0.1986500 0.13702239 0.15446541  0.1771844 0.15932916
# [5,]    4 -0.0006029704 0.1272102 0.12120821 0.14414218  0.1719591 0.14570226
# [6,]    5  0.2735284252 0.1531630 0.11588034 0.13756376  0.1664252 0.13642955
# 
# Results for correlation parameters 
# Alpha   Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.01507489 0.005703503 0.006104103 0.006534946 0.006111479
# [2,]     1 0.25887872 0.030962559 0.033930093 0.037199138 0.033191034
# [3,]     2 0.18385469 0.017456037 0.018854982 0.020406698 0.018696151
# [4,]     3 0.15404411 0.015696802 0.017041337 0.018538306 0.016805009


