# complete case analysis of incomplete HHN

library(geeCRT)

##################################################################
# Screened for Smoking
##################################################################

HHN_smoking_screened_c <- read.xlsx("HHN_smoking_screened_c.xlsx", colNames=TRUE)

# cluster identifier (217 clusters in total)
id <- (HHN_smoking_screened_c$site_id)

# cluster period sizes
m <- HHN_smoking_screened_c$smoking_screened_denom
summary(m)
summary(as.numeric(tapply(m,id,sum)))

# responses
y <- HHN_smoking_screened_c$smoking_screened_num/HHN_smoking_screened_c$smoking_screened_denom
summary(y)

# period identifier
period <- as.numeric(factor(HHN_smoking_screened_c$quarter))

# cluster size (across all periods)
n <- as.numeric(tapply(period,id,length))
range(n)


#################################################
# Average intervention effects model
#################################################

# marginal mean design matrix
trt <- as.numeric(HHN_smoking_screened_c$phase>0)
X <- NULL
for(i in 1:nrow(HHN_smoking_screened_c)){
  time <- rep(0,max(period))
  time[period[i]] <- 1
  X <- rbind(X,c(time, trt[i]))
}

###### cluster-period mean analysis

### NE

# GEE/UEE
fit_aie_nex_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_aie_nex_uee
# GEE and UEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta  Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.4490162 0.09662820 0.10313292 0.10338117 0.10363002 0.10334532
# [2,]    1 0.4298590 0.09664277 0.09999372 0.10025345 0.10051397 0.10026163
# [3,]    2 0.4408298 0.09751382 0.10339021 0.10367520 0.10396110 0.10372579
# [4,]    3 0.4229450 0.10024561 0.10661288 0.10694210 0.10727269 0.10708633
# [5,]    4 0.3873187 0.10263020 0.10844576 0.10880398 0.10916401 0.10901449
# [6,]    5 0.3279662 0.10862785 0.11283582 0.11322668 0.11361964 0.11354766
# [7,]    6 0.3234389 0.10862500 0.10805948 0.10843616 0.10881494 0.10873448
# [8,]    7 0.2522880 0.10801497 0.10512505 0.10550326 0.10588360 0.10582370
# [9,]    8 0.2098761 0.10792918 0.10363695 0.10401650 0.10439821 0.10435456
# [10,]    9 0.1902485 0.10825931 0.10382042 0.10420244 0.10458664 0.10454455
# [11,]   10 0.1682365 0.10849768 0.10401178 0.10439630 0.10478299 0.10473535
# [12,]   11 0.2364589 0.05253848 0.07174799 0.07214152 0.07253733 0.07217055
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4698677 0.02421176 0.02427091 0.02433020 0.02427157
# [2,]     1 0.3914685 0.02692738 0.02699681 0.02706643 0.02699904

# GEE/MAEE
fit_aie_nex_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_aie_nex_maee
# GEE and MAEE for Cluster-Period Summaries
# Number of Iterations: 9
# Results for marginal mean parameters
# Beta  Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.4490266 0.09686769 0.10313279  0.1033810 0.10362985 0.10334513
# [2,]    1 0.4298471 0.09688236 0.09999393  0.1002537 0.10051416 0.10026170
# [3,]    2 0.4407922 0.09775643 0.10338909  0.1036741 0.10395994 0.10372446
# [4,]    3 0.4229005 0.10049916 0.10661339  0.1069426 0.10727318 0.10708657
# [5,]    4 0.3872855 0.10289322 0.10844529  0.1088035 0.10916352 0.10901365
# [6,]    5 0.3279071 0.10891321 0.11283446  0.1132253 0.11361824 0.11354581
# [7,]    6 0.3233742 0.10891031 0.10805921  0.1084359 0.10881463 0.10873375
# [8,]    7 0.2522215 0.10829897 0.10512491  0.1055031 0.10588342 0.10582308
# [9,]    8 0.2097717 0.10821269 0.10363575  0.1040153 0.10439696 0.10435286
# [10,]    9 0.1901851 0.10854495 0.10381947  0.1042015 0.10458564 0.10454308
# [11,]   10 0.1681664 0.10878419 0.10401101  0.1043955 0.10478216 0.10473408
# [12,]   11 0.2365330 0.05270116 0.07174880  0.0721423 0.07253809 0.07217129
# 
# Results for correlation parameters
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4722098 0.02433030 0.02438973 0.02444931 0.02439039
# [2,]     1 0.3933123 0.02706339 0.02713317 0.02720314 0.02713540


### ED

# GEE/UEE
fit_aie_ed_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_aie_ed_uee
# GEE and UEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.40847176 0.09823782 0.10088324 0.10112304 0.10136342 0.10118102
# [2,]    1 0.42468154 0.09806761 0.09860833 0.09884363 0.09907951 0.09892296
# [3,]    2 0.46905290 0.09861658 0.10151059 0.10175777 0.10200558 0.10184848
# [4,]    3 0.51527459 0.10088745 0.10171224 0.10195769 0.10220375 0.10197461
# [5,]    4 0.51222710 0.10253265 0.09955306 0.09979579 0.10003917 0.09983748
# [6,]    5 0.49727885 0.10635341 0.09811784 0.09836145 0.09860578 0.09842099
# [7,]    6 0.49016617 0.10646178 0.09660981 0.09684883 0.09708857 0.09691624
# [8,]    7 0.41107346 0.10579230 0.09211824 0.09234694 0.09257633 0.09240947
# [9,]    8 0.36909969 0.10570726 0.08942413 0.08964570 0.08986795 0.08970003
# [10,]    9 0.33394747 0.10592659 0.08906782 0.08928965 0.08951215 0.08934196
# [11,]   10 0.31324021 0.10668569 0.08896909 0.08919246 0.08941652 0.08924142
# [12,]   11 0.06467265 0.04050582 0.03074304 0.03088774 0.03103315 0.03084416
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4864894 0.024325095 0.02438288 0.024440802 0.024375197
# [2,]     1 0.9398190 0.009743582 0.00977124 0.009798978 0.009778093

# GEE/MAEE
fit_aie_ed_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_aie_ed_maee
# GEE and MAEE for Cluster-Period Summaries
# Number of Iterations: 9
# Results for marginal mean parameters
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.40851819 0.09846858 0.10088495 0.10112474 0.10136511 0.10118270
# [2,]    1 0.42470659 0.09829813 0.09861072 0.09884601 0.09908190 0.09892532
# [3,]    2 0.46903691 0.09884824 0.10151260 0.10175978 0.10200759 0.10185046
# [4,]    3 0.51527670 0.10112517 0.10171381 0.10195926 0.10220531 0.10197618
# [5,]    4 0.51223728 0.10277426 0.09955350 0.09979623 0.10003961 0.09983790
# [6,]    5 0.49728564 0.10660401 0.09811809 0.09836170 0.09860602 0.09842123
# [7,]    6 0.49017048 0.10671257 0.09661019 0.09684922 0.09708895 0.09691661
# [8,]    7 0.41107166 0.10604145 0.09211846 0.09234716 0.09257655 0.09240968
# [9,]    8 0.36908084 0.10595588 0.08942396 0.08964553 0.08986777 0.08969986
# [10,]    9 0.33393234 0.10617582 0.08906844 0.08929026 0.08951277 0.08934259
# [11,]   10 0.31322584 0.10693670 0.08896996 0.08919334 0.08941739 0.08924230
# [12,]   11 0.06468432 0.04060090 0.03074321 0.03088790 0.03103329 0.03084432
# 
# Results for correlation parameters
# Alpha  Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4887967 0.024442013 0.024500075 0.024558277 0.024492358
# [2,]     1 0.9398064 0.009747049 0.009774716 0.009802464 0.009781562

# formatting results
library(xtable)
beta_row <- function(fit){
  a <- round(fit$beta,3)
  b <- round(sqrt(diag(fit$BC1)[1:12]),3)
  paste0(a," (",b,")")
}

tab1 <- data.frame(c1=beta_row(fit_aie_nex_uee),c2=beta_row(fit_aie_nex_maee),c3=beta_row(fit_aie_ed_uee),c4=beta_row(fit_aie_ed_maee))
xtable(tab1)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllll}
# \hline
# & c1 & c2 & c3 & c4 \\ 
# \hline
# 1 & 0.449 (0.103) & 0.449 (0.103) & 0.408 (0.101) & 0.409 (0.101) \\ 
# 2 & 0.43 (0.1) & 0.43 (0.1) & 0.425 (0.099) & 0.425 (0.099) \\ 
# 3 & 0.441 (0.104) & 0.441 (0.104) & 0.469 (0.102) & 0.469 (0.102) \\ 
# 4 & 0.423 (0.107) & 0.423 (0.107) & 0.515 (0.102) & 0.515 (0.102) \\ 
# 5 & 0.387 (0.109) & 0.387 (0.109) & 0.512 (0.1) & 0.512 (0.1) \\ 
# 6 & 0.328 (0.113) & 0.328 (0.113) & 0.497 (0.098) & 0.497 (0.098) \\ 
# 7 & 0.323 (0.108) & 0.323 (0.108) & 0.49 (0.097) & 0.49 (0.097) \\ 
# 8 & 0.252 (0.106) & 0.252 (0.106) & 0.411 (0.092) & 0.411 (0.092) \\ 
# 9 & 0.21 (0.104) & 0.21 (0.104) & 0.369 (0.09) & 0.369 (0.09) \\ 
# 10 & 0.19 (0.104) & 0.19 (0.104) & 0.334 (0.089) & 0.334 (0.089) \\ 
# 11 & 0.168 (0.104) & 0.168 (0.104) & 0.313 (0.089) & 0.313 (0.089) \\ 
# 12 & 0.236 (0.072) & 0.237 (0.072) & 0.065 (0.031) & 0.065 (0.031) \\ 
# \hline
# \end{tabular}
# \end{table}


###### CIC
Ustar <- function(y,m,X,beta){
  u <- 1/(1+exp(c(-X%*%beta)))
  v <- u*(1-u)
  t(X)%*%(X*v*m)
}

sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_nex_uee$beta) %*% (fit_aie_nex_uee$BC1[1:12,1:12]) )) # 9639.852
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_nex_maee$beta) %*% (fit_aie_nex_maee$BC1[1:12,1:12]) )) # 9639.852
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_ed_uee$beta) %*% (fit_aie_ed_uee$BC1[1:12,1:12]) )) # 9133.806
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_ed_maee$beta) %*% (fit_aie_ed_maee$BC1[1:12,1:12]) )) # 9133.928



#################################################
# Extended incremental intervention effects model
#################################################

# marginal mean design matrix
cohort <- HHN_smoking_screened_c$cohort
trt <- matrix(c(rep(0, 1), 1/4, 2/4, 3/4, 1, rep(1, 6),
                rep(0, 2), 1/4, 2/4, 3/4, 1, rep(1, 5),
                rep(0, 3), 1/4, 2/4, 3/4, 1, rep(1, 4),
                rep(0, 3), 1/4, 2/4, 3/4, 1, rep(1, 4),
                rep(0, 4), 1/4, 2/4, 3/4, 1, rep(1, 3),
                rep(0, 5), 1/4, 2/4, 3/4, 1, rep(1, 2)), 
              nrow=6, ncol=11, byrow = TRUE)
X <- NULL
for(i in 1:nrow(HHN_smoking_screened_c)){
  time <- rep(0,max(period))
  time[period[i]] <- 1
  X <- rbind(X,c(time, trt[cohort[i], period[i]]))
}

###### cluster-period mean analysis

### NE

# GEE/UEE
fit_eiie_nex_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_eiie_nex_uee
# GEE and UEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.45085438 0.09675719 0.10324439  0.1034928  0.1037418  0.1034373
# [2,]    1 0.45134425 0.09661622 0.09981357  0.1000583  0.1003037  0.1000352
# [3,]    2 0.46545231 0.09715691 0.10205718  0.1023208  0.1025851  0.1023815
# [4,]    3 0.46353705 0.09881554 0.10332554  0.1036176  0.1039106  0.1038586
# [5,]    4 0.39404451 0.10306862 0.10906568  0.1094364  0.1098088  0.1099794
# [6,]    5 0.31336481 0.11121910 0.11952009  0.1200124  0.1205080  0.1208366
# [7,]    6 0.23474613 0.12034364 0.13220903  0.1328451  0.1334859  0.1338448
# [8,]    7 0.12161111 0.12604063 0.14017145  0.1409013  0.1416367  0.1420733
# [9,]    8 0.05188059 0.13035513 0.14597132  0.1467652  0.1475652  0.1480580
# [10,]    9 0.03165376 0.13081000 0.14671189  0.1475163  0.1483269  0.1488186
# [11,]   10 0.01029343 0.13082371 0.14625399  0.1470604  0.1478731  0.1483545
# [12,]   11 0.39391163 0.09145325 0.12585320  0.1266694  0.1274917  0.1273263
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4708834 0.02429371 0.02435287 0.02441217 0.02435387
# [2,]     1 0.3921114 0.02698992 0.02705932 0.02712890 0.02706192

# GEE/MAEE
fit_eiie_nex_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_eiie_nex_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.45086466 0.09700444 0.10324418  0.1034926  0.1037416  0.1034371
# [2,]    1 0.45133743 0.09686287 0.09981414  0.1000589  0.1003043  0.1000357
# [3,]    2 0.46542606 0.09740534 0.10205699  0.1023206  0.1025849  0.1023810
# [4,]    3 0.46351617 0.09907090 0.10332524  0.1036173  0.1039103  0.1038578
# [5,]    4 0.39402983 0.10334103 0.10906371  0.1094344  0.1098068  0.1099766
# [6,]    5 0.31332635 0.11152214 0.11951729  0.1200096  0.1205051  0.1208326
# [7,]    6 0.23468636 0.12067984 0.13220690  0.1328429  0.1334836  0.1338413
# [8,]    7 0.12154117 0.12639790 0.14016921  0.1408989  0.1416343  0.1420694
# [9,]    8 0.05176591 0.13072745 0.14596778  0.1467615  0.1475614  0.1480528
# [10,]    9 0.03158033 0.13118486 0.14670880  0.1475131  0.1483235  0.1488138
# [11,]   10 0.01021341 0.13119861 0.14625130  0.1470576  0.1478702  0.1483502
# [12,]   11 0.39399574 0.09174193 0.12584934  0.1266654  0.1274876  0.1273219
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4733041 0.02441392 0.02447337 0.02453296 0.02447438
# [2,]     1 0.3940211 0.02712933 0.02719908 0.02726902 0.02720168


### ED

# GEE/UEE
fit_eiie_ed_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_eiie_ed_uee
# GEE and UEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta  Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.4090376 0.09823890 0.10100314 0.10124259 0.10148262 0.10119984
# [2,]    1 0.4258826 0.09796937 0.09803415 0.09826593 0.09849826 0.09821662
# [3,]    2 0.4618862 0.09863894 0.10040397 0.10064792 0.10089249 0.10059669
# [4,]    3 0.4928908 0.10205320 0.10170672 0.10196491 0.10222382 0.10189088
# [5,]    4 0.4553882 0.10936367 0.10219507 0.10247972 0.10276534 0.10243264
# [6,]    5 0.4044489 0.12211127 0.10467155 0.10499357 0.10531698 0.10496701
# [7,]    6 0.3516454 0.13568544 0.10874911 0.10911679 0.10948634 0.10915787
# [8,]    7 0.2465991 0.14418493 0.10955317 0.10994766 0.11034435 0.11003481
# [9,]    8 0.1880969 0.15019931 0.11030585 0.11071793 0.11113240 0.11084331
# [10,]    9 0.1532741 0.15026845 0.10951133 0.10992167 0.11033443 0.11004731
# [11,]   10 0.1326388 0.15083982 0.10900851 0.10941831 0.10983051 0.10954249
# [12,]   11 0.2440114 0.11467006 0.07394022 0.07435207 0.07476651 0.07430162
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4864495 0.024401200 0.024459083 0.024517106 0.024452380
# [2,]     1 0.9398609 0.009737481 0.009765163 0.009792926 0.009770119

# GEE/MAEE
fit_eiie_ed_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_eiie_ed_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta  Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.4090841 0.09847822 0.10100484 0.10124429 0.10148431 0.10120153
# [2,]    1 0.4259077 0.09820818 0.09803681 0.09826859 0.09850092 0.09821927
# [3,]    2 0.4618698 0.09887936 0.10040641 0.10065036 0.10089492 0.10059911
# [4,]    3 0.4928958 0.10230360 0.10170802 0.10196621 0.10222511 0.10189215
# [5,]    4 0.4553994 0.10963412 0.10219508 0.10247972 0.10276533 0.10243261
# [6,]    5 0.4044552 0.12241622 0.10467194 0.10499395 0.10531735 0.10496735
# [7,]    6 0.3516445 0.13602677 0.10874957 0.10911723 0.10948677 0.10915826
# [8,]    7 0.2465902 0.14454885 0.10955417 0.10994865 0.11034532 0.11003572
# [9,]    8 0.1880687 0.15057913 0.11030662 0.11071869 0.11113316 0.11084400
# [10,]    9 0.1532502 0.15064852 0.10951343 0.10992377 0.11033651 0.11004934
# [11,]   10 0.1326151 0.15122132 0.10901099 0.10942078 0.10983298 0.10954492
# [12,]   11 0.2440307 0.11496696 0.07394126 0.07435310 0.07476753 0.07430262
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4888403 0.024520033 0.024578199 0.024636505 0.024571460
# [2,]     1 0.9398360 0.009741725 0.009769419 0.009797193 0.009774377

# formatting results
library(xtable)
beta_row <- function(fit){
  a <- round(fit$beta,3)
  b <- round(sqrt(diag(fit$BC1)[1:12]),3)
  paste0(a," (",b,")")
}

tab2 <- data.frame(c1=beta_row(fit_eiie_nex_uee),c2=beta_row(fit_eiie_nex_maee),c3=beta_row(fit_eiie_ed_uee),c4=beta_row(fit_eiie_ed_maee))
xtable(tab2)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllll}
# \hline
# & c1 & c2 & c3 & c4 \\ 
# \hline
# 1 & 0.451 (0.103) & 0.451 (0.103) & 0.409 (0.101) & 0.409 (0.101) \\ 
# 2 & 0.451 (0.1) & 0.451 (0.1) & 0.426 (0.098) & 0.426 (0.098) \\ 
# 3 & 0.465 (0.102) & 0.465 (0.102) & 0.462 (0.101) & 0.462 (0.101) \\ 
# 4 & 0.464 (0.104) & 0.464 (0.104) & 0.493 (0.102) & 0.493 (0.102) \\ 
# 5 & 0.394 (0.109) & 0.394 (0.109) & 0.455 (0.102) & 0.455 (0.102) \\ 
# 6 & 0.313 (0.12) & 0.313 (0.12) & 0.404 (0.105) & 0.404 (0.105) \\ 
# 7 & 0.235 (0.133) & 0.235 (0.133) & 0.352 (0.109) & 0.352 (0.109) \\ 
# 8 & 0.122 (0.141) & 0.122 (0.141) & 0.247 (0.11) & 0.247 (0.11) \\ 
# 9 & 0.052 (0.147) & 0.052 (0.147) & 0.188 (0.111) & 0.188 (0.111) \\ 
# 10 & 0.032 (0.148) & 0.032 (0.148) & 0.153 (0.11) & 0.153 (0.11) \\ 
# 11 & 0.01 (0.147) & 0.01 (0.147) & 0.133 (0.109) & 0.133 (0.109) \\ 
# 12 & 0.394 (0.127) & 0.394 (0.127) & 0.244 (0.074) & 0.244 (0.074) \\ 
# \hline
# \end{tabular}
# \end{table}


###### CIC
Ustar <- function(y,m,X,beta){
  u <- 1/(1+exp(c(-X%*%beta)))
  v <- u*(1-u)
  t(X)%*%(X*v*m)
}

sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_nex_uee$beta) %*% (fit_eiie_nex_uee$BC1[1:12,1:12]) )) # 9858.545
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_nex_maee$beta) %*% (fit_eiie_nex_maee$BC1[1:12,1:12]) )) # 9858.485
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_ed_uee$beta) %*% (fit_eiie_ed_uee$BC1[1:12,1:12]) )) # 9262.172
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_ed_maee$beta) %*% (fit_eiie_ed_maee$BC1[1:12,1:12]) )) # 9262.298


