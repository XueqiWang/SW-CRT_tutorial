# complete case analysis of incomplete HHN

library(openxlsx)
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
# Constant intervention effects model
#################################################

# marginal mean design matrix
trt <- as.numeric(HHN_smoking_screened_c$phase>0)
stra <- as.numeric(HHN_smoking_screened_c$cohort<4)
X <- NULL
for(i in 1:nrow(HHN_smoking_screened_c)){
  time <- rep(0,max(period))
  time[period[i]] <- 1
  X <- rbind(X,c(time, trt[i], stra[i]))
}

###### cluster-period mean analysis

### NE

# GEE/UEE
fit_aie_nex_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_aie_nex_uee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.44321450 0.12220868 0.12525543 0.12571794 0.12618225 0.12553624
# [2,]    1 0.42415453 0.12168313 0.12288734 0.12336152 0.12383761 0.12323497
# [3,]    2 0.43517322 0.12193542 0.12609149 0.12658693 0.12708434 0.12649573
# [4,]    3 0.41733394 0.12367733 0.12939332 0.12992977 0.13046853 0.12992474
# [5,]    4 0.38170453 0.12565394 0.13221029 0.13277264 0.13333761 0.13283887
# [6,]    5 0.32231997 0.13054154 0.13727736 0.13786500 0.13845542 0.13805087
# [7,]    6 0.31779431 0.13052513 0.13389817 0.13447741 0.13505943 0.13463189
# [8,]    7 0.24664913 0.13003134 0.13112695 0.13170382 0.13228347 0.13188659
# [9,]    8 0.20424268 0.12989825 0.12987947 0.13045459 0.13103249 0.13066316
# [10,]    9 0.18463219 0.13010440 0.12968836 0.13026204 0.13083851 0.13048224
# [11,]   10 0.16260145 0.13037420 0.12992515 0.13049883 0.13107529 0.13072222
# [12,]   11 0.23633480 0.05261927 0.07163796 0.07203121 0.07242676 0.07200682
# [13,]   12 0.01382499 0.17988313 0.17520062 0.17605881 0.17692145 0.17619680
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4699155 0.02423229 0.02428825 0.02434434 0.02429200
# [2,]     1 0.3914478 0.02699209 0.02705732 0.02712272 0.02706459

# GEE/MAEE
fit_aie_nex_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_aie_nex_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.44330669 0.12274633 0.12525595 0.12571844 0.12618272 0.12553651
# [2,]    1 0.42421691 0.12221857 0.12288713 0.12336129 0.12383735 0.12323469
# [3,]    2 0.43520125 0.12247117 0.12608824 0.12658363 0.12708100 0.12649250
# [4,]    3 0.41734738 0.12421771 0.12939021 0.12992660 0.13046531 0.12992196
# [5,]    4 0.38173617 0.12619924 0.13220732 0.13276961 0.13333452 0.13283637
# [6,]    5 0.32232718 0.13109813 0.13727290 0.13786048 0.13845084 0.13804711
# [7,]    6 0.31779435 0.13108145 0.13389475 0.13447394 0.13505589 0.13462914
# [8,]    7 0.24664374 0.13058502 0.13112319 0.13170000 0.13227959 0.13188356
# [9,]    8 0.20418251 0.13045007 0.12987447 0.13044952 0.13102736 0.13065897
# [10,]    9 0.18462657 0.13065710 0.12968406 0.13025768 0.13083408 0.13047874
# [11,]   10 0.16259564 0.13092706 0.12992149 0.13049510 0.13107149 0.13071937
# [12,]   11 0.23642409 0.05279352 0.07163739 0.07203059 0.07242609 0.07200616
# [13,]   12 0.01368887 0.18069840 0.17519028 0.17604840 0.17691095 0.17618644
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4740437 0.02447765 0.02453420 0.02459088 0.02453797
# [2,]     1 0.3950418 0.02726293 0.02732885 0.02739493 0.02733613


### ED

# GEE/UEE
fit_aie_ed_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_aie_ed_uee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.37407020 0.1220949 0.12105176 0.12150721 0.12196444 0.12146456
# [2,]    1 0.39062801 0.1216538 0.11987515 0.12033123 0.12078910 0.12033448
# [3,]    2 0.43527888 0.1217948 0.12219568 0.12265917 0.12312447 0.12266357
# [4,]    3 0.48176718 0.1233309 0.12227974 0.12274182 0.12320570 0.12263875
# [5,]    4 0.47868602 0.1246749 0.12151882 0.12198271 0.12244844 0.12191729
# [6,]    5 0.46366942 0.1277899 0.12144941 0.12191729 0.12238706 0.12187706
# [7,]    6 0.45655196 0.1277949 0.12076914 0.12123370 0.12170014 0.12120740
# [8,]    7 0.37752913 0.1272274 0.11728553 0.11774149 0.11819932 0.11772007
# [9,]    8 0.33569016 0.1270412 0.11538114 0.11583071 0.11628214 0.11580610
# [10,]    9 0.30080067 0.1269880 0.11484200 0.11528859 0.11573704 0.11526746
# [11,]   10 0.28016571 0.1275279 0.11499763 0.11544406 0.11589234 0.11542661
# [12,]   11 0.06448893 0.0406822 0.03073299 0.03087822 0.03102416 0.03083382
# [13,]   12 0.08213392 0.1738169 0.17226160 0.17311891 0.17398072 0.17335189
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr BC2-stderr  BC3-stderr
# [1,]     0 0.4869837 0.024316880 0.024372336 0.02442792 0.024364633
# [2,]     1 0.9394228 0.009826846 0.009854319 0.00988187 0.009862648

# GEE/MAEE
fit_aie_ed_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_aie_ed_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta   Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.37421297 0.12260921 0.12105458 0.12151001 0.12196722 0.12146781
# [2,]    1 0.39072955 0.12216756 0.11987877 0.12033483 0.12079269 0.12033872
# [3,]    2 0.43531032 0.12230913 0.12219712 0.12266059 0.12312586 0.12266564
# [4,]    3 0.48181231 0.12384855 0.12227914 0.12274118 0.12320502 0.12263838
# [5,]    4 0.47874432 0.12519410 0.12151859 0.12198245 0.12244814 0.12191741
# [6,]    5 0.46373055 0.12831207 0.12144928 0.12191714 0.12238688 0.12187739
# [7,]    6 0.45661625 0.12831678 0.12077002 0.12123456 0.12170098 0.12120876
# [8,]    7 0.37758567 0.12774649 0.11728619 0.11774213 0.11819994 0.11772120
# [9,]    8 0.33571817 0.12755875 0.11538185 0.11583141 0.11628282 0.11580725
# [10,]    9 0.30082679 0.12750468 0.11484402 0.11529060 0.11573904 0.11526992
# [11,]   10 0.28020442 0.12804498 0.11500128 0.11544770 0.11589597 0.11543069
# [12,]   11 0.06450493 0.04078284 0.03073301 0.03087822 0.03102415 0.03083384
# [13,]   12 0.08206108 0.17459470 0.17225878 0.17311603 0.17397778 0.17334943
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr BC2-stderr  BC3-stderr
# [1,]     0 0.4910444 0.024553389 0.024609388 0.02466552 0.024601595
# [2,]     1 0.9396006 0.009801319 0.009828731 0.00985622 0.009837066

# formatting results
library(xtable)

beta_est_row <- function(fit){
  a <- round(fit$beta,3)
}
beta_se_row <- function(fit, BC){
  if (BC==0){
    b <- round(sqrt(diag(fit$BC0)[1:13]),3)
  } else if (BC==1){
    b <- round(sqrt(diag(fit$BC1)[1:13]),3)
  }
}

alpha_est_row <- function(fit){
  a <- round(fit$alpha,3)
}
alpha_se_row <- function(fit, BC){
  if (BC==0){
    b <- round(sqrt(diag(fit$BC0)[14:15]),3)
  } else if (BC==2){
    b <- round(sqrt(diag(fit$BC2)[14:15]),3)
  }
}

tab1_beta <- data.frame(c1=beta_est_row(fit_aie_nex_uee), c1s=beta_se_row(fit_aie_nex_uee,0), 
                        c2=beta_est_row(fit_aie_nex_maee), c2s=beta_se_row(fit_aie_nex_maee,1),
                        c3=beta_est_row(fit_aie_ed_uee), c3s=beta_se_row(fit_aie_ed_uee,0),
                        c4=beta_est_row(fit_aie_ed_maee), c4s=beta_se_row(fit_aie_ed_maee,1))
xtable(tab1_beta, digits=3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & c1 & c1s & c2 & c2s & c3 & c3s & c4 & c4s \\ 
# \hline
# 1 & 0.443 & 0.125 & 0.443 & 0.126 & 0.374 & 0.121 & 0.374 & 0.122 \\ 
# 2 & 0.424 & 0.123 & 0.424 & 0.123 & 0.391 & 0.120 & 0.391 & 0.120 \\ 
# 3 & 0.435 & 0.126 & 0.435 & 0.127 & 0.435 & 0.122 & 0.435 & 0.123 \\ 
# 4 & 0.417 & 0.129 & 0.417 & 0.130 & 0.482 & 0.122 & 0.482 & 0.123 \\ 
# 5 & 0.382 & 0.132 & 0.382 & 0.133 & 0.479 & 0.122 & 0.479 & 0.122 \\ 
# 6 & 0.322 & 0.137 & 0.322 & 0.138 & 0.464 & 0.121 & 0.464 & 0.122 \\ 
# 7 & 0.318 & 0.134 & 0.318 & 0.134 & 0.457 & 0.121 & 0.457 & 0.121 \\ 
# 8 & 0.247 & 0.131 & 0.247 & 0.132 & 0.378 & 0.117 & 0.378 & 0.118 \\ 
# 9 & 0.204 & 0.130 & 0.204 & 0.130 & 0.336 & 0.115 & 0.336 & 0.116 \\ 
# 10 & 0.185 & 0.130 & 0.185 & 0.130 & 0.301 & 0.115 & 0.301 & 0.115 \\ 
# 11 & 0.163 & 0.130 & 0.163 & 0.130 & 0.280 & 0.115 & 0.280 & 0.115 \\ 
# 12 & 0.236 & 0.072 & 0.236 & 0.072 & 0.064 & 0.031 & 0.065 & 0.031 \\ 
# 13 & 0.014 & 0.175 & 0.014 & 0.176 & 0.082 & 0.172 & 0.082 & 0.173 \\ 
# \hline
# \end{tabular}
# \end{table}

tab1_alpha <- data.frame(c1=alpha_est_row(fit_aie_nex_uee), c1s=alpha_se_row(fit_aie_nex_uee,0), 
                         c2=alpha_est_row(fit_aie_nex_maee), c2s=alpha_se_row(fit_aie_nex_maee,2),
                         c3=alpha_est_row(fit_aie_ed_uee), c3s=alpha_se_row(fit_aie_ed_uee,0),
                         c4=alpha_est_row(fit_aie_ed_maee), c4s=alpha_se_row(fit_aie_ed_maee,2))
xtable(tab1_alpha, digits=3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & c1 & c1s & c2 & c2s & c3 & c3s & c4 & c4s \\ 
# \hline
# 1 & 0.470 & 0.024 & 0.474 & 0.025 & 0.487 & 0.024 & 0.491 & 0.025 \\ 
# 2 & 0.391 & 0.027 & 0.395 & 0.027 & 0.939 & 0.010 & 0.940 & 0.010 \\ 
# \hline
# \end{tabular}
# \end{table}


###### CIC
Ustar <- function(y,m,X,beta){
  u <- 1/(1+exp(c(-X%*%beta)))
  v <- u*(1-u)
  t(X)%*%(X*v*m)
}

sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_nex_uee$beta) %*% (fit_aie_nex_uee$BC1[1:13,1:13]) )) # 16955.37
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_nex_maee$beta) %*% (fit_aie_nex_maee$BC1[1:13,1:13]) )) # 16954.5
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_ed_uee$beta) %*% (fit_aie_ed_uee$BC1[1:13,1:13]) )) # 16194.54
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_aie_ed_maee$beta) %*% (fit_aie_ed_maee$BC1[1:13,1:13]) )) # 16194.49



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
  X <- rbind(X,c(time, trt[cohort[i], period[i]], stra[i]))
}

###### cluster-period mean analysis

### NE

# GEE/UEE
fit_eiie_nex_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_eiie_nex_uee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta     Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0  0.452210798 0.12253349  0.1257022  0.1261662  0.1266321  0.1259932
# [2,]    1  0.452690034 0.12214589  0.1231042  0.1235708  0.1240393  0.1234384
# [3,]    2  0.466783107 0.12208146  0.1256226  0.1261057  0.1265907  0.1260573
# [4,]    3  0.464845932 0.12260816  0.1271149  0.1276238  0.1281348  0.1277507
# [5,]    4  0.395334382 0.12536325  0.1330968  0.1336677  0.1342412  0.1341113
# [6,]    5  0.314649573 0.13162381  0.1435202  0.1441907  0.1448647  0.1449224
# [7,]    6  0.236028402 0.13910982  0.1555760  0.1563662  0.1571613  0.1572849
# [8,]    7  0.122894306 0.14399721  0.1627631  0.1636290  0.1645005  0.1647342
# [9,]    8  0.053164849 0.14765952  0.1680291  0.1689475  0.1698719  0.1701829
# [10,]    9  0.032933346 0.14799384  0.1684153  0.1693407  0.1702721  0.1705906
# [11,]   10  0.011577155 0.14807443  0.1680860  0.1690117  0.1699436  0.1702555
# [12,]   11  0.393973682 0.09158902  0.1252758  0.1260905  0.1269113  0.1265831
# [13,]   12 -0.003274473 0.18024588  0.1742061  0.1750611  0.1759204  0.1752575
# 
# Results for correlation parameters 
# Alpha Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.470856 0.02423431 0.02428994 0.02434570 0.02429338
# [2,]     1 0.392101 0.02697596 0.02704083 0.02710585 0.02704788

# GEE/MAEE
fit_eiie_nex_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="nest_exch", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_eiie_nex_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta     Estimate  MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0  0.452283051 0.12307686  0.1257030  0.1261670  0.1266329  0.1259939
# [2,]    1  0.452738129 0.12268823  0.1231051  0.1235718  0.1240402  0.1234393
# [3,]    2  0.466803634 0.12262384  0.1256211  0.1261041  0.1265891  0.1260559
# [4,]    3  0.464864222 0.12315277  0.1271117  0.1276205  0.1281315  0.1277481
# [5,]    4  0.395359896 0.12591436  0.1330931  0.1336639  0.1342374  0.1341087
# [6,]    5  0.314647603 0.13219038  0.1435148  0.1441851  0.1448590  0.1449185
# [7,]    6  0.235996634 0.13969536  0.1555702  0.1563604  0.1571553  0.1572809
# [8,]    7  0.122845383 0.14459459  0.1627562  0.1636220  0.1644933  0.1647293
# [9,]    8  0.053051991 0.14826559  0.1680202  0.1689384  0.1698627  0.1701762
# [10,]    9  0.032875000 0.14860137  0.1684074  0.1693325  0.1702638  0.1705848
# [11,]   10  0.011518461 0.14868174  0.1680787  0.1690043  0.1699360  0.1702503
# [12,]   11  0.394095741 0.09190085  0.1252645  0.1260791  0.1268997  0.1265715
# [13,]   12 -0.003361002 0.18106810  0.1742001  0.1750550  0.1759143  0.1752517
# 
# Results for correlation parameters 
# Alpha  Estimate BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]     0 0.4750329 0.02448051 0.02453673 0.02459308 0.02454020
# [2,]     1 0.3957289 0.02725013 0.02731568 0.02738140 0.02732277


### ED

# GEE/UEE
fit_eiie_ed_uee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=FALSE)
fit_eiie_ed_uee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 8 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.37770219 0.1222181 0.12130178 0.12175752 0.12221505 0.12155007
# [2,]    1 0.39485693 0.1216980 0.11961504 0.12006873 0.12052421 0.11987806
# [3,]    2 0.43123197 0.1218071 0.12172048 0.12218167 0.12264466 0.12198355
# [4,]    3 0.46266016 0.1240194 0.12295846 0.12343131 0.12390603 0.12320060
# [5,]    4 0.42532229 0.1297517 0.12490414 0.12540125 0.12590042 0.12521808
# [6,]    5 0.37437314 0.1404729 0.12871832 0.12924856 0.12978112 0.12909702
# [7,]    6 0.32151860 0.1523974 0.13352178 0.13408948 0.13465986 0.13401764
# [8,]    7 0.21646147 0.1601402 0.13498735 0.13557441 0.13616437 0.13556301
# [9,]    8 0.15801125 0.1655678 0.13609969 0.13669868 0.13730070 0.13673436
# [10,]    9 0.12341369 0.1654603 0.13521935 0.13581448 0.13641265 0.13585613
# [11,]   10 0.10283254 0.1659145 0.13499734 0.13559115 0.13618800 0.13563937
# [12,]   11 0.24354137 0.1153808 0.07397152 0.07438448 0.07480003 0.07431964
# [13,]   12 0.07515126 0.1739439 0.17184566 0.17270114 0.17356111 0.17289770
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4874790 0.024277022 0.024331696 0.024386500 0.024321842
# [2,]     1 0.9394711 0.009813721 0.009841213 0.009868783 0.009847223

# GEE/MAEE
fit_eiie_ed_maee <- cpgeeSWD(y=y, X=X, id=id, m=m, corstr="exp_decay", family="binomial", epsilon=1e-8, alpadj=TRUE)
fit_eiie_ed_maee
# GEE and MAEE for Cluster-Period Summaries 
# Number of Iterations: 9 
# Results for marginal mean parameters 
# Beta   Estimate MB-stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,]    0 0.37783079 0.1227397  0.1213045 0.12176020 0.12221770 0.12155240
# [2,]    1 0.39494451 0.1222193  0.1196185 0.12007218 0.12052764 0.11988124
# [3,]    2 0.43124934 0.1223292  0.1217219 0.12218308 0.12264604 0.12198473
# [4,]    3 0.46269362 0.1245467  0.1229568 0.12342959 0.12390427 0.12319857
# [5,]    4 0.42536486 0.1302874  0.1249021 0.12539917 0.12589829 0.12521588
# [6,]    5 0.37441469 0.1410256  0.1287160 0.12924620 0.12977870 0.12909464
# [7,]    6 0.32155554 0.1529709  0.1335192 0.13408688 0.13465719 0.13401522
# [8,]    7 0.21648719 0.1607267  0.1349846 0.13557158 0.13616147 0.13556049
# [9,]    8 0.15800523 0.1661641  0.1360963 0.13669527 0.13729724 0.13673134
# [10,]    9 0.12340629 0.1660560  0.1352174 0.13581249 0.13641061 0.13585455
# [11,]   10 0.10283701 0.1665109  0.1349966 0.13559032 0.13618711 0.13563895
# [12,]   11 0.24357513 0.1156982  0.0739683 0.07438121 0.07479671 0.07431644
# [13,]   12 0.07511541 0.1747292  0.1718452 0.17270063 0.17356055 0.17289737
# 
# Results for correlation parameters 
# Alpha  Estimate  BC0-stderr  BC1-stderr  BC2-stderr  BC3-stderr
# [1,]     0 0.4916078 0.024514168 0.024569378 0.024624719 0.024559428
# [2,]     1 0.9396378 0.009789165 0.009816598 0.009844109 0.009822631

# formatting results
library(xtable)

beta_est_row <- function(fit){
  a <- round(fit$beta,3)
}
beta_se_row <- function(fit, BC){
  if (BC==0){
    b <- round(sqrt(diag(fit$BC0)[1:13]),3)
  } else if (BC==1){
    b <- round(sqrt(diag(fit$BC1)[1:13]),3)
  }
}
alpha_est_row <- function(fit){
  a <- round(fit$alpha,3)
}
alpha_se_row <- function(fit, BC){
  if (BC==0){
    b <- round(sqrt(diag(fit$BC0)[14:15]),3)
  } else if (BC==2){
    b <- round(sqrt(diag(fit$BC2)[14:15]),3)
  }
}

tab2_beta <- data.frame(c1=beta_est_row(fit_eiie_nex_uee), c1s=beta_se_row(fit_eiie_nex_uee,0), 
                        c2=beta_est_row(fit_eiie_nex_maee), c2s=beta_se_row(fit_eiie_nex_maee,1),
                        c3=beta_est_row(fit_eiie_ed_uee), c3s=beta_se_row(fit_eiie_ed_uee,0),
                        c4=beta_est_row(fit_eiie_ed_maee), c4s=beta_se_row(fit_eiie_ed_maee,1))
xtable(tab2_beta, digits=3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & c1 & c1s & c2 & c2s & c3 & c3s & c4 & c4s \\ 
# \hline
# 1 & 0.452 & 0.126 & 0.452 & 0.126 & 0.378 & 0.121 & 0.378 & 0.122 \\ 
# 2 & 0.453 & 0.123 & 0.453 & 0.124 & 0.395 & 0.120 & 0.395 & 0.120 \\ 
# 3 & 0.467 & 0.126 & 0.467 & 0.126 & 0.431 & 0.122 & 0.431 & 0.122 \\ 
# 4 & 0.465 & 0.127 & 0.465 & 0.128 & 0.463 & 0.123 & 0.463 & 0.123 \\ 
# 5 & 0.395 & 0.133 & 0.395 & 0.134 & 0.425 & 0.125 & 0.425 & 0.125 \\ 
# 6 & 0.315 & 0.144 & 0.315 & 0.144 & 0.374 & 0.129 & 0.374 & 0.129 \\ 
# 7 & 0.236 & 0.156 & 0.236 & 0.156 & 0.322 & 0.134 & 0.322 & 0.134 \\ 
# 8 & 0.123 & 0.163 & 0.123 & 0.164 & 0.216 & 0.135 & 0.216 & 0.136 \\ 
# 9 & 0.053 & 0.168 & 0.053 & 0.169 & 0.158 & 0.136 & 0.158 & 0.137 \\ 
# 10 & 0.033 & 0.168 & 0.033 & 0.169 & 0.123 & 0.135 & 0.123 & 0.136 \\ 
# 11 & 0.012 & 0.168 & 0.012 & 0.169 & 0.103 & 0.135 & 0.103 & 0.136 \\ 
# 12 & 0.394 & 0.125 & 0.394 & 0.126 & 0.244 & 0.074 & 0.244 & 0.074 \\ 
# 13 & -0.003 & 0.174 & -0.003 & 0.175 & 0.075 & 0.172 & 0.075 & 0.173 \\ 
# \hline
# \end{tabular}
# \end{table}

tab2_alpha <- data.frame(c1=alpha_est_row(fit_eiie_nex_uee), c1s=alpha_se_row(fit_eiie_nex_uee,0), 
                         c2=alpha_est_row(fit_eiie_nex_maee), c2s=alpha_se_row(fit_eiie_nex_maee,2),
                         c3=alpha_est_row(fit_eiie_ed_uee), c3s=alpha_se_row(fit_eiie_ed_uee,0),
                         c4=alpha_est_row(fit_eiie_ed_maee), c4s=alpha_se_row(fit_eiie_ed_maee,2))
xtable(tab2_alpha, digits=3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & c1 & c1s & c2 & c2s & c3 & c3s & c4 & c4s \\ 
# \hline
# 1 & 0.471 & 0.024 & 0.475 & 0.025 & 0.487 & 0.024 & 0.492 & 0.025 \\ 
# 2 & 0.392 & 0.027 & 0.396 & 0.027 & 0.939 & 0.010 & 0.940 & 0.010 \\ 
# \hline
# \end{tabular}
# \end{table}


###### CIC
Ustar <- function(y,m,X,beta){
  u <- 1/(1+exp(c(-X%*%beta)))
  v <- u*(1-u)
  t(X)%*%(X*v*m)
}

sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_nex_uee$beta) %*% (fit_eiie_nex_uee$BC1[1:13,1:13]) )) # 17198.82
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_nex_maee$beta) %*% (fit_eiie_nex_maee$BC1[1:13,1:13]) )) # 17198.29
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_ed_uee$beta) %*% (fit_eiie_ed_uee$BC1[1:13,1:13]) )) # 16416.02
sum(diag( Ustar(y=y,m=m,X=X,beta=fit_eiie_ed_maee$beta) %*% (fit_eiie_ed_maee$BC1[1:13,1:13]) )) # 16416.06



#################################################
# Combined results
#################################################

beta_row <- function(fit, BC){
  a <- format(round(fit$beta,3), nsmall=3)
  if (BC==0){
    b <- format(round(sqrt(diag(fit$BC0)[1:13]),3), nsmall=3)
  } else if (BC==1){
    b <- format(round(sqrt(diag(fit$BC1)[1:13]),3), nsmall=3)
  }
  paste0(a," (",b,")")
}
beta_p <- function(fit, BC){
  a <- fit$beta
  if (BC==0){
    b <- sqrt(diag(fit$BC0)[1:13])
  } else if (BC==1){
    b <- sqrt(diag(fit$BC1)[1:13])
  }
  p <- format(round(2*pnorm(-abs(a/b)), 3), nsmall=3)
}
beta_ci <- function(fit, BC){
  a <- fit$beta
  if (BC==0){
    b <- sqrt(diag(fit$BC0)[1:13])
  } else if (BC==1){
    b <- sqrt(diag(fit$BC1)[1:13])
  }
  ci_l <- format(round(a-1.96*b, 3), nsmall=3)
  ci_u <- format(round(a+1.96*b, 3), nsmall=3)
  ci <- paste0("(", ci_l, ", ", ci_u, ")")
}

alpha_row <- function(fit, BC){
  a <- format(round(fit$alpha,3), nsmall=3)
  if (BC==0){
    b <- format(round(sqrt(diag(fit$BC0)[14:15]),3), nsmall=3)
  } else if (BC==2){
    b <- format(round(sqrt(diag(fit$BC2)[14:15]),3), nsmall=3)
  }
  paste0(a," (",b,")")
}
alpha_p <- function(fit, BC){
  a <- fit$alpha
  if (BC==0){
    b <- sqrt(diag(fit$BC0)[14:15])
  } else if (BC==2){
    b <- sqrt(diag(fit$BC2)[14:15])
  }
  p <- format(round(2*pnorm(-abs(a/b)), 3), nsmall=3)
}
alpha_ci <- function(fit, BC){
  a <- fit$alpha
  if (BC==0){
    b <- sqrt(diag(fit$BC0)[14:15])
  } else if (BC==2){
    b <- sqrt(diag(fit$BC2)[14:15])
  }
  ci_l <- format(round(a-1.96*b, 3), nsmall=3)
  ci_u <- format(round(a+1.96*b, 3), nsmall=3)
  ci <- paste0("(", ci_l, ", ", ci_u, ")")
}

tabm_beta <- data.frame(c1=beta_row(fit_eiie_ed_uee,0), c2=beta_p(fit_eiie_ed_uee,0), c3=beta_ci(fit_eiie_ed_uee,0), 
                    c4=beta_row(fit_aie_ed_uee,1), c5=beta_p(fit_aie_ed_uee,1), c6=beta_ci(fit_aie_ed_uee,1))
xtable(tabm_beta)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllllll}
# \hline
# & c1 & c2 & c3 & c4 & c5 & c6 \\ 
# \hline
# 1 & 0.378 (0.121) & 0.002 & ( 0.140, 0.615) & 0.374 (0.122) & 0.002 & ( 0.136, 0.612) \\ 
# 2 & 0.395 (0.120) & 0.001 & ( 0.160, 0.629) & 0.391 (0.120) & 0.001 & ( 0.155, 0.626) \\ 
# 3 & 0.431 (0.122) & 0.000 & ( 0.193, 0.670) & 0.435 (0.123) & 0.000 & ( 0.195, 0.676) \\ 
# 4 & 0.463 (0.123) & 0.000 & ( 0.222, 0.704) & 0.482 (0.123) & 0.000 & ( 0.241, 0.722) \\ 
# 5 & 0.425 (0.125) & 0.001 & ( 0.181, 0.670) & 0.479 (0.122) & 0.000 & ( 0.240, 0.718) \\ 
# 6 & 0.374 (0.129) & 0.004 & ( 0.122, 0.627) & 0.464 (0.122) & 0.000 & ( 0.225, 0.703) \\ 
# 7 & 0.322 (0.134) & 0.016 & ( 0.060, 0.583) & 0.457 (0.121) & 0.000 & ( 0.219, 0.694) \\ 
# 8 & 0.216 (0.135) & 0.109 & (-0.048, 0.481) & 0.378 (0.118) & 0.001 & ( 0.147, 0.608) \\ 
# 9 & 0.158 (0.136) & 0.246 & (-0.109, 0.425) & 0.336 (0.116) & 0.004 & ( 0.109, 0.563) \\ 
# 10 & 0.123 (0.135) & 0.361 & (-0.142, 0.388) & 0.301 (0.115) & 0.009 & ( 0.075, 0.527) \\ 
# 11 & 0.103 (0.135) & 0.446 & (-0.162, 0.367) & 0.280 (0.115) & 0.015 & ( 0.054, 0.506) \\ 
# 12 & 0.244 (0.074) & 0.001 & ( 0.099, 0.389) & 0.064 (0.031) & 0.037 & ( 0.004, 0.125) \\ 
# 13 & 0.075 (0.172) & 0.662 & (-0.262, 0.412) & 0.082 (0.173) & 0.635 & (-0.257, 0.421) \\ 
# \hline
# \end{tabular}
# \end{table}

tabm_alpha <- data.frame(c1=alpha_row(fit_eiie_ed_uee,0), c2=alpha_p(fit_eiie_ed_uee,0), c3=alpha_ci(fit_eiie_ed_uee,0), 
                         c4=alpha_row(fit_aie_ed_uee,2), c5=alpha_p(fit_aie_ed_uee,2), c6=alpha_ci(fit_aie_ed_uee,2))
xtable(tabm_alpha)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllllll}
# \hline
# & c1 & c2 & c3 & c4 & c5 & c6 \\ 
# \hline
# 1 & 0.487 (0.024) & 0.000 & (0.440, 0.535) & 0.487 (0.024) & 0.000 & (0.439, 0.535) \\ 
# 2 & 0.939 (0.010) & 0.000 & (0.920, 0.959) & 0.939 (0.010) & 0.000 & (0.920, 0.959) \\ 
# \hline
# \end{tabular}
# \end{table}

tabs_beta <- data.frame(c1=beta_row(fit_eiie_nex_uee,0), c2=beta_p(fit_eiie_nex_uee,0), c3=beta_ci(fit_eiie_nex_uee,0), 
                    c4=beta_row(fit_aie_nex_uee,1), c5=beta_p(fit_aie_nex_uee,1), c6=beta_ci(fit_aie_nex_uee,1))
xtable(tabs_beta)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllllll}
# \hline
# & c1 & c2 & c3 & c4 & c5 & c6 \\ 
# \hline
# 1 &  0.452 (0.126) & 0.000 & ( 0.206, 0.699) & 0.443 (0.126) & 0.000 & ( 0.197, 0.690) \\ 
# 2 &  0.453 (0.123) & 0.000 & ( 0.211, 0.694) & 0.424 (0.123) & 0.001 & ( 0.182, 0.666) \\ 
# 3 &  0.467 (0.126) & 0.000 & ( 0.221, 0.713) & 0.435 (0.127) & 0.001 & ( 0.187, 0.683) \\ 
# 4 &  0.465 (0.127) & 0.000 & ( 0.216, 0.714) & 0.417 (0.130) & 0.001 & ( 0.163, 0.672) \\ 
# 5 &  0.395 (0.133) & 0.003 & ( 0.134, 0.656) & 0.382 (0.133) & 0.004 & ( 0.121, 0.642) \\ 
# 6 &  0.315 (0.144) & 0.028 & ( 0.033, 0.596) & 0.322 (0.138) & 0.019 & ( 0.052, 0.593) \\ 
# 7 &  0.236 (0.156) & 0.129 & (-0.069, 0.541) & 0.318 (0.134) & 0.018 & ( 0.054, 0.581) \\ 
# 8 &  0.123 (0.163) & 0.450 & (-0.196, 0.442) & 0.247 (0.132) & 0.061 & (-0.011, 0.505) \\ 
# 9 &  0.053 (0.168) & 0.752 & (-0.276, 0.383) & 0.204 (0.130) & 0.117 & (-0.051, 0.460) \\ 
# 10 &  0.033 (0.168) & 0.845 & (-0.297, 0.363) & 0.185 (0.130) & 0.156 & (-0.071, 0.440) \\ 
# 11 &  0.012 (0.168) & 0.945 & (-0.318, 0.341) & 0.163 (0.130) & 0.213 & (-0.093, 0.418) \\ 
# 12 &  0.394 (0.125) & 0.002 & ( 0.148, 0.640) & 0.236 (0.072) & 0.001 & ( 0.095, 0.378) \\ 
# 13 & -0.003 (0.174) & 0.985 & (-0.345, 0.338) & 0.014 (0.176) & 0.937 & (-0.331, 0.359) \\ 
# \hline
# \end{tabular}
# \end{table}

tabs_alpha <- data.frame(c1=alpha_row(fit_eiie_nex_uee,0), c2=alpha_p(fit_eiie_nex_uee,0), c3=alpha_ci(fit_eiie_nex_uee,0), 
                         c4=alpha_row(fit_aie_nex_uee,2), c5=alpha_p(fit_aie_nex_uee,2), c6=alpha_ci(fit_aie_nex_uee,2))
xtable(tabs_alpha)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rllllll}
# \hline
# & c1 & c2 & c3 & c4 & c5 & c6 \\ 
# \hline
# 1 & 0.471 (0.024) & 0.000 & (0.423, 0.518) & 0.470 (0.024) & 0.000 & (0.422, 0.518) \\ 
# 2 & 0.392 (0.027) & 0.000 & (0.339, 0.445) & 0.391 (0.027) & 0.000 & (0.338, 0.445) \\ 
# \hline
# \end{tabular}
# \end{table}


