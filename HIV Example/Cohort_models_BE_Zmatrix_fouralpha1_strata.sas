/*
	file : h:\PCORI-SW\Paper3-Tutorial\Cohortdata\Cohort_models_BE_Zmatrix_fouralpha1_strata.sas
	name : John Preisser
	date: 11/13/2024
    Modified: John Preisser
    what: run GEE/MAEE for cohort example in tutorial paper - HIV analysis Table 6, with strata

4 parameters in correlation model, ICC for different persons within period, and 3-Dependence ICC within persons

    input: h:\PCORI-SW\Paper3-Tutorial\Cohortdata\hivtesting.sas7bdat
    output: h:\PCORI-SW\Paper3-Tutorial\Cohortdata\Cohort_models_BE_Zmatrix_fouralpha1.rtf
*/

%include 'h:\PCORI-SW\Paper-MAEE_Macro\MAEEV2.04.sas';
LIBNAME dat "h:\PCORI-SW\Paper3-Tutorial\Cohortdata";
ods exclude enginehost;


ods rtf style=journal file="H:\PCORI-SW\Paper3-Tutorial\Cohortdata\Cohort_models_BE_Zmatrix_fouralpha1_strata.rtf";

data hivtesting;
    set dat.hivtesting;
run;
proc sort data=hivtesting out = hivtest;
   by clusternum ID time;
run;
proc print data=hivtest(obs=35);
    title "data is hivtesting, first 35 observations";
	id clusternum ID time;
run;
* hivt: has individual undergone HIV testing in the previous three months;
proc freq data=hivtest;
    table Shandong*cluster/nopercent nocol;
run;

data a;
    set hivtest;
	by clusternum ID time;
	if first.clusternum then do; personcnt=0; obsnum=0; end;
	if first.ID then personcnt+1;
	obsnum+1;
*	if personcnt le 75 then output;
	drop condition;
run;

Proc print data=a(obs=35);
    id Shandong cluster obsnum personcnt ID;
run;

*alpha1 (z1) - the ICC for different subjects within the same period; 
*alpha2 (z2) - the ICC for different subjects across different time periods; 
*alpha3 (z3) - the ICC for the same subject across different time points;

proc sql;
	create table zd as
	select
		a1.obsnum as zpair1,
		a2.obsnum as zpair2,
		a1.hivt as hivt1,
		a2.hivt as hivt2,
		a1.clusternum as clusternum1,
        a2.clusternum as clusternum2,
        a1.ID as ID1,
		a2.ID as ID2,
		a1.time as time1,
		a2.time as time2,
		(a1.clusternum = a2.clusternum and a1.ID ne a2.ID and a1.time = a2.time) as z1,
		(a1.clusternum = a2.clusternum and a1.ID = a2.ID and (a2.time - a1.time = 1)) as z2,
		(a1.clusternum = a2.clusternum and a1.ID = a2.ID and (a2.time - a1.time = 2)) as z3,
		(a1.clusternum = a2.clusternum and a1.ID = a2.ID and (a2.time - a1.time = 3)) as z4
		
	from
		a as a1,
		a as a2
	where
		a1.obsnum < a2.obsnum and
                a1.clusternum eq a2.clusternum;
quit;

proc sort data=zd;
    by clusternum1 zpair1 zpair2;
run;
proc print data=zd(obs=100);
    title "data is zmatrix for hiv testing, first 300 observations, in proper sort order";
    id clusternum1 zpair1 zpair2;
run;
proc summary data=zd;
    var z1 z2 z3 z4;
	output out=numcolsz
	   sum = sumz1 sumz2 sumz3 sumz4;
run;
proc print data=numcolsz;
run;


* add PRINTRANGE=NO to suppress printing range violations except when variance computed;

 title "GEE/MAEE estimation of logistic model for binary outcome hiv testing. Correlation model with 3-dependence user-defined z-matrix and identity link";
%GEEMAEE(xydata=a, yvar= hivt, ytype = binary, link =logit, xvar= Shandong period1 period2 period3 period4 intervention, maxiter=100, PRINTRANGE=NO, corrlink = identity, df_choice=2,
    clusterID= clusternum, Zdata=zd, zvar = z1 z2 z3 z4, zpair = zpair1 zpair2, makevone=NO, makephione=YES, alpadj=MAEE, ESTOUT=parmests2, VAR_BETA=BC1, VAR_ALPHA=BC2);

proc print data=parmests2;
run;
proc print data=Var_beta;
run;
proc print data=Var_alpha;
run;

* Var_beta contains the BC1 covariance matrix for the beta estimates;
* create data of variances on one observation after removing covariances;
data varbeta;
   set Var_beta;
   if _n_=1 then do; variable="shandong    ";     variance=shandong; end;
   if _n_=2 then do; variable="period1     ";     variance=period1;  end;
   if _n_=3 then do; variable="period2     ";     variance=period2;  end;
   if _n_=4 then do; variable="period3     ";     variance=period3;  end;
   if _n_=5 then do; variable="period4     ";     variance=period4;  end;
   if _n_=6 then do; variable="intervention"; variance=intervention; end;
   keep variable variance;
run;
*proc print data=varbeta;
*run; 
data varalpha;
   set Var_alpha;
   if _n_=1 then do; variable="icc_parm1"; variance=icc_parm1; end;
   if _n_=2 then do; variable="icc_parm2"; variance=icc_parm2; end;
   if _n_=3 then do; variable="icc_parm3"; variance=icc_parm3; end;
   if _n_=4 then do; variable="icc_parm4"; variance=icc_parm4; end;
   keep variable variance;
run;
*proc print data=varalpha;
*run; 
proc transpose data=parmests2 out=estimates;
run;
*proc print data=estimates;
*run;
data allvariances;
   set varbeta varalpha;
run;
*proc print data=allvariances;
*run;
* Use t with I-2 = 6 degrees of freedom - there are 8 clusters;
data conflimits;
  merge estimates(rename=(Col1=estimate)) allvariances;
  tcrit = tinv(0.975,6); 
  se_est = sqrt(variance);
  lower_CI = estimate - tcrit*se_est;
  upper_CI = estimate + tcrit*se_est; 
  if not(_name_="ICC_parm1") and not(_name_="ICC_parm2") and not(_name_="ICC_parm3") and not(_name_="ICC_parm4") then do;
     OR_est = exp(estimate);
	 lower_OR = exp(lower_CI);
	 upper_OR = exp(upper_CI);
  end;
run;
proc print data=conflimits;
  title "95% confidence limits based on t-distribution with I-2=6 df, BC1 and BC2 standard errors for betas and ICCs, respectively";
run;
ods rtf close;
