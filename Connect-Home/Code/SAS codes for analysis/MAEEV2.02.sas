/**********************************************************************************************************

SAS Macro GEEMAEE version 2

filename: MAEEV2.sas

A SAS macro for matrix-adjusted estimating equations (MAEE), which are finite-sample bias-corrected estimating equations for intra-cluster correlation parameters used in
iteratively reweighted least squares estimation with the standard GEE for marginal mean model parameters, generalizing the uncorrected estimating equations for correlated 
binary outcomes of Prentice (1988), with extensions to correlated continuous and count outcomes.

By Ying Zhang and John Preisser, 2022 Sep
Original version by Bing Lu and John Preisser 2005

Primary references: 

Zhang Y, Preisser JS, Li F, Turner EL, Toles M, Rathouz PJ. GEEMAEE: A SAS macro for the analysis of correlated outcomes based on GEE and finite-sample 
adjustments with application to cluster randomized trials. Computer Methods and Programs in Biomedicine. 2023 Mar;230:107362.
doi: 10.1016/j.cmpb.2023.107362. Epub 2023 Jan 20. 

Preisser JS, Lu B, Qaqish BF. Finite sample adjustments in estimating equations and covariance estimators for 
intracluster correlations. Statistics in Medicine 2008; 27:5764-5785. 

Other references:

Kalema G, Molenberghs G, Kasshun W. Second-order generalized estimating equations for correlated count data. Comput Stat 2016; 31: 749-770.

Li F, Turner EL, Preisser JS. Sample size determination for GEE analyses of stepped wedge cluster randomized trials. Biometrics 2018; 74(4): 1450-1458.

Prentice RL. Correlated binary regression with covariates specific to each binary observation. Biometrics 1988; 44:1033-1048.

Zhang Y, Preisser JS, Turner EL, Rathouz PJ, Toles M, Li F. A general method for calculating power for GEE analysis of complete and incomplete stepped wedge cluster randomized trials. 
Statistical Methods in Medical Research 2023. 32(1), 71-87.



***********************************************************************************************************
Aim: The GEEMAEE SAS macrto Fits choice among three types of estimating equations for 
correlated binary,normal, count random variables when the association structure is of interest. 
***********************************************************************************************************
Edits
***********************************************************************************************************
# Edits by Ying Zhang on May 2021:
# add options of type of the outcomes,(binary,normal and count) and available link functions for the marginal mean model and marginal model for correlation
# add modules to calculate the marginal mean (linkid), marginal variance(CREATEL) and the derivative to beta(CREATEC) of three outcomes with different link functions. 
Range checks
# for binary outcomes, range of the marginal mean will be checked to be within [0-1]
# for count outcomes, range of the marginal mean for count outcomes will be positive.

# add modules to calculate the marginal mean (CORRLINKID), the derivative to alpha(deri2) of correlation parameters with different link functions. 
# add common scale parameter (PHI) estimated with method of moments for continuous outcome and extra-Poisson variance. For Poisson, default is to estimate phi, but it can be set to 1
# for continuous outcomes, use Identity matrix as default working covariance matrix (W in equation (1), Preisser, Lu and Qaqish) in correlation estimation equations


# for count outcome with assumed Poisson distribution, default is working covariance matrix W for correlations based the bivariate Poisson distribution and equation (8) of Kalema, Molenberghs and Kassahun (2016), setting (default) phi=1
# for simplified count analysis, we recommend to use option makevone=1 with estimated dispersion parameter PHI
# Add the finite-sample variance adjustment FG/BC3 in the version
# add option for creating an output data set containing parameter estimates and covariance matrices for BC0, KC and MD, FG

# delete the option of wdata and wvar in the MAEEV1.01.sas, instead calculate the cluster sample size vector, n by the cluster id column in the XYdata

# change the name of the SAS macro from PRENTICEcorrJSP to GEEMAEE
# edit the output into more long table format


# Version 1.07 Edits on Dec,2021 by Ying Zhang
# add inner correlation structures in the SAS macro
# Z data is optional in the version
# period ID is required for built in correlation structures and subject ID is required for the cohort design.
# add the (y1,y2) pair check into the user-defined Zdata, in the check_Z_pairs module;
# add the diagnostics in the SAS macro on the cluster level and subject level;

# Edits on Feb,2022 by Ying Zhang
# Version 1.08 Edits on May,2022 by Ying Zhang
# add the diagnostics in the SAS macro on the cluster period level and examine it on the real dataset;
# Change the inverse function (INV) into the generalized inverse to reduce the convergence error due to the ICCs parameters out of the convergence;


# Version 2.0  Edits on Sep,2022 by Ying Zhang
# Version MAEEV2 Edits on Sep,2022 by Ying Zhang
# Output the results of the estimates and BC0-BC3 in datasets
# All deletion diagnostics use the bias-corrected variance estimator(MB/BC0-BC3) for beta and alpha 

# Version 2.01 Edits on Jan, 2023 by Ying Zhang
# Edit the codes for the BC3 calculation for the senario with cluster size=1
# Version 2.02 Edits on March, 2023 by Ying Zhang
# Solve bugs about the corr specification error due to %else %if in the correlation structures part;
# replace all %else %if into %if
**************************************************************************************************************
All Inputs of the SAS macro
**************************************************************************************************************


    xydata: The data set containing the outcome and marginal mean covariates 
      yvar: The outcome in xydata
	  ytype: binary, normal or count outcomes
	  link: link function in the population-averaged model,identity /log/ logit/ probit
	  corrlink: link function for the correlation-averaged model ,Identity/  Log/  Logit/ Fishersz(Fisher's Z transformation)
      xvar: A list of marginal mean covariates in xydata. Must all be numeric
      clusterID: Cluster ID (consecutive integers) in xydata, need to be continuous and successive,
      PeriodID: Period ID (consecutive integers) in xydata
      SubjectID: Subject ID (consecutive integers) in xydata
      Corr: Correlation structures in the GEE analysis of multi-period GRTs
      zdata: The data set containing the covariates for the marginal correlation regression (Optional)
      zvar: A list of marginal correlation covariates in zdata. Must all be numeric  (Optional)
      ypair: The Z var should be all pairs j < k for each cluster i, ie (1,2) (1,3).... (1,n)... (2,3).... (n-1, n) for each cluster 
   maxiter: Maximum number of iterations for Fisher Scoring,optional, default at 50
   epsilon: Tolerance for convergence,optional, default at 0.001 
printrange: Print details of range violations? (Should initially check) (YES, NO)
            Range is always checked, however. Range is always checked and printed
            when variance computed.
printfinal: Print final results; 1=yes, 0=no, optional, default is yes.           
    shrink: THETA or ALPHA for binary outcomes. See Readme file for details.
  makevone: Sets VEE = 1 if YES instead of w_ijk (if NO) from Prentice (1988). 
            See Sharples and Breslow (1992).
  makephione: Set the dispersion parameter as 1 
    alpadj: Bias adjustment to alpha estimating equations: 1 - no adjustment (UEE,  Prentice, 1988;
            2 - scalar adjustment (SAEE, Sharples and Breslow (1992); 3 -
        Matrix-based I-H adjustment (MAEE, Preisser, Lu and Qaqish 2008), Default is UEE.       
    clsout: output data set of cluster deletion diagnostics for beta and alpha (not required)
    obsout: output data set of observation deletion diagnostics for beta and alpha (not required)
    cpsout: output data set of cluster period deletion diagnostics for beta and alpha (not required)
    ESTOUT: output the beta and alpha estimates into a dataset 
    VAR_BETA: output the variance estimators for parameters in the marginal mean model into a dataset, MB/BC0/BC1/BC2/BC3
    VAR_ALPHA: output the sandwich variance estimators for parameters in the marginal correlation model into a dataset, BC0/BC1/BC2/BC3

    

Cluster ID's should be integers from 1 to K. 

BIAS-CORRECTED VARIANCE ESTIMATOR OF M & D (2001), K and C (2001) and choice of three different bias
corrections to alpha estimating equations including Preisser and Lu's matrix-based adjustment



*******************************************************************************************************/



%macro GEEMAEE(xydata=, yvar=, ytype =, link = ,corrlink=, xvar=, clusterID=, PeriodID=, SubjectID=,Corr =, zdata=, zvar=, zpair=,
maxiter=50, epsilon=0.0001, printrange=YES, 
printfinal = 1, shrink=ALPHA, makevone=, makephione=, alpadj=, CLSOUT=, OBSOUT=, CPOUT=, ESTOUT = ,VAR_BETA=,VAR_ALPHA=); 


/**************************************************************************************

Aim:     calculate the cluster size for each cluster

**************************************************************************************/
                                                                                                                                                                                                    
   proc means noprint data = &xydata;
      var &yvar;
      by &clusterID;
      output out = _sample(keep = n) n = n;
      run;

   proc iml worksize =4000000;
    show space;

  
/**************************************************************************************

MODULE: GRABDATA
Aim:     Reads the data xydata dataset for analysis

OUTPUT
	      y: The outcomes
	      X: Marginal mean covariates
	     id: Cluster identifier
	      n: Vector of cluster sample sizes
	      Z: Marginal correlation covariates
	maxiter: Max number of iterations
	epsilon: Tolerence for convergence
	   nobs: The total count of obs in the dataset xydata
  clusterid: The cluster id indicator column in the xydata
# update remove the wdata and w vector in the SAS macro GEECORR.
**************************************************************************************/
   start GRABDATA(y, X, n, Z,  maxiter, epsilon ,nobs,clusterid,Ypair_error);
      maxiter = &maxiter;
	  epsilon = &epsilon;
      use &xydata NOBS nobs;
         read all var{&yvar} into y;
         read all var{&xvar} into X;
         read all var{&clusterID} into clusterid;
         *read all var{&periodID} into periodid;
      close &xydata;

      nobs = nrow(X);
      use _sample;
         read all var{n} into n;
      close _sample;
         
 /*If user has their user-defined correlation structures with Z data, read zdata in the reader */
 %if %sysevalf(%superq(Corr)=,boolean) %then %do;
	      use &zdata;
	         read all var{&zvar} into Z;
	         read all var{&zpair} into Ypairs_user;
	      close &zdata;
	  /*Check the Ypairs_user with the real used pair of observations in the correlations*/
	      run check_Z_pairs(n,Ypairs_user,Ypair_error);            
      %end;/*the user provide */
      
       %else %do;
          run  Generate_Z_corr(n, Z); /*Z is the SAS macro defined Z matrix*/
       %end;
     
      finish GRABDATA;
      
      

/**************************************************************************************
MODULE: GRAB_periodsID
Aim: Reads in the period ID for deletion diagnostics at cluster period level 
    and generate the period pair for the Z matrix for the built-in Correlation structures definition

Input:
		n: The cluster size vector
OUTPUT
      periodid: Output period ID
      Zpair_periods: Pairs of the periods for the row of Z matrix
     
**************************************************************************************/
   start GRAB_periodsID(n, periodid, Zpair_periods);

    run BEGINEND(firstx, lastx, n);
* Generate the pair of observations of the Z matrix;
     use &xydata NOBS nobs;
         read all var{&periodID} into periodid;
     close &xydata;
     Zpair_periods = {};
      
     do i = 1 to nrow(n); 
     period_id_i = periodid[firstx[i]:lastx[i],];
 	 Zpair_periods_i = {};/*when n[i]^=1*/
   	  if n[i]>1 then do;/*empty Z pair when n[i]=1*/
	   	       l =1;
	   	       Zpair_periods_i = J(NCHOOSE2(n[i]),2,1);/*the z vector for pairwise correlation in the ith cluster has size n[i]*(n[i]-1)/2  */
	    	   do j = 1 to (n[i]-1);
	            	do k = (j+1) to n[i];	
	            	      Zpair_periods_i[l,1] = period_id_i [j];/*the period of the jth observation of the lth pair of correlation in the ist cluster*/
	                      Zpair_periods_i[l,2] = period_id_i [k];/*the period of the jth observation of the lth pair of correlation in the ist cluster*/	 
	            	      l = l+1;
	            	 end;
	            end; /*j*/
	     end;/*n[i]>1 */
       Zpair_periods = Zpair_periods//Zpair_periods_i;    
   	  end;/*i=1 loop*/ 

     
      finish GRAB_periodsID;



   
/**************************************************************************************

MODULE: check_Z_pairs
	Check the user-provided pair of observations in Ypair for the Z matrix with the defined pair of observations in the SAS macro
	if Ypairs_user = Ypair built-in then the user provides the correct order of the pair of correlation.
	Use when Corr is blank

INPUT
      n: Vector of cluster sample sizes
      Ypairs_user: the Z pair provided by the user

OUTPUT       
       Ypair_error: 1: user-provided pair of observations in Ypair for the Z matrix is the same 
       with the defined pair of observations in the SAS macro;
       The software will be terminated if the Ypair_error=1;
       0- errors in the definition of the Ypairs_user
 
**************************************************************************************/


START  check_Z_pairs(n, Ypairs_user,Ypair_error);  
    Ypair_error=0;
    run BEGINEND(firstx, lastx, n);
     Ypair = {};
     do i = 1 to nrow(n);   	    
   	       l =1;
   	       Ypair_i = J(NCHOOSE2(n[i]),2,1);
    	   do j = 1 to (n[i]-1);
            	do k = (j+1) to n[i];
            	      Ypair_i[l,1] = j;/*the jth observation of the lth pair of correlation*/
            	      Ypair_i[l,2] = k;/*the kth observation of the lth pair of correlation*/ 
            	      l = l+1;
            	 end;
            end; /*j*/
           Ypair = Ypair//Ypair_i;    
   	  end;/*i=1 loop*/ 

   	
   if min(Ypair = Ypairs_user) = 0 then do;
       print "Error: The order of Ypair is incorrect in at least one cluster compared to orders in the SAS macro";
       Ypair_error = 1;
   end;
   
  
  
finish check_Z_pairs;


   
/**************************************************************************************
MODULE: Generate_Z_corr
Generate ZvaR matrix for the five built-in correlation structures:
EX:  Simple Exchangeable 
NE: Nested Exchangeable
ED: Exponential Decay
BE: Block Exchangeable 
PD: Proportional Decay

INPUT
      xydata: The data set containing the outcome and marginal mean covariates 
      n : the number of cluster size in each  cluster


OUTPUT
       Z: The ZVAR matrix for all the observations by a built-in correlation structures
       
 
**************************************************************************************/


START  Generate_Z_corr(n, Z);

    run BEGINEND(firstx, lastx, n);
    Z = {};
    /*For simple exchangeable correlation, there is only one ICC with Z as one */
   
   %if %upcase(&Corr)= EX %then %do;
     do i = 1 to nrow(n);
      Z_i = {}; /*empty Z pair when n[i]=1*/
      if n[i]>1 then do; 
	   	        Z_i = J(NCHOOSE2(n[i]),1,1);
	   	        Z = Z//Z_i;
	   end;/*n[i]>1*/
   	  end;/*i=1 loop*/ 
   	  %if %sysevalf(%superq(corrlink)=,boolean) %then %let corrlink = IDENTITY;

   %end; /*Corr="EX"*/
   	
   	/*cross-section design correlation*/
   	/*For exponential decay correlation, ICCs decay with the periods between obs*/
   %if %upcase(&Corr)=ED %then %do;
      use &xydata NOBS nobs;
         read all var{&periodID} into periodid;
      close &xydata;
    
     do i = 1 to nrow(n);  
        period_id_i = periodid[firstx[i]:lastx[i],];
        Z_i = {};
   	      if n[i]>1 then do;
	   	       l =1;
	   	       Z_i = J(NCHOOSE2(n[i]),2,1);
	    	   do j = 1 to (n[i]-1);
	            	do k = (j+1) to n[i];
	            	      Z_i[l,2] = abs(period_id_i[j] - period_id_i[k]);
	            	      l = l+1;
	            	 end;
	            end; /*do loop of j*/
	        end; /*n[i]>1*/
           Z = Z//Z_i;	        
   	  end;/*i=1 to n loop*/ 
   	 %if %sysevalf(%superq(corrlink)=,boolean) %then %let corrlink = LOG;
   	 *print "ED" ;
   %end; /*Corr="ED"*/
   	
 /*For nested correlation structures, ICCs differs from within cluster and the across clusters */  	
    %if %upcase(&Corr)=NE %then %do;
      use &xydata NOBS nobs;
         read all var{&periodID} into periodid;
      close &xydata;
    
     do i = 1 to nrow(n);
        Z_i = {};
        period_id_i = periodid[firstx[i]:lastx[i],];
        *print  i period_id_i;
     if n[i]>1 then do;	  	
   	       l =1;
   	       Z_i = J(NCHOOSE2(n[i]),2,1);
   	       *print Z_i;
    	   do j = 1 to (n[i]-1);
            	do k = (j+1) to n[i];
            	      Z_i[l,1] = (period_id_i[j] = period_id_i[k]);
            	      Z_i[l,2] = (period_id_i[j] ^= period_id_i[k]);
            	      l = l+1;
            	 end;
            end; /*j*/
      end;  /*n[i]>1*/ 
      Z = Z//Z_i;
    end;/*i=1 loop*/ 
   	 %if %sysevalf(%superq(corrlink)=,boolean) %then %let corrlink = Identity;
    *print "NE";  
   	%end; /*Corr="NE"*/
   	
 /*Cohort design*/  
/*For Proportional decay correlation structures, ICCs differs from within cluster and the across clusters */  	  
    %if %upcase(&Corr)=PD %then %do;
      use &xydata NOBS nobs;
         read all var{&periodID} into periodid;
         read all var{&subjectID} into subjectid;
      close &xydata;
    
     do i = 1 to nrow(n);
     
        period_id_i = periodid[firstx[i]:lastx[i],];
        subjectid_i = subjectid[firstx[i]:lastx[i],];
    	if i=1 then do;
    	   l =1;
    	   Z = J(NCHOOSE2(n[1]),2,1);
    	   do j = 1 to (n[1]-1);
            	do k = (j+1) to n[1];
            	      Z[l,1] = (subjectid_i[j] ^= subjectid_i[k]); /*exist for all the obeservations from different subjects*/
            	      Z[l,2] = abs(period_id_i[j] - period_id_i[k]);
            	      l = l+1;
            	 end;
            end; /*j*/
    	end;/*i=1*/
    	
    	
   	    else do; /*when i^=1*/
   	        l =1;
   	       Z_i = J(NCHOOSE2(n[i]),2,1);
    	   do j = 1 to (n[i]-1);
            	do k = (j+1) to n[i];
            	      Z_i[l,1] = (subjectid_i[j] ^= subjectid_i[k]);
            	      Z_i[l,2] = abs(period_id_i[j] - period_id_i[k]);
            	      l = l+1;
            	 end;
            end; /*j*/
           Z = Z//Z_i;
   	    end;/*when i^=1*/
   	    
   	    
   	  end;/*i=1 loop*/ 
   	%if %sysevalf(%superq(corrlink)=,boolean) %then %let corrlink = LOG;
   	%end; /*Corr="PD"*/
/*For Block exchangeable correlation structures, ICCs differs from within cluster and the across clusters,and within subject */  	  
   
    %if %upcase(&Corr)=BE %then %do;
      use &xydata NOBS nobs;
         read all var{&periodID} into periodid;
         read all var{&subjectID} into subjectid;
      close &xydata;
    
     do i = 1 to nrow(n);
     
        period_id_i = periodid[firstx[i]:lastx[i],];
        subjectid_i = subjectid[firstx[i]:lastx[i],];
    	if i=1 then do;
    	   l=1;
    	   Z = J(NCHOOSE2(n[1]),3,1);
    	   do j = 1 to (n[1]-1);
            	do k = (j+1) to n[1];
            	      Z[l,1] = (period_id_i[j] = period_id_i[k]);
            	      Z[l,2] = (subjectid_i[j] ^= subjectid_i[k])&(period_id_i[j] ^= period_id_i[k]);
            	      Z[l,3] = (subjectid_i[j] = subjectid_i[k])&(period_id_i[j] ^= period_id_i[k]);
            	      l = l+1;
            	 end;
            end; /*j*/
    	end;/*i=1*/
    	
    	
   	    else do; /*when i^=1*/
   	       l =1;
   	       Z_i = J(NCHOOSE2(n[i]),3,1);
    	   do j = 1 to (n[i]-1);
            	do k = (j+1) to n[i];
            	      Z_i[l,1] = (period_id_i[j] = period_id_i[k]);
            	      Z_i[l,2] = (subjectid_i[j] ^= subjectid_i[k])&(period_id_i[j] ^= period_id_i[k]);
            	      Z_i[l,3] = (subjectid_i[j] = subjectid_i[k])&(period_id_i[j] ^= period_id_i[k]);
            	      l = l+1;
            	 end;
            end; /*j*/
           Z = Z//Z_i;
   	    end;/*when i^=1*/
   	    
   	    
   	  end;/*i=1 loop*/ 
   	%if %sysevalf(%superq(corrlink)=,boolean) %then %let corrlink = LOG;
   	%end; /*Corr="BE"*/
 
  finish Generate_Z_corr;


           
/**************************************************************************************

MODULE: BEGINEND
Creates two vectors that have the start and end points for each cluster

INPUT
n: Vector of cluster sample sizes

OUTPUT
first: Vector with starting row for cluster i
 last: Vector with ending row for cluster i

**************************************************************************************/

   start BEGINEND(first, last, n);
      last = cusum(n);
      first = last - n + 1;
   finish BEGINEND;
   


/**************************************************************************************
Module: IS_POS_DEF
    A = symmetric matrix
    returns 1 if A is positive definite
              0 otherwise
**************************************************************************************/

start is_pos_def (A);
	ev = eigval(A);
	*print ev;
    return (min(ev[,1]) >=0);
finish;


/**************************************************************************************

MODULE: NCHOOSE2
For a value n, returns n choose 2

INPUT
n: Vector of cluster sample sizes

OUTPUT
n choose 2

*************************************************************************************/

   start NCHOOSE2(n);
      return((n#n - n)/2);
      finish NCHOOSE2;
      
/************************************************************************************
MODULE: FINDINV
Calculate the inverse for symmetric and positive finite matrix
************************************************************************************/
   start FINDINV (A);
   AHALF=ROOT(A);
   GINV=ginv(AHALF);
   AINV=GINV*GINV`;
   return (AINV);
   finish;


/**************************************************************************************

MODULE: GETROWB
Generates a row of the -E(d(corr)/dbeta) (ie Q) matrix

INPUT
   mu: marginal means for cluster i
gamma: pairwise correlation for outcome j and k for cluster i
    j: indicator for mean 
    k: indicator for mean 
    X: covariate matrix for cluster i
    y: The outcome: binary, count or normal;
    s2: The dispersion parameters, assumed 1 for binary outcomes and not restricted for normal and count outcomes

OUTPUT
Row of -E(d(corr)/dbeta) (ie Q) matrix
  #      Note that the original expression omits a term in taking the derivatives, 
  #      hence not exactly the same as the 1st equation in Page 1040
  #      The current module added back the term in Page 1040
  #      Match the calculation of GETROWB in R package geeCRT
# added with the option of three type of responses
# after V2, edit the row calculation for count outcomes with the consideration of overdisperstion parameter, S2
**************************************************************************************/

   start GETROWB(mu, gamma, j, k, X, y,s2);
 /*The calculation of the C_j, V_j, C_k, V_k has be taken into account the link function and outcome types*/  			
		C_j = CREATEC(X[j,], mu[j]);
		C_k = CREATEC(X[k,], mu[k]);
		V_j = CREATEL(mu[j],s2);
		V_k = CREATEL(mu[k],s2);
   		%IF &ytype = binary  %THEN %DO; 
   			row = gamma/2*((1-2*mu[j])*C_j/V_j+(1-2*mu[k])*C_k/V_k+
                     2*C_j/(y[j]-mu[j])+
                     2*C_k/(y[k]-mu[k])); 
        * row = (gamma/2) * ( (1-2*mu[j])*X[j,] + (1-2*mu[k])*X[k,] );/*The one in Trace's codes, not sure which one to use*/
        %end;
		%IF &ytype = normal  %THEN %DO;
		    if (y[j]-mu[j])=0 | (y[k]-mu[k])=0 then row =J(1,%sysfunc(countw(&xvar)),0);/*if there is the case with y=mu, the code could not run and we will set rowb=0*/
			else row = gamma*(C_j/(y[j]-mu[j])+C_k/(y[k]-mu[k])); 
			
		%end;
		%IF &ytype = count  %THEN %DO;

			row = gamma/2*(s2*C_j/V_j+s2*C_k/V_k+
                     2*C_j/(y[j]-mu[j])+
                     2*C_k/(y[k]-mu[k])); %end; /* edition after V2,consider the overdispersion parameter into the count outcomes;*/
                                                /* only impact the variance estimations of correlation parameters */
		
      return(row);
    finish GETROWB;



/**************************************************************************************
MODULE: LINKID
Define link function.  This module defines a function that applies the appropiate link
to the linear predictor.  The linear predictor (lp), sometimes also called eta in various
textbooks, defines the relationship to the explanatory variables and regression
coefficients (X\beta).

Meanwhile, the link function describes on how the mean (mu=E[Y|X]) is related to the linear
predictor.

INPUT
link:	1 - Identity
		2 - Logarithm
		3 - Logit
	  	4 - Probit
  lp: Linear predictor/eta/Xbeta

OUTPUT
  the expected mean, by the transformation of the linear predictor with the selected link applied
**************************************************************************************/
start linkid(lp); 
    %IF %upcase(&link) = IDENTITY %THEN %DO;
            return(lp); 
    %END;
    
    %IF %upcase(&link) = LOG %THEN %DO;
            return(exp(lp)); 
    %END;
    
    %IF %upcase(&link) = LOGIT %THEN %DO;
            return(exp(lp)/(1+exp(lp))); 
    %END;

    %IF %upcase(&link) = PROBIT %THEN %DO; 
            return(probnorm(lp));
    %END;
finish linkid;



/**************************************************************************************

MODULE: CREATEA
Creates residual for beta estimating equation, (Y - mu)

INPUT
mu: Vector of n_i marginal means
 y: Outcome vector for ith cluster

OUTPUT
residuals for beta estimating equation

**************************************************************************************/

   start CREATEA(mu, y);
      return(y - mu);
      finish CREATEA;

/**************************************************************************************

MODULE: CREATEC
Creates derivative matrix for beta estimating equation, dmu/dbeta

INPUT
 X: Covariate matrix for cluster i
mu: Vector of n_i marginal means


&link:	1 - Identity
		2 - Logarithm
		3 - Logit
	  	4 - Probit
  lp: Linear predictor/eta/Xbeta

OUTPUT
derivative matrix for beta estimating equation

**************************************************************************************/

   start CREATEC(X, mu);

  	%IF %upcase(&link) = IDENTITY %THEN %DO;
  		nobs = nrow(X);
        return(X); /*identity link*/
    %END;
    
    %IF %upcase(&link) = LOG %THEN %DO; /*log link*/
        return(X#mu);  
    %END;
    
    %IF %upcase(&link)= LOGIT %THEN %DO; /*logit link*/;
        return(X#(mu#(1-mu)));  
    %END;

	%IF %upcase(&link) = PROBIT %THEN %DO; /*probit link*/;
		lp = Quantile('NORMAL',mu);
        return(X#pdf('NORMAL', lp));  
    %END;
      
   finish CREATEC;


/**************************************************************************************

MODULE: CREATEL
	Creates variance function for different outcome for uncorrected variable

INPUT
 	ytype: the outcome type: binary, normal, and count 
	mu: Vector of n_i marginal means
	s2: dispersion parameters

OUTPUT
	derivative matrix for beta estimating equation

**************************************************************************************/

   start CREATEL(mu, s2);

  	%IF %upcase(&ytype) = BINARY %THEN %DO;
  		if min(mu)<0 | max(mu)>1 then do; minmu = min(mu);maxmu = max(mu);
  		print minmu maxmu  "Error: The mean of the binary outcome out of bound [0,1]"; end;
  		else do; s2=1; return(diag(mu#(1-mu))); end;
    %END;
    
    %IF %upcase(&ytype) = NORMAL %THEN %DO;
    	nobs = nrow(mu);
        return(diag(J(nobs, 1, 1)#s2)); 
    %END;
    
    %IF %upcase(&ytype) = COUNT %THEN %DO;
    	if min(mu)<0  then do; print "The mean of the count outcome is negative"; end;
  		else do; return(diag(mu#s2)); end; 
    %END;
    *print mu s2;
      
   finish CREATEL;

/**************************************************************************************

MODULE: CREATEB
Creates covariance matrix for beta estimating equation, var(Y)

INPUT
   mu: Vector of n_i marginal means
    n: Sample size (scalar) for cluster i
   s2: Dispersion parameter 
gamma: Vector of rho_jk between outcome j and k for cluster i

OUTPUT
covariance matrix for beta estimating equation

**************************************************************************************/

   start CREATEB(mu, gamma, n, s2);
	  %IF %upcase(&ytype) = BINARY  %THEN %DO;
		  d = s2*mu # (1 - mu);
	      B = diag(d);
	      l = 1;
	      do j = 1 to n - 1;
	         do k = j + 1 to n;
	            B[j,k] = s2*sqrt(mu[j] * mu[k] * (1-mu[j]) * (1-mu[k])) * gamma[l] ;
	            l = l + 1;
	            end;
	         end;
	      return(sqrsym(symsqr(B`)));
	  %end;

	  %IF %upcase(&ytype) = NORMAL  %THEN %DO;
	  	  d = s2*J(n,1,1);
	      B = diag(d);
	      l = 1;
	      do j = 1 to n - 1;
	         do k = j + 1 to n;
	            B[j,k] = s2 * gamma[l] ;
	            l = l + 1;
	            end;
	         end;
	      return(sqrsym(symsqr(B`)));
			     
	  %END;

	  %IF %upcase(&ytype) = COUNT  %THEN %DO;
	  	  d =  s2*mu;
	      B = diag(d);
	      l = 1;
	      do j = 1 to n - 1;
	         do k = j + 1 to n;
	            B[j,k] = s2 *sqrt(mu[j] * mu[k]) * gamma[l] ;
	            l = l + 1;
	            end;
	         end;
	      return(sqrsym(symsqr(B`)));
			     
	  %END;
      finish CREATEB;



/**************************************************************************************
MODULE: CORRLINKID
Define correlation link function.  This module defines a function that applies the appropiate link
to the correlation linear predictor.  The linear predictor (lp2), sometimes also called eta2 in various
textbooks, defines the relationship to the explanatory variables and regression
coefficients (Z\alpha).

Meanwhile, the link function describes on how the correlation (rho=E[R|Z]) is related to the linear
predictor.

INPUT
link:	1 - Identity
		2 - Logarithm
		3 - Logit
		4 - Fisher's Z Transformation
  lp2: Linear predictor2/eta2/Zalpha

OUTPUT
  Transformation of the linear predictor with the selected link applied
**************************************************************************************/
start corrlinkid(lp2); 
    %IF %upcase(&corrlink) = IDENTITY %THEN %DO;
            return(lp2); 
    %END;
    
    %IF %upcase(&corrlink) = LOG %THEN %DO;
            return(exp(lp2)); 
    %END;

    %IF %upcase(&corrlink) = LOGIT %THEN %DO;
            rhoi = exp(lp2); 
            return(rhoi/(1+rhoi)); 
    %END;

	%IF %upcase(&corrlink) = FISHERSZ %THEN %DO;
            return((exp(lp2)-1)/(exp(lp2)+1)); 
    %END;
finish;


/**************************************************************************************
MODULE: DERI2
This module defines a function that applies the derivative of the correlation (rho_i) with 
respect to the linear predictor (lp2 or eta2_i) based on the selected link for the correlation.

\frac{\partial(\rho_i)}{\partial(\lp2_i)}
	= Z (\frac{\partial(\lp2_i)}{\partial(\rho_i)})^{-1}

INPUT
link:	1 - Identity
		2 - Logarithm
		3 - Logit
		4 - Fisher's Z Transformation
OUTPUT
  Derivative of the correlation (rho_i) with respect to the linear predictor (lp2 or eta2_i)
  based on the link selected (see above)		
**************************************************************************************/
start deri2(Z , rhoi); 
    %IF %upcase(&corrlink) = IDENTITY %THEN %DO;
        return(Z ); 
    %END;
    
    %IF %upcase(&corrlink) = LOG %THEN %DO;
        return(Z # rhoi);  
    %END;

    %IF %upcase(&corrlink) = LOGIT %THEN %DO;
        return(Z # (rhoi#(1-rhoi)));  
    %END;

    %IF %upcase(&corrlink) = FISHERSZ %THEN %DO;
        return(Z # ((1+rhoi)#(1-rhoi))/2);  
    %END;
finish;




/**************************************************************************************

MODULE: INITBETA
Generates initial values for beta. Approximates regression without clustering effects
also need to consider the link function

INPUT
y: The outcome
X: Marginal mean covariates
n: Vector of cluster sample sizes


OUTPUT
beta: Vector of marginal mean parameters without assuming no correlation among the group
ustar: Approximate information matrix
s2: Initiate the dispersion parameter
# add the normal outcome and count outcome with canonical link
**************************************************************************************/

   start INITBETA(beta, ustar, y, X, n,s2);

 		  beta = solve(X` * X, X` * y);
	      do i = 1 to 2;
		  	 lp = X*beta;
		  	 mu  = linkid(lp);
		  	 %IF &ytype = normal %THEN %DO;
			      s2 = sum((y-mu)`*(y-mu))/(sum(n)-ncol(X));   
			  %END;
			  
			 %else %DO; s2 =1;%end;
			 *print s2;
		  	 L = CREATEL(mu, s2);
			 C = CREATEC(X, mu); /*is the derivative matrix D , dmu/dbeta*/

	         U = C`*inv(L) * (y - mu);
	         ustar = C` * inv(L) * C;
	         d = solve(ustar, U);
	         beta = beta + d;
	      end;


       finish INITBETA; 



/**************************************************************************************
MODULE: INITVALS
Generates initial values for beta and correlation parameter alpha. 
Uses 0.1 for initial value for alphas. Beta uses approximate parms from a logistic
regression on all observations, treating them independently

INPUT
y: The outcome
X: Marginal mean covariates
n: Vector of cluster sample sizes
q: Number of marginal correlation parameters

OUTPUT
 beta: Vector of marginal mean parameters
alpha: Vector of marginal correlation parameters
ustar: Approximate information matrix
**************************************************************************************/

   start INITVALS(beta, ustar, alpha, y, X, n, q,s2);
      
      
    %IF %upcase(&corrlink) = IDENTITY %THEN %DO;
        alpha = J(q, 1, 0.01);
    %END;
    
    %IF %upcase(&corrlink) = LOG %THEN %DO;
    	alpha = J(q,1,0);
         alpha[1] = -1;
         alpha[2] = -0.1;
    %END;

    %IF %upcase(&corrlink) = LOGIT %THEN %DO;
         alpha = J(q,1,log(0.01/0.99));   
    %END;

	%IF %upcase(&corrlink) = FISHERSZ %THEN %DO;
         alpha = J(q,1,log(0.99/1.01));   
    %END;
      run INITBETA(beta, ustar, y, X, n, s2);
   finish INITVALS;

/**************************************************************************************
MODULE: VEERDB
Embedded within SCORE MODULE to do the following:
(i) create VEE, R, and DB matrices in the estimating equations for beta and alpha
(ii) check to see if correlation is within the range based on marginal means for binary outcomes

INPUT
    gamma_c: correlation vector for cluster
    mu_c:    mean vector for cluster
    X_c:     Covariates for beta estimating equations
    y_c:     outcome vector in the cluster 
    flag:    Performs an Eigenanalysis of variance estimation of Y (B) to see if positive definite. Only
                 called when computing the variance at the end (0 = no, 1 = yes). Prints
                 warning for each cluster violation.
    n:   cluster size in the ith cluster
    i:   index for cluster
    p:    dimension of beta
    s2:   Dispersion parameter
    Hi1,Ci,Gi: The matrix in the SAEE and MAEE for finite sample adjustment

OUTPUT
    VEE: working variances for alpha estimating equations
    R:   residual vector for alpha estimating equations
    DB:  -E(d(corr)/dbeta) (i.e. Q) matrix
    rangeflag: Checks to see if correlation is within Frechet bound based on marginal means for binary outcomes
           (0 = in range, 1 = out of range). See Prentice (1988).
**************************************************************************************/
start VEERDB(VEE,R,DB,rangeflag,mu_c,gamma_c,X_c,y_c,flag,n,p, s2,Hi1,Ci,Gi,E);

	VEE = J(NCHOOSE2(n), 1, 0);
	R = J(NCHOOSE2(n), 1, 0);
	DB = J(NCHOOSE2(n), p, 0);

	l = 1;
	do j = 1 to n - 1;/*the empirical covariance is calculated according to the order of y in cluster and period, note that z should be calculated followed the order of y*/
	   do k = j + 1 to n;

	%if &ytype = binary %then %do; 	
		if mu_c[j]<0 | mu_c[j]>1 | mu_c[k]<0 | mu_c[k]>1 then do;
			margmeanj = mu_c[j];
	        margmeank = mu_c[k];
			print "Marginal Mean for binary outcomes outside of 0-1 Range Detected for Pair:",  i j k margmeanj margmeank;
			stop;
		end;
	
	   rho_upper_limit = min(sqrt((mu_c[j]*(1-mu_c[k]))/(mu_c[k]*(1-mu_c[j]))), 
	                       sqrt((mu_c[k]*(1-mu_c[j]))/(mu_c[j]*(1-mu_c[k])))); 
	
	   rho_lower_limit = max(-sqrt((mu_c[j]*mu_c[k])/((1-mu_c[j])*(1-mu_c[k]))), 
	                       -sqrt(((1-mu_c[j])*(1-mu_c[k]))/(mu_c[j]*mu_c[k]))); 
	
	  * RANGE CHECKS (Frechet bound) ;
		*when flag = 0, these are range checks that are printed during the estimation iterations - this is only printed when printrange = YES;
	   if ((gamma_c[l] >= rho_upper_limit) | (gamma_c[l] <= rho_lower_limit))
	        & flag = 0 then do;
	       rangeflag = 1;
	
	       %if %upcase(&printrange) = YES %then %do; 
	           print "Range Violation Detected for Pair:",    j k;
		       CovarJ = X_c[j,];
		       CovarK = X_c[k,];
		       print covarj covark;
	
	           corrjk = gamma_c[l];
	           margmeanj = mu_c[j];
	           margmeank = mu_c[k];
	           print rho_lower_limit corrjk rho_upper_limit margmeanj margmeank;
	       %end;/*%upcase(&printrange) = YES*/
	
	   end;/*flag = 0 and Frechet bound is violated*/
	
	   *when flag = 1, these range checks are only printed for the final estimation iteration - this is always printed;
	   if ((gamma_c[l] >= rho_upper_limit) | (gamma_c[l] <= rho_lower_limit))
	        & flag = 1 then do;
			rangeflag = 1;
	
	       print "Range Violation Detected for Pair:",    j k;
	
	   end;/*flag = 1 and Frechet bound is violated*/
	   
	 %end;/*Frechet bound check for %if = binary outcomes*/
	  
	   %if %upcase(&makevone) = YES %then %do;
	       VEE[l] = 1; 
	   %end; /*makeone = 1*/
	
	   %else %do; /*makeone ^= 1, generate the exact variance of residual*/
		   	%if  %upcase(&ytype) = BINARY %then %do;
		       		VEE[l] = 1 + ((1-2*mu_c[j])*(1-2*mu_c[k])*gamma_c[l]) /sqrt(mu_c[j]*(1-mu_c[j])*mu_c[k]*(1-mu_c[k])) - gamma_c[l]**2; 	* variance of sample correlation - Prentice (1988) p. 1039; 
		     %end;
		   	
		   	%IF %upcase(&ytype) = NORMAL  %THEN %DO;
	               VEE[l] = 1 + gamma_c[l]**2;
	  		%end;
	
			%IF %upcase(&ytype) = COUNT  %THEN %DO;
	         		VEE[l] = Corr_Var_Poisson(mu_c[j], mu_c[k],gamma_c[l]);
	 		%end;

		%end;/*else makeone =NO*/
		
		/*generate the R matrix in the correlation estimating equations*/
		 %if  %upcase(&ytype) = BINARY %then %do;
	   			  %if %upcase(&alpadj) = UEE %then %do;
                      R[l] = ((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) / sqrt(mu_c[j]*(1-mu_c[j])*mu_c[k]*(1-mu_c[k])) - gamma_c[l]; *UEE;
                  %end;       

                  %if %upcase(&alpadj) = SAEE %then %do; *       Scalar-adjustment of Sharples and Breslow (1992);
                      R[l] = (1/((1- Hi1[j,j])*(1- Hi1[k,k])))*((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) / sqrt(mu_c[j]*(1-mu_c[j])*mu_c[k]*(1-mu_c[k])) - gamma_c[l];
                  %end;  

                   %if %upcase(&alpadj) = MAEE %then %do; *       Matrix-based multiplicative correction, (I - H_i)^{-1}, Preisser, Lu and Qaqish, 2008;
                      R[l] = Ci[j, ]*Gi[ ,k] - gamma_c[l];            
                   %end;
	   	%end;
	   	
	   	%IF %upcase(&ytype) = NORMAL %THEN %DO;
                *       Uncorrected estimating equations (UEE) of Prentice (1988);
                       %if %upcase(&alpadj) = UEE %then %do;
                              R[l] = ((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) / s2 - gamma_c[l];
                       %end;       
				*       Scalar-adjustment of Sharples and Breslow (1992);
                       %if %upcase(&alpadj) = SAEE %then %do;
                              R[l] = (1/((1- Hi1[j,j])*(1- Hi1[k,k])))*
                                     ((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) / s2 - gamma_c[l];
                       %end;  
				*       Matrix-based multiplicative correction, (I - H_i)^{-1}, Preisser, Lu and Qaqish, 2008;
                       %if %upcase(&alpadj) = MAEE %then %do;
                              R[l] = Ci[j, ]*Gi[ ,k] - gamma_c[l];      
                        %end;

		%end;

		%IF %upcase(&ytype) = COUNT  %THEN %DO;
           *       Uncorrected estimating equations (UEE) of Prentice (1988);
                       %if %upcase(&alpadj) = UEE %then %do;
                              R[l] = ((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) /  sqrt(mu_c[j]*mu_c[k])/s2 - gamma_c[l];
                       %end;       
		  *       Scalar-adjustment of Sharples and Breslow (1992);
                       %if %upcase(&alpadj) = SAEE %then %do;
                              R[l] = (1/((1- Hi1[j,j])*(1- Hi1[k,k])))*((y_c[j] - mu_c[j])*(y_c[k] - mu_c[k])) / sqrt(mu_c[j]*mu_c[k])/s2 - gamma_c[l];
                       %end;  
		  *       Matrix-based multiplicative correction, (I - H_i)^{-1}, Preisser, Lu and Qaqish, 2008;
                       %if %upcase(&alpadj) = MAEE %then %do;
                              R[l] = Ci[j, ]*Gi[ ,k] - gamma_c[l];            
                       %end;

		%end;

		
	    /* insert check that variance is nonnegative */
              if VEE[l] <= 0 then do;
                     VEEFLAG=1; 
                     print "The variance of the empirical variance estimator is negative for Pair:",    j k;
                     goto leaveVEE;
              end;
         DB[l,] = GETROWB(mu_c, gamma_c[l], j, k, X_c,y_c,s2); 							* row of the -E(d(corr)/dbeta);
	   	 l = l + 1;
	   end;/*j*/
	end;/*i*/

leaveVEE:finish VEERDB;


/**************************************************************

MODULE: INVBIG
compute (A - mm`)^{-1}c without performing the inverse directly
INPUT
   ainvc: inverse of matrix A times vector c
   ainvm: inverse of matrix A times matrix (with low no. of columns) M
      M: matrix of eigen column vectors m1,m2, ..
      c: vector 
  start: of do loop
  end:   of do loop, rank of X
**************************************************************/

Start INVBIG(ainvc,ainvm,m,c,start,end);
	do i=start to end;
	   b = ainvm[,i];
	   bt = b`;
	   btm = bt*m;
	   btmi = btm[,i];
	   gam = 1 - btmi;
	   bg = b/gam;
	   ainvc = ainvc + bg*(bt*c);
	   if i< end then ainvm = ainvm  + bg*btm;  
	end;
finish INVBIG;

/**************************************************************

MODULE: PHI
Generates moment estimates for dispersion parameter with finite sample  adjustment
INPUT
   Ustar: inverse of matrix A times vector c
	beta: Vector of marginal mean parameters estimated in the kth step 
 alpha: Vector of marginal correlation parameters estimated in the kth step 
     y: The outcome
     X: Marginal mean covariates
     Z: Marginal correlation covariates
     n: Vector of cluster sample sizes
     p: Number of marginal mean parameters
     q: Number of marginal correlation parameters
    s2: The estimated dispersion parameter in the kth step 
   

 
  start: of do loop
  end:   of do loop, rank of X
**************************************************************/



START  PHI(Ustar, beta, alpha, y, X, Z, n, p,q,s2);
    RSS=0;
    naive=ginv(Ustar[1:p,1:p]);
    
    run BEGINEND(firstx, lastx, n);
    run BEGINEND(firstz, lastz, NCHOOSE2(n));
    
    do i = 1 to nrow(n);
           X_c = X[firstx[i]:lastx[i],];
           y_c = y[firstx[i]:lastx[i]];
           lp_c = X_c*beta;
           *print Lp_c;
		   mu_c = linkid(lp_c);
		   Z_c = Z[firstz[i]:lastz[i],];
           lp_alpha = Z_c * alpha;
           gamma_c = corrlinkid(lp_alpha);
		    
           if n[i] = 1 then do;
              B = CREATEL(mu_c, s2);
              R_i= ginv(B)*(y_c - mu_c);
              RSS = RSS + s2*R_i**2;

           end;  /* do loop, n[i]=1 */

		else do;/* do loop, n[i]^=1 */
		   C = CREATEC(X_c, mu_c);
           B = CREATEB(mu_c, gamma_c, n[i],s2);
           L_c = CREATEL(mu_c, s2);
           A = CREATEA(mu_c, y_c);

           *INVB=FINDINV(B);
           INVB=GINV(B);
           CtinvB = C`*invB; 
           Hi1=C*naive*CtinvB;
		   R_i = inv(sqrt(L_c))*(y_c - mu_c);/*uncorrected residual vector for the uth clusters, just related to be assumed mean and variance (L matrix)*/
			
			%if %upcase(&alpadj) = UEE %then %do; /*For UEE*/			
				 
				 RSS = RSS + s2*R_i`*R_i; 
			%end;
			
			%if %upcase(&alpadj) = SAEE %then %do; /*For SAEE*/
				do j=1 to n[i];
                    	RSS = RSS + s2*(R_i[j]/(1- Hi1[j,j]))**2;
                end;
		
			%end;
			
              %if %upcase(&alpadj) = MAEE %then %do;
					IminusH = diag(J(n[i],1,1)) - Hi1;
                    Ci = inv(sqrt(L_c))*ginv(IminusH)*sqrt(L_c);
                    Gi=  R_i*R_i`; 
                    do j=1 to n[i];
                    		RSS = RSS + s2*Ci[j, ]*Gi[ ,j];
                    end;
                 %end;
		end;/* do loop, n[i]^=1 */
  end;/*i in 1:nrow(n)*/

    
    if RSS<0 then do;/*only happen with MAEE, */
    RSS= 0;
	           do i = 1 to nrow(n);
	           X_c = X[firstx[i]:lastx[i],];
	           y_c = y[firstx[i]:lastx[i]];
	           lp_c = X_c*beta;
			   mu_c = linkid(lp_c);
			   Z_c = Z[firstz[i]:lastz[i],];
	           lp_alpha = Z_c * alpha;
	           gamma_c = corrlinkid(lp_alpha);
			    
	           if n[i] = 1 then do;
	              B = CREATEL(mu_c, s2);
	              R_i= ginv(B)*(y_c - mu_c);
	              RSS = RSS + s2*R_i**2;
	
	           end;  /* do loop, n[i]=1 */
	
			else do;/* do loop, n[i]^=1 */
			   C = CREATEC(X_c, mu_c);
	           B = CREATEB(mu_c, gamma_c, n[i],s2);
	           L_c = CREATEL(mu_c, s2);
	           A = CREATEA(mu_c, y_c);
	
	           *INVB=FINDINV(B);
	           INVB=GINV(B);
	           CtinvB = C`*invB; 
	           Hi1=C*naive*CtinvB;
			   R_i = inv(sqrt(L_c))*(y_c - mu_c);/*uncorrected residual vector for the uth clusters, just related to be assumed mean and variance (L matrix)*/
				

						IminusH = diag(J(n[i],1,1)) - Hi1;
	                    Ci = inv(sqrt(L_c))*ginv(IminusH)*sqrt(L_c);
	                    Gi=  R_i*R_i`; 
	                    do j=1 to n[i];
	                    		RSS = RSS + s2*abs(Ci[j, ]*Gi[ ,j]);
	                    end;/*s2>0*/

			end;/* do loop, n[i]^=1 */
	  end;/*i in 1:nrow(n)*/

 end;
   s2=RSS/(sum(n)-p);      
   *print s2;
    leavephi:finish PHI;



  /**************************************************************
  # MODULE: Corr_Var_Poisson
  # AIM:  calculate the variance of residual of count responses by Kalema, assumed to have no dispersion, s2=1 in the case 
  # INPUT :
  # mui: the ith mean of the poisson response
  # muj: the jth mean of the poisson response
  # gamma_ij: the working correlation correspinding to the pair of ith and jth subjects.
  # Output: v_rirj, the exact variance of empirical pair-wise correlation vector for Poisson outcomes
  **************************************************************/
	Start  Corr_Var_Poisson (mui, muj,gamma_ij);
	    muij = gamma_ij* sqrt(mui* muj);/* The covariance between Yi and Yj;*/
	    E_yiyj_square = mui**2*muj**2  +  mui*muj**2 +
	      mui**2* muj + 2*muij**2 + mui* muj*(1+4*muij)+2*muij*(mui+ muj)+ muij;
	    E_yi_square_yj = (mui**2+mui)*muj + 2*muij*mui + muij;
	    E_yj_square_yi =  mui*(muj**2+muj) + 2*muij*muj + muij;
	    E_yiyj =   gamma_ij*sqrt(mui*muj) + mui*muj;
	    E_yi_square = mui +mui**2;
	    E_yj_square = muj +muj**2;
	    E_rirj_square = (E_yiyj_square - 2*muj*E_yi_square_yj- 2*mui*E_yj_square_yi + muj**2*E_yi_square+
	                       mui**2*E_yj_square +4*mui*muj*E_yiyj - 3*mui**2*muj**2)/(mui*muj);
	    v_rirj = E_rirj_square - gamma_ij**2;
	    
	return (v_rirj);
	  
	  
    finish Corr_Var_Poisson;

/**************************************************************************************

MODULE: SCORE
Generates the score matrix for each cluster and approximate information
to be used to estimate parameters and generate standard errors

INPUT
     beta: Vector of marginal mean parameters
    alpha: Vector of marginal correlation parameters
        y: The outcome
        X: Marginal mean covariates
        Z: Marginal correlation covariates
        n: Vector of cluster sample sizes
        p: Number of marginal mean parameters
        q: Number of marginal correlation parameters
	   s2: Dispersion parameter, s2 = 1 at default for binary outcomes but need to be estimated for normal outcome and count outcome.
     flag: Performs an Eigenanalysis of B to see if positive definite. Only
           called when computing the variance at the end (0 = no, 1 = yes). Prints
           warning for each cluster violation.
rangeflag: Checks to see if correlation is within range based on marginal means
           (0 = in range, 1 = out of range). See Prentice (1988).
VEEFLAG:  Checks if all variances for alpha e.e. are positive, terminates if not
NPSDFLAG: Checks if B is positive definite, terminates if not
NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
clusterid : CLuster indicator

OUTPUT
     U: Score vector
UUtran: Sum of U_i*U_i` across all clusters
 Ustar: Approximate information matrix

 # need changes

**************************************************************************************/
   start SCORE(U, UUtran, Ustar, ustarold, beta, alpha, y, X, Z, n, p,
   q,s2 , flag, rangeflag, VEEFLAG, NPSDFLAG, NPSDADJFLAG,clusterid);
      
      U = J(p+q, 1, 0);
      UUtran = J(p+q, p+q, 0);
      Ustar = J(p+q, p+q, 0);
      naiveold =ginv(Ustarold[1:p,1:p]);  /* for Hi1 below, the A matrix in the Preisser 2008 based on the initial beta estimation */
		
      run BEGINEND(firstx, lastx, n);
      run BEGINEND(firstz, lastz, NCHOOSE2(n));
      
      do i = 1 to nrow(n);
           *print U;
           X_c = X[firstx[i]:lastx[i],];
           y_c = y[firstx[i]:lastx[i]];
           if n[i] = 1 then Z_c = 0;
              else Z_c = Z[firstz[i]:lastz[i],];

			lp_c = X_c*beta;
			mu_c =linkid(lp_c);
			L_c = CREATEL(mu_c, s2);
			lp_alpha = Z_c * alpha;
            gamma_c = corrlinkid(lp_alpha);
            C = CREATEC(X_c, mu_c); /*is the derivative vector D in the Prentice paper*/
            E = deri2(Z_c, gamma_c);/*is the derivative vector E in the Prentice paper*/
           if n[i] = 1 then do;/*only update the beta information and scores*/
              U_c = C`*inv(L_c) * (y_c - mu_c);
              UUtran_c = U_c * U_c`;
              Ustar_c = C` * inv(L_c) * C; /*D`V^-1D*/
              
              U[1:p] = U[1:p] + U_c ;
              UUtran[1:p, 1:p] = UUtran[1:p, 1:p] + UUtran_c;
              Ustar[1:p, 1:p] = Ustar[1:p, 1:p] + Ustar_c ;
           end; /*n[i]=1*/

           else do;/*n[i]^=1*/
                U_c = J(p+q, 1, 0);
                Ustar_c = J(p+q, p+q, 0);
	           
				
                 VEE = J(NCHOOSE2(n[i]), 1, 0);
                 R = J(NCHOOSE2(n[i]), 1, 0);
                 DB = J(NCHOOSE2(n[i]), p, 0);

                 B = CREATEB(mu_c, gamma_c, n[i],s2);/*Creates covariance matrix for var(Y), considering correlation structures*/
                 A = CREATEA(mu_c, y_c);

                 CtinvB = C`*ginv(B); 
                 Hi1=C*naiveold*CtinvB;  

                 %if &alpadj = MAEE %then %do;
					 SQVARFUN = sqrt(L_c);                    
                     INVSQVAR = inv(SQVARFUN);
*                IH=(I(n[i])-Hi1);
*                    Ci = FINDINV(IH);
* IH above needs to be symmetric, but we can use the following;
                     CT = t(C);
                     omega = C*naiveold*CT;
                     vminomega = B - omega; /*approximate Gi^-1*/

                     psd_vmin = is_pos_def(vminomega);
                     mineig = min(eigval(vminomega));
                     Ci = B*ginv(vminomega); /*Gi in the paper*/
                     RX = INVSQVAR*(y_c - mu_c);
                     Gi=RX*RX`; 
                     if psd_vmin=0 then do;                      
                       NPSDADJFLAG=1;     
                     end;   
                 %end;
           *print U;      
                 
			run VEERDB(VEE,R,DB,rangeflag,mu_c,gamma_c,X_c,y_c,1,n[i],p ,s2,Hi1,Ci,Gi,E);

/*			
* Check for pos def of B;
* score module will allow the non psd matrix but with warning;*/ ;
                if flag = 1 then do;
                    if min(eigval(B)) <= 0 then do;

                       print"Warning: in cluster:" i,
                            "Var(Y) is not Positive-Definite";
							NPSDFLAG = 1;                     
                    end;/*min(eigval(B)) <= 0 */
                end;

                 INVB=ginv(B);
                 U_c[1:p] = C` * INVB * A;
                 U_c[p+1:p+q] = E` * (R # (1 / VEE)); /*need changes in the different link functions*/
                 UUtran_c = U_c * U_c`;
                 Ustar_c[1:p, 1:p]  = C` * INVB * C;
                 
* Next line is commented out to give updating equations (17) in Prentice (1988);
* A more detailed procedure is given by un-commenting to get procedure based upon (15) of Prentice;
				
                 Ustar_c[p+1:p+q, 1:p] = E` * (DB # (1 / VEE));/*need changes in the different link functions*/
                 Ustar_c[p+1:p+q, p+1:p+q] = E` * (E # (1 / VEE));/*need changes in the different link functions*/
                

                 U = U + U_c ;
                 UUtran = UUtran + UUtran_c ;
                 Ustar = Ustar + Ustar_c;
                 *print U_c Ustar_c;
             end;  /* else do (cluster size > 1) */
      end;  /* do i = 1 to nrow(n) */
     *print Ustar;
      rangeflag = 0;

leaveScore:finish SCORE;





/**************************************************************************************
MODULE: CIC_Criteria
Generates the Correlation information criteria (CIC) for the specified mean and correlation model

INPUT
        n: Vector of cluster sample sizes
        y: The outcome
        X: Marginal mean covariates
     beta: Vector of marginal mean parameters
   robust: Robust covariance matrix for beta and alpha
        S2: DISPERSION parameters

OUTPUT
	CIC:  correlation information criteria (CIC)
**************************************************************************************/
start CIC_Criteria(CIC, n, y, X, beta, robust, s2);

      run BEGINEND(firstx, lastx, n);
      p =nrow(beta); 
	  invcovind =J(p,p,0) ;
      do i = 1 to nrow(n);
           X_c = X[firstx[i]:lastx[i],];
           y_c = y[firstx[i]:lastx[i]];

		   	   lp = X_c * beta;
			   mu_c = linkid(lp);
			   D_c = CREATEC(X_c, mu_c);		/* partial mu / partial beta*/
			   A_c = CREATEL(mu_c,s2); 		        /* variance for cluster of size 1 */

			   invcovind_c = D_c` * inv(A_c) * D_c;

			   invcovind = invcovind + invcovind_c;
	  end;

	  robustbeta = robust[1:p,1:p]; *empirical covariance matrix for beta, BC0;
	  CIC = trace(invcovind * robustbeta);

finish CIC_Criteria;



/****************************************************************************************
MODULE: MAKEVAR_WT_DIAG
Creates covariance matrix of beta and alpha and output any deletion diagnostics .

INPUT
  beta: Vector of marginal mean parameters
 alpha: Vector of marginal correlation parameters
     y: The 0/1 binary outcome
     X: Marginal mean covariates
     Z: Marginal correlation covariates
     n: Vector of cluster sample sizes
     p: Number of marginal mean parameters
     q: Number of marginal correlation parameters

OUTPUT
robust: Robust covariance matrix for beta and alpha
 naive: Naive (Model-Based) covariance matrix for beta
 varMD, varKC, varFG: BC2, BC1 and BC3 for the bias-corrected sandwich variance 
 # need changes
**************************************************************************************/

   start MAKEVAR_WT_DIAG(robust, naive, varMD, varKC,varFG, beta, alpha, Ustarold, y, X, Z, n, p, q,s2 ,
                 VEEFLAG,ROBFLAG,NPSDADJFLAG, nobs,clusterid);
      
      run SCORE(U, UUtran, Ustar, Ustarold, beta, alpha, y, X, Z, n, p, q,s2, 0, 0,VEEFLAG,NPSDFLAG,NPSDADJFLAG,clusterid);
* U from module SCORE is not needed because the Scores are updated with the lastes estimates;
      naiveboth = GINV(Ustar);
      naive = naiveboth[1:p,1:p];   /* if B=0, this is upper left matrix of inustar, for beta */ 
      naivealp = naiveboth[p+1:p+q,p+1:p+q];  /* for use in cluster diagnostics */
	  robust = GINV(Ustar) * UUtran * GINV(Ustar); *UUtran = LAMBDA matrix from Prentice 15;
      Hbeta    = J(nrow(n),1,0); /* cluster leverage for beta */
      Halpha   = J(nrow(n),1,0); /* cluster leverage for alpha */

* new commands to compute INV(I - H1);
      Call Eigen(evals1,evecs1,naive);
      sqrevals1 = sqrt(evals1);
      sqe1 = evecs1*diag(sqrevals1);

* new commands to compute INV(I - H2);
      Call Eigen(evals2,evecs2,naivealp);
      sqrevals2 = sqrt(evals2);
      sqe2 = evecs2*diag(sqrevals2);

     /***************************
     * Bias-corrected variance estimators BC1-BC3 in GEEMAEE paper, as well in Li 2018*
     ****************************/

      UUtran = J(p+q, p+q, 0);
      UUbc = J(p+q, p+q, 0);
      UUbc2 = J(p+q, p+q, 0);
      UUbc3 = J(p+q, p+q, 0);
      Ustar = J(p+q, p+q, 0);
      inustar = J(p+q, p+q, 0);
      Ustar_c_array=J(nrow(n), (p+q)*(p+q), 0);
      UUtran_c_array=J(nrow(n), (p+q)*(p+q), 0);
 
      run BEGINEND(firstx, lastx, n);
      run BEGINEND(firstz, lastz, NCHOOSE2(n));
 
      do i = 1 to nrow(n);
           X_c = X[firstx[i]:lastx[i],];
           y_c = y[firstx[i]:lastx[i]];
           lp_c = X_c*beta;
		   mu_c = linkid(lp_c);
		   
        if n[i] = 1 then do;
              B = CREATEL(mu_c, s2);
              C = CREATEC(X_c, mu_c);
              U_i= C`*inv(B)*(y_c - mu_c);
              UUtran_c = U_i * U_i`;
              UUtran[1:p, 1:p] = UUtran[1:p, 1:p] + UUtran_c;
              *print U_i UUtran_c;
              Ustar_c = C` * inv(B) * C;;     
              Ustar[1:p, 1:p] = Ustar[1:p, 1:p] + Ustar_c ;

              Hi1=C*naive*C`*inv(B);
              U_c = C` *inv(B)* (y_c - mu_c)/(1 - Hi1);
              UUbc_c = U_c * U_c`;      
              UUbc[1:p, 1:p] = UUbc[1:p, 1:p] + UUbc_c ;/*BC2*/
              UUbc_ic = U_c * U_i`;      
              UUbc2[1:p, 1:p] = UUbc2[1:p, 1:p] + UUbc_ic ;/*BC1*/
		      Ustar_c_array[i,1:p*p]=rowvec(Ustar_c);/*for BC3, only update variance estimator for beta*/   
	          UUtran_c_array[i,1:p*p]=rowvec(UUtran_c);/*for BC3, only update variance estimator for beta*/
            
        end;  /* do loop, n[i]=1 */
        
        else do; /* do loop, n[i] ^=1 */

           U_i=J(p+q, 1, 0);
           U_c = J(p+q, 1, 0);
           Ustar_c = J(p+q, p+q, 0); 
           Z_c = Z[firstz[i]:lastz[i],];
           lp_alpha = Z_c * alpha;
           gamma_c = corrlinkid(lp_alpha);
		   E = deri2(Z_c, gamma_c);
* commands for beta;
           C = CREATEC(X_c, mu_c);
           B = CREATEB(mu_c, gamma_c, n[i],s2);
           L_c = CREATEL(mu_c, s2);
           A = CREATEA(mu_c, y_c);
           INVB=ginv(B);
           U_i[1:p]=C`*invB*A;
           
           CtinvB = C`*invB; 
           Hi1=C*naive*CtinvB;
           Hbeta[i] = trace(Hi1);  /* cluster leverage for beta */
*           IH=(I(n[i])-Hi1);
*           Ci=FINDINV(IH);

* commands for generalized inverse - beta;
           ai1 = invB;
           mm1 = C*sqe1;
           ai1A=ai1*A;
           ai1m1=ai1*mm1;
           run INVBIG(ai1A,ai1m1,mm1,A,1,p);/*calculate the (I - H^(0.5)H^(0.5))^(-1)*A*/
           U_c[1:p] = C`*ai1A;

* commands for alpha;
           VEE = J(NCHOOSE2(n[i]), 1, 0);
           R = J(NCHOOSE2(n[i]), 1, 0);
           DB = J(NCHOOSE2(n[i]), p, 0);
           RX = J(n[i], 1, 0);

                 %if %upcase(&alpadj) = MAEE %then %do;
                     SQVARFUN = sqrt(L_c);
                     INVSQVAR = inv(SQVARFUN);
*                IH=(I(n[i])-Hi1);
*                    Ci = FINDINV(IH);
* IH above needs to be symmetric, but we can use the following;
                     CT = t(C);
                     vminomega = B - C*naive*CT;
                     psd_vmin = is_pos_def(vminomega);
                     mineig = min(eigval(vminomega));
                     Ci = B*GINV(vminomega);
                     RX = INVSQVAR*(y_c - mu_c);
                     Gi=RX*RX`;

                        if psd_vmin=0 then do; 
	                       NPSDADJFLAG=1; 
	                    end;   
                 %end;

		  run VEERDB(VEE,R,DB,rangeflag,mu_c,gamma_c,X_c,y_c,1,n[i],p ,s2,Hi1,Ci,Gi,E);
 			     
           U_i[p+1:p+q]=E`*(R#(1/VEE));/*All Z_c need changes*/
           mm2 = E*sqe2;
           ai2R = R/VEE; 
           ai2m2 = mm2#(1/VEE);
           run INVBIG(ai2R,ai2m2,mm2,R,1,q);
           U_c[p+1:p+q] = E` *ai2R;
           
			 %IF (%length(&CLSOUT) ^= 0) %THEN %DO; 
			/* cluster leverage for alpha */                     
					  Zctvee = E`#(1/VEE`);
			          Halpha[i] = trace(E*naivealp*Zctvee);  /* preserves memory */
			* note trace(ABC) = trace(CAB), Searle (1982), p. 45;
			
			          free Zctvee;            
			          free Hi2;
			%end;

           Ustar_c[1:p, 1:p]  = C` * INVB * C;
           Ustar_c[p+1:p+q, 1:p] = E` * (DB # (1 / VEE));
           Ustar_c[p+1:p+q, p+1:p+q] = E` * (E # (1 / VEE));
        
           Ustar = Ustar + Ustar_c;

           UUtran_c = U_i * U_i`;
           UUtran = UUtran + UUtran_c ;
           UUbc_c = U_c * U_c`;
           UUbc = UUbc + UUbc_c ;
           UUbc_ic=U_c*U_i`;
           UUbc2 = UUbc2 + UUbc_ic ;
		   Ustar_c_array[i,]=rowvec(Ustar_c);/*for BC3,record the Fi information in the ith cluster*/   
	       UUtran_c_array[i,]=rowvec(UUtran_c);/*for BC3*/
       
        end;  /*     do loop n[i] > 1  */
        end;     /*   do i = 1 to nrow(n) */

*The former  loop is for the varKC/Bias corrected variance calculation;
     inustar[1:p,1:p] = GINV(Ustar[1:p,1:p]);
     inustar[p+1:p+q, p+1:p+q] = GINV(Ustar[p+1:p+q, p+1:p+q]);
     inustar[p+1:p+q,1:p] = -inustar[p+1:p+q, p+1:p+q]*Ustar[p+1:p+q, 1:p]*inustar[1:p,1:p];
     inustartr = inustar`;
*print inustar UUtran;
** BC0 or usual Sandwich estimator of Prentice (1988);     
     robust = inustar * UUtran * inustartr;
*print UUbc2;
** BC2 or Variance estimator that extends Mancl and DeRouen (2001);
     varMD = inustar * UUbc * inustartr;

** BC1 or Variance estimator that extends Kauermann and Carroll (2001);
     varKC = inustar * (UUbc2+UUbc2`) * inustartr/2;
     *print robust varMD varKC;
** BC3 or Variance estimator that extends Kauermann and Carroll (2001); 
    * refer to Fan Li
    * calculating adjustment factor for BC3;
      
      do i =1 to nrow(n);
       Fi = J((p+q),(p+q),0);
         if n[i]>1 then do ;
              
		      DVD_omegainv = diag(shape(Ustar_c_array[i,], (p+q), (p+q))*inustar) ;	  
				do j =1 to (p+q);
					Fi[j,j] = 1/sqrt(1-min(0.75,DVD_omegainv[j,j])); 
				end;
              UUbc3=UUbc3+Fi*shape(UUtran_c_array[i,], (p+q), (p+q))*Fi;
          end;
          else if n[i]=1 then do ;/*only provide info on beta estimation*/
		       
                 DVD_omegainv = diag(shape(Ustar_c_array[i,1:p*p], p, p)*inustar[1:p,1:p]) ;	  
				do j =1 to p;
					Fi[j,j] = 1/sqrt(1-min(0.75,DVD_omegainv[j,j])); 
				end;
              UUbc3[1:p,1:p]=UUbc3[1:p,1:p]+Fi[1:p,1:p]*shape(UUtran_c_array[i,1:p*p], p, p)*Fi[1:p,1:p];/*just update the beta part for BC3*/
              
          end;
      end;
     varFG=inustar*UUbc3*inustartr;


        if min(vecdiag(robust))  <= 0 then ROBFLAG = 1;
        if min(vecdiag(varMD))  <= 0 then ROBFLAG = 1;
        if min(vecdiag(varKC))  <= 0 then ROBFLAG = 1;
        if min(vecdiag(varFG))  <= 0 then ROBFLAG = 1;
    
*For deletion diagnostics in the second loop;

      /*for any level deletion diagnostics*/
      *calculate the size of deletion statistics;
   Run_DIAG = (%length(&CLSOUT)^= 0 | %length(&CPOUT)^= 0| %length(&OBSOUT)^= 0);
   IF Run_DIAG=1 then do; 
    

    %if %length(&periodID)^=0 %then %do;
             run GRAB_periodsID(n, periodid, Zpair_periods);
             *print periodid;
    %end;/*finish the deletion diagnostics calculation*/

      dfbetcls = J(nrow(n),p,0);  /* dfbeta for a cluster */
      dfalpcls = J(nrow(n),q,0);  /* dfalpha for a cluster */
      cookdcls = J(nrow(n),1,0); /* cooks distance overall */ 
      cookbetcls = J(nrow(n),1,0); /* cooks distance for beta */
      cookalpcls = J(nrow(n),1,0);  /* cooks distance for alpha */
      naivecookdcls = J(nrow(n),1,0);       /* cooks distance overall */ 
      naivecookbetcls = J(nrow(n),1,0);     /* cooks distance for beta */
      naivecookalpcls = J(nrow(n),1,0);     /* cooks distance for alpha */

      /*for observation level deletion diagnostics*/
      Dfbetobs=J(nobs,p,0);
      Dfalpobs = J(nobs,q,0);
      cookdobs = J(nobs,1,0);
      CookBetObs = J(nobs,1,0);
      CookAlpObs = J(nobs,1,0);
      naivecookdobs = J(nobs,1,0);  
      naivecookbetobs = J(nobs,1,0);
      naivecookalpobs = J(nobs,1,0);

      cluster_periods_count = J(nrow(n),1,0);
      do i = 1 to nrow(n);
         cluster_periods_count[i] = ncol(unique(periodid[firstx[i]:lastx[i]]));
      end;
      
      run BEGINEND(firstcp, lastcp, cluster_periods_count);

      
      Dfbetcp=J(lastcp[nrow(n)],p,0);
      Dfalpcp = J(lastcp[nrow(n)],q,0);
      cookdcp = J(lastcp[nrow(n)],1,0);
      CookBetcp = J(lastcp[nrow(n)],1,0);
      CookAlpcp = J(lastcp[nrow(n)],1,0);
      naivecookdcp = J(lastcp[nrow(n)],1,0);  
      naivecookbetcp = J(lastcp[nrow(n)],1,0);
      naivecookalpcp = J(lastcp[nrow(n)],1,0);
      clusterid_CP = J(lastcp[nrow(n)],1,0);
      periodid_CP = J(lastcp[nrow(n)],1,0);
      CP_size= J(lastcp[nrow(n)],1,0);
 
      
 


* norms for diagnostics;
/*Deletion diagnostics in the cluster level*/

      id_Cluster = J(nobs,1,0);
      ClusterN = J(nobs,1,0);
      Levbetobs = J(nobs,1,0); /* observation leverage */ 
 
	  invBC1 = GINV(varKC); 
	  invBC1bet = GINV(varKC[1:p,1:p]);
	  invBC1alp = GINV(varKC[p+1:p+q,p+1:p+q]); 

      do i = 1 to nrow(n);
           X_c = X[firstx[i]:lastx[i],];
           y_c = y[firstx[i]:lastx[i]];
           lp_c = X_c*beta;
		   mu_c = linkid(lp_c);
		   *print i ;
		   
        if n[i] = 1 then do;
              B = CREATEL(mu_c, s2);
              C = CREATEC(X_c, mu_c);
              U_i= C`*inv(B)*(y_c - mu_c);
              Hi1=C*naive*C`*inv(B);
              U_c = C` *inv(B)* (y_c - mu_c)/(1 - Hi1);

             /* Halpha is set to 0 by default */
            %IF (%length(&CLSOUT) ^= 0) %THEN %DO; 
                 dfbetclsi = naive*U_c; 
                 Dfbetcls[i,] = dfbetclsi`;
                 
                 /* dfalpcls[i,] value is 0 by default */
                 dfalpclsi = dfalpcls[i,]`;
             %end;
                 

                 /*observation level cook's d*/
                %IF (%length(&OBSOUT) ^= 0) %THEN %DO; 
                 Dfbetobs[firstx[i],] = dfbetclsi`;/* dfbetobs = dfbetcls when n[i]=1*/
                 Dfobsboth = dfbetclsi//dfalpclsi;
                 cookdobs[firstx[i]] =   Dfobsboth`*invBC1*Dfobsboth/(p+q);
                 cookbetobs[firstx[i]] = Dfbetcls[i,]*invBC1bet*dfbetclsi/p;
                 cookalpobs[firstx[i]] = Dfalpcls[i,]*invBC1alp*dfalpclsi/q;
                 levbetobs[firstx[i]] = Hbeta[i];
                 dfalpobs[firstx[i],] = dfalpclsi`;/* dfalpobs = dfalpcls when n[i]=1*/
                %end;
                  /*CP level cook's d*/
                 %IF (%length(&CPOUT) ^= 0) %THEN %DO; 
				 dfbetCP[firstcp[i],] = dfbetclsi`;
				 dfalpCP[firstcp[i],] = dfalpclsi`;
	             clusterid_CP[firstcp[i]] =i;
	             periodid_CP[firstcp[i]] =periodid[firstx[i]];
                 cookbetcp[firstcp[i]] =   dfbetCP[firstcp[i],]*invBC1bet*dfbetCP[firstcp[i],]`/p;
                 cookalpcp[firstcp[i]] = dfalpCP[firstcp[i],]*invBC1alp*dfalpCP[firstcp[i],]`/q;
                 cookdcp[firstcp[i]] = Dfobsboth`*invBC1*Dfobsboth/(p+q);
                 CP_size[firstcp[i]] = 1;
                 %end;
                 
                 
                 
        end;  /* do loop, n[i]=1 */
        
        else do; /* do loop, n[i] ^=1 */

           U_i=J(p+q, 1, 0);
           U_c = J(p+q, 1, 0);
           Ustar_c = J(p+q, p+q, 0); 
           Z_c = Z[firstz[i]:lastz[i],];
           lp_alpha = Z_c * alpha;
           gamma_c = corrlinkid(lp_alpha);
		   E = deri2(Z_c, gamma_c);
* commands for beta;
           C = CREATEC(X_c, mu_c);
           B = CREATEB(mu_c, gamma_c, n[i],s2);
           L_c = CREATEL(mu_c, s2);
           A = CREATEA(mu_c, y_c);
           INVB=GINV(B);
           U_i[1:p]=C`*invB*A;
           
           CtinvB = C`*invB; 
           Hi1=C*naive*CtinvB;
           Hbeta[i] = trace(Hi1);  /* cluster leverage for beta */
*           IH=(I(n[i])-Hi1);
*           Ci=FINDINV(IH);

* commands for generalized inverse - beta;
           ai1 = invB;
           mm1 = C*sqe1;
           ai1A=ai1*A;
           ai1m1=ai1*mm1;
           run INVBIG(ai1A,ai1m1,mm1,A,1,p);/*calculate the (I - H^(0.5)H^(0.5))^(-1)*A*/
           U_c[1:p] = C`*ai1A;

* commands for alpha;
           VEE = J(NCHOOSE2(n[i]), 1, 0);
           R = J(NCHOOSE2(n[i]), 1, 0);
           DB = J(NCHOOSE2(n[i]), p, 0);
           RX = J(n[i], 1, 0);

                 %if %upcase(&alpadj) = MAEE %then %do;
                     SQVARFUN = sqrt(L_c);
                     INVSQVAR = inv(SQVARFUN);
*                IH=(I(n[i])-Hi1);
*                    Ci = FINDINV(IH);
* IH above needs to be symmetric, but we can use the following;
                     CT = t(C);
                     vminomega = B - C*naive*CT;
                     psd_vmin = is_pos_def(vminomega);
                     mineig = min(eigval(vminomega));
                     Ci = B*ginv(vminomega);
                     RX = INVSQVAR*(y_c - mu_c);
                     Gi=RX*RX`;

                if psd_vmin=0 then do;
	                     NPSDADJFLAG=1; 
	*                      print vminomega psd_vmin mineig;
	             end;   
                %end;/*MAEE*/

		  run VEERDB(VEE,R,DB,rangeflag,mu_c,gamma_c,X_c,y_c,1,n[i],p ,s2,Hi1,Ci,Gi,E);
 			     
           U_i[p+1:p+q]=E`*(R#(1/VEE));/*All Z_c need changes*/
           mm2 = E*sqe2;
           ai2R = R/VEE; 
           ai2m2 = mm2#(1/VEE);
           run INVBIG(ai2R,ai2m2,mm2,R,1,q);
           U_c[p+1:p+q] = E` *ai2R;

/* cluster leverage for alpha */                     
		  Zctvee = E`#(1/VEE`);
          Halpha[i] = trace(E*naivealp*Zctvee);  /* preserves memory */
* note trace(ABC) = trace(CAB), Searle (1982), p. 45;

          free Zctvee;            
          free Hi2;

        %IF %length(&CLSOUT)^= 0 %then %DO;
			/*begin GEECORR DIAG part*/
			/*cluster level diagnostics*/
              dfbetclsi = naive*U_c[1:p]; 
              dfalpclsi = naivealp*U_c[p+1:p+q];
              dfbetcls[i,] = dfbetclsi`;
              dfalpcls[i,] = dfalpclsi`;
              
	          Dfboth = dfbetclsi//dfalpclsi;
	          cookdcls[i] =   Dfboth`*invBC1*Dfboth/(p+q);
	          cookbetcls[i] = Dfbetcls[i,]*invBC1bet*dfbetclsi/p;
	          cookalpcls[i] = Dfalpcls[i,]*invBC1alp*dfalpclsi/q;
	
        %end;

    %IF %length(&OBSOUT)^= 0 %then %DO;
        /* cluster deletion diagnostics */
         run DFOBS(p,q,n, nobs, i, B,invB, C, naive, A,firstx,robust,VEE,Z_c,E,
                naivealp,R,id_cluster, clustern,invBC1,invBC1bet,invBC1alp,
                dfbetobs,cookbetobs,dfalpobs,cookalpobs,cookdobs);
         levbetobs[firstx[i]:lastx[i]] = vecdiag(Hi1);          
    %END;
       
    %IF %length(&CPOUT)^= 0 %then %DO;  
    *print firstz lastz;
    	run DFCP(p,q,n,i,B,invB,C,naive,A,firstcp[i] ,periodid[firstx[i]:lastx[i]],VEE,E,naivealp,R,clusterid_CP,periodid_CP,Zpair_periods[firstz[i]:lastz[i],],CP_size,
    	invBC1,invBC1bet,invBC1alp,dfbetCP,cookbetcp,Dfalpcp,cookalpcp,cookdcp);
    %END;

        end;  /*     do loop n[i] > 1  */
        end;     /*   do i = 1 to nrow(n) */




    /*for cluster level diagnostics dataset output, out of i iteration loop*/
     %IF (%length(&CLSOUT) ^= 0) %THEN %DO;
         id = (unique(id_cluster))`;
         if type(id)="N" then id=char(id);
         %if %sysevalf(%superq(zvar)^=,boolean) %then %do;
		         C1 = {"I" "NI" "Hbeta" "Halpha" "cookbetcls" "cookalpcls" &xvar &zvar};
		 %end;
	     %else %do;
	        if q=1 then C1 = {"I" "NI" "Hbeta" "Halpha" "cookbetcls" "cookalpcls" &xvar "ICC_parm1"};
			else if q=2 then C1 = {"I" "NI" "Hbeta" "Halpha" "cookbetcls" "cookalpcls" &xvar "ICC_parm1" "ICC_parm2"};;
		    else if q=3 then C1 = {"I" "NI" "Hbeta" "Halpha" "cookbetcls" "cookalpcls" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3"};;
			else  C1 = {"I" "NI" "Hbeta" "Halpha" "cookbetcls" "cookalpcls" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3" "ICC_parm4"};;
		 %end;
        
         xvar = dfbetcls;
         zvar = dfalpcls;
         i = (1:nrow(n))`;
         clout = i||n||Hbeta||Halpha||cookbetcls||cookalpcls ||xvar||zvar;  
         create &CLSOUT from clout [ rowname=id  colname=c1];
         append from clout [rowname=id] ;
     	 print "Cluster Deletion Diagnostics located in &CLSOUT";
      %END;

    %IF (%length(&OBSOUT) ^= 0) %THEN %DO; 
         obsID = (1:nobs)`;
         if type(obsID)="N" then obsID=char(obsID);
         xvar = dfbetobs;
         zvar = dfalpobs;
         clobsout = id_cluster||clustern||levbetobs||cookbetobs||cookalpobs||xvar||zvar;  
         *print clobsout;
          %if %sysevalf(%superq(zvar)^=,boolean) %then %do;
		         C2 = {"CLUSTER" "NI" "levbetobs" "cookbetobs" "cookalpobs" &xvar &zvar};
		 %end;
	     %else %do;
	        if q=1 then 
		    	C2 = {"CLUSTER" "NI" "levbetobs" "cookbetobs" "cookalpobs" &xvar "ICC_parm1"};
			else if q=2 then C2 = {"CLUSTER" "NI" "levbetobs" "cookbetobs" "cookalpobs" &xvar "ICC_parm1" "ICC_parm2" };
			else if q=2 then C2 = {"CLUSTER" "NI" "levbetobs" "cookbetobs" "cookalpobs" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3" };
			else q=4 then C2 = {"CLUSTER" "NI" "levbetobs" "cookbetobs" "cookalpobs" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3" "ICC_parm4"};

		 %end;
         create &OBSOUT from clobsout [rowname=obsID  colname=c2];
         append from clobsout [rowname=obsID ];
         print "Observation Deletion Diagnostics located in &OBSOUT";
    %END;


    %IF (%length(&CPOUT) ^= 0) %THEN %DO; 
         CPID = (1:lastcp[nrow(n)])`;
         if type(CPID)="N" then CPID=char(CPID);
         xvar = dfbetCP;
         zvar = dfalpCP;
         CPout = clusterid_CP||periodid_CP||CP_size||cookbetCP||cookalpCP||xvar||zvar;  

          %if %sysevalf(%superq(zvar)^=,boolean) %then %do;
		         C3 = {"CLUSTERID" "PERIODID" "ClusterPeriodSize" "cookbetCP" "cookalpCP" &xvar &zvar};%end;
	     %else %do;
	        if q=1 then 
		    	C3 = {"CLUSTERID" "PERIODID" "ClusterPeriodSize" "cookbetCP" "cookalpCP" &xvar "ICC_parm1"};
			else if q=2 then C3 = {"CLUSTERID" "PERIODID" "ClusterPeriodSize" "cookbetCP" "cookalpCP" &xvar "ICC_parm1" "ICC_parm2" };
		 	else if q=3 then C3 = {"CLUSTERID" "PERIODID" "ClusterPeriodSize" "cookbetCP" "cookalpCP" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3" };
			else if q=4 then C3 = {"CLUSTERID" "PERIODID" "ClusterPeriodSize" "cookbetCP" "cookalpCP" &xvar "ICC_parm1" "ICC_parm2" "ICC_parm3" "ICC_parm4"};
		%end;
         create &CPOUT from CPout [rowname=CPID  colname=c3];
         append from CPout [rowname=CPID ];
     print "Cluster Period Deletion Diagnostics located in &CPOUT";
    %END;
 end;   /*end GEECORR DIAG part*/

		/*
		print  Z;
         Z1 = {"Z1","Z2"};
         create Zdata from Z [colname=Z1  ];
         append from Z ;*/



          
 
       /******************************************/

     leavemvar:finish MAKEVAR_WT_DIAG;



/**************************************************************************************
MODULE: DFOBS
Calculates observation level deletion diagnostics for the mean and covariance parameters.

OUTPUT
deletion diagnostics for alpha and beta parameters, 
approximated by the one step iteration.  Location specified by user.

NOTES:
i indexes cluster
n[i]=cluster size
nrow(n) = number of clusters
p= number of mean parameters (beta)
q= number of covariance parameters (alpha)
naive = (D' Vinv D)inv = inxwx, for beta

B = covariance matrix for beta estimating equation, var(Y)

X_c = X covariates for cluster
obs_i = within cluster variable = vector of 1:n[i]

A = y-u
C = d u / d Beta = d u / d eta X = Linv X (=D from preisser notes)
B = V
W = Vinv

naivealp = (S'Tinv S)inv = M_2 inv
mu: Vector of n_i marginal means
n: Sample size (scalar) for cluster i
gamma: Vector of rho_jk between outcome j and k for cluster i
VEE:     working variances for alpha estimating equatinons = L` T L
R:   residual vector for alpha estimating equations = E 2
DB:  -E(d(corr)/dbeta) (i.e. Q) matrix 
E = d rho / d alpha = d rho / d eta Z = Linv Z (=E from preisser notes)
***********************************************************************************/
start DFOBS(p,q,n,nobs,i,B,invB,C,naivebet,A,firstx,robust,VEE,Z_c,E,naivealp,R,id_cluster,
clustern,invBC1,invBC1bet,invBC1alp,dfbetobs,cookbetobs,dfalpobs,cookalpobs,cookdobs);

            do ii = 1 to n[i];

            obs_i = 1:n[i];     *vector;
            comp = loc(obs_i^= ii);

                    /*for beta*/
                        /**************************************/

                    vii = B[comp, ii];
                    wii = invB[comp, ii];

                    Wti = invB[comp, comp]-(wii*wii`)/invB[ii, ii]; * equal to inv(B[comp, comp]) directly, 
                    inv(B[comp,comp]) = invB[comp,comp] - invB[comp,ii]*(invB[ii,ii])^-1*invB[ii,comp];
                    WV = wti*vii;
                    
                    Ctilda = C[ii, ]`-C[comp, ]`*WV; 
                    Htilda = Ctilda`*naivebet*Ctilda;   
                    Btilda = B[ii, ii]-vii`*WV;     /*V_ij in the perin's paper*/ 
                    Atilda = A[ii]-A[comp]`*WV;      

                    hdobs = Ctilda*Atilda/(Btilda-Htilda);

                    dfbetii = naivebet*hdobs;
                    dfbetobs[firstx[i]- 1 + ii, ] = dfbetii`;

                    /*Cluster information*/
                    id_cluster[firstx[i]- 1 + ii] =i;
                    clustern[firstx[i]- 1 + ii] =n[i];

                    /**************************************/

					/* Cook's D normalized by robust variance estimate  */
					cbetobs = (dfbetii`*invBC1bet*dfbetii)/p; 
					cookbetobs[firstx[i] -1 + ii] = cbetobs;


                    /*for alpha*//* if n[i]=2 special treatment is needed because there is only one row in z*/
                        /**************************************/
                    dimobs = n[i]-1;
            
                    if n[i] > 2 then do;
                        obs =  J(dimobs,1,0);
                        comp = J(NCHOOSE2(dimobs),1,0); end;
                    else if n[i] =2 then do; comp = .; obs = 1;end;                     
                
                    count=0;
                    ncomp=1;
                    nii=1;

/*routine to make list of observations to be removed*/ 
/*The correlation vector is defined with the order (1,2), (1,3),...*/                               
                    do obsj = 1 to n[i] - 1;
                        do obsk =obsj+1 to n[i];
                            count=count+1;
                                if ( ii^=obsj & ii^=obsk) then do;
                                    comp[ncomp]=count;
                                    ncomp=ncomp+1;
                                end;
                                else do;    obs[nii]=count;
                                        nii=nii+1;
                                end;
                        end;    end;/*end routine*/ 
                            

            VEEDIAG=diag(VEE);
            if n[i] > 2  then do;
                    veeDCOMPinv = diag(1/VEE[comp]);/*since. the VEE is diagnoal*/
                    veeveeinv = VEEDIAG[obs,comp]*veeDCOMPinv;  

                    Ztilda = E[obs,]-veeveeinv*E[comp,];
                    VEEtilda = VEEDIAG[obs,obs]-veeveeinv*VEEDIAG[comp,obs];
                    H2tilda = Ztilda*naivealp*Ztilda`;
                    Rtilda = R[obs]-veeveeinv*R[comp];
            end; /*  if n[i] > 2  */
            else if n[i] = 2 then do;
                    Ztilda = E[obs,];
                    VEEtilda = VEEDIAG[obs,obs];
                    H2tilda = Ztilda*naivealp*Ztilda`;
                    Rtilda = R[obs];
            end; /*  if n[i] = 2  */

                    ginvv = VEEtilda - H2tilda;
                    hadobs = Ztilda`*ginv(ginvv)*Rtilda;

                    dfalphaii = naivealp*hadobs; 
                    Dfalpobs[ firstx[i] -1 + ii, ] = Dfalphaii`;

        calpobs = (dfalphaii`*invBC1alp*dfalphaii)/q; 
        cookalpobs[firstx[i] -1 + ii] = calpobs;

        Dfboth = dfbetii//dfalphaii;
        cookdobs[firstx[i]-1+ii] =   Dfboth`*invBC1*Dfboth/(p+q);

end;

finish DFOBS;



/**************************************************************************************
MODULE: DFCP
Calculates cluster period level deletion diagnostics 
for the mean and covariance parameters.

OUTPUT
deletion diagnostics for alpha and beta parameters, 
approximated by the one step iteration.  Location specified by user.

Inputs:
i: indexes cluster
n[i]: cluster size
X_periodid_i: The period id vector in the X matrix in the ith cluster
periodid_CP: period number for each observation in the ith cluster
Zpair_periodid: period number for each pair-wise correlation in the ith cluster
output_start_point_i: The start points of the output vector;

p= number of mean parameters (beta)
q= number of covariance parameters (alpha)
naivebet = (D' Vinv D)inv = inxwx combing for all the clusters

B = covariance matrix for beta estimating equation, var(Y) from working correlation structure

X_c = X covariates for cluster
obs_i = within cluster variable = vector of 1:n[i]

In the ith cluster specific
A = y-u
C = d u / d Beta = d u / d eta X = Linv X (=D from preisser notes)
B = V
W = Vinv

naivealp = (S'Tinv S)inv = M_2 inv
mu: Vector of n_i marginal means
n: Sample size (scalar) for cluster i
gamma: Vector of rho_jk between outcome j and k for cluster i
VEE:     working variances for alpha estimating equatinons = L` T L
R:   residual vector for alpha estimating equations = E 2
DB:  -E(d(corr)/dbeta) (i.e. Q) matrix 
E = d rho / d alpha = d rho / d eta Z = Linv Z (=E from preisser notes)
dfbetCP,cookbetcp,Dfalpcp,cookalpcp,cookdcp, clusterid,periodid,CP_size: 
Diagnostics Outputs in the module
***********************************************************************************/
start DFCP(p,q,n,i,B,invB,C,naivebet,A,output_start_point_i ,X_periodid_i,VEE,E,naivealp,R,clusterid_CP,periodid_CP,Zpair_periodid,CP_size,
invBC1,invBC1bet,invBC1alp,dfbetCP,cookbetcp,Dfalpcp,cookalpcp,cookdcp);
/*Have a estimation of the number of the period number in the i th cluster*/

			unique_periodsi = unique(X_periodid_i);
			n_unique_periodsi = ncol(unique_periodsi);
			*print i unique_periodsi X_periodid_i;
            do loc = 1 to n_unique_periodsi;	/*loop through all periods in the cluster*/
				ij = loc(X_periodid_i = unique_periodsi[loc]);
	            comp = loc(X_periodid_i^= unique_periodsi[loc]);
	
	                    /*Diagnostics for beta*/
	                        /**************************************/
	
	                    vij = B[comp, ij]; 
	                    inv_vij = inv(B[comp, comp]); * equal to invB[comp,comp] - invB[comp,ij]*(invB[ij,ij])^-1*invB[ij,comp];
	                    WV = inv_vij*vij;
	                    *print WV;
	                    Ctilda = C[ij, ]`-C[comp, ]`*WV; /*transpose version of D_ijtilda mentioned ;*/
	                   *print Ctilda;
	                    Htilda = Ctilda`*naivebet*Ctilda;   /*Htilda need name changes)*/
	                   *print Htilda;
	                    Btilda = B[ij, ij]-vij`*WV;     /*V_ij in the perin's paper*/ 
	                   *print Btilda ;
	                    Atilda = A[ij]-WV`*A[comp];      /*uncorrected residuals*/
	                   *print Atilda;
	
	
	                    dfbetij = naivebet*Ctilda*ginv(Btilda-Htilda)*Atilda; /*the diagnostics of the beta of deleting cluster period ij to periods*/
	                    *print dfbetij;
	                    /*The location needs changes, i-1 +loc*/
	                    
	                    dfbetCP[output_start_point_i -1 + loc, ] = dfbetij`;
	
	                    /*Cluster period information*/
	                    clusterid_CP[output_start_point_i -1 + loc] =i;
	                    periodid_CP[output_start_point_i -1 + loc] =unique_periodsi[loc];
	                    CP_size[output_start_point_i -1 + loc] = sum(X_periodid_i=unique_periodsi[loc]);
	
	                    /**************************************/
	
						/* Cook's D normalized by robust variance estimate  */
						/* Do we need to discuss the Cook's D by different variance estimators?*/
						cbetcp = (dfbetij`*invBC1bet*dfbetij)/p; 
						cookbetcp[output_start_point_i -1 + loc] = cbetcp;
	
	                      /**************************************/
	                    /*for alpha*/
	                   /* if n[i]=2 special treatment is needed because there is only one row in z*/
	                  /**************************************/
	            
	            if n[i] > 2 then do;
	                    	cp = loc((Zpair_periodid[,1] = unique_periodsi[loc])|(Zpair_periodid[,2] = unique_periodsi[loc]));/*cp is the location for ij CP for the deletion of covariance matrix*/
	            end; /*  if n[i] > 2  */
	            else if n[i] =2 then do; cp = 1;
	            end;  /*1*1 matrix for the ICCs estimating equations*/                   
	           	
	           			VEEDIAG=diag(VEE);
	                    Etilda = E[cp,];
	                    VEEtilda = VEEDIAG[cp,cp];
	                    H2tilda = Etilda*naivealp*Etilda`;
	                    Rtilda = R[cp];

	
	                    ginvv = VEEtilda - H2tilda;
	                    hadcp = Etilda`*ginv(ginvv)*Rtilda;
	
	                    dfalphaij = naivealp*hadcp; 
	                    Dfalpcp[ output_start_point_i -1 + loc, ] = Dfalphaij`;
	
				        calpcp = (dfalphaij`*invBC1alp*dfalphaij)/q; 
				        cookalpcp[output_start_point_i -1 + loc] = calpcp;
	
	        Dfboth = dfbetij//dfalphaij;
	        cookdcp[output_start_point_i -1 + loc] =   Dfboth`*invBC1*Dfboth/(p+q);
end;/*for period do loop*/	


finish DFCP;



/***********************************************************

MODULE: FITPRENTICE
Performs Prentice (1988) Method

INPUT
      y: The outcome, 0/1 for binary outcome, -Inf to Inf for continuous outcomes and 0-Inf for count outcomes 
      X: Marginal mean covariates
      Z: Marginal correlation covariates
      n: Vector of cluster sample sizes
maxiter: Max number of iterations
epsilon: Tolerence for convergence

OUTPUT
    beta: A p x 1 vector of marginal mean parameter estimates  
   alpha: A q x 1 vector of marginal correlation parameter estimates  
  robust: Robust covariance matrix for beta and alpha
   naive: Naive (Model-Based) covariance matrix for beta
   niter: Number of iterations required for convergence
converge: Did the algorithm converge (0 = no, 1 = yes)
   
**************************************************************************************/

   start FITPRENTICE(beta, alpha, robust, naive, varMD, varKC,varFG, niter,
   converge, y, X, Z, n,s2, maxiter, epsilon, VEEFLAG, SINGFLAG, ROBFLAG,
   ALPFLAG, NPSDFLAG, NPSDADJFLAG,nobs,clusterid);
      p = ncol(X);
      q = ncol(Z);
      delta = J(p+q, 1, 2*epsilon);
      max_modi = 20;
      converge = 0;
     
      rangeflag=0;
      run INITVALS(beta, ustar, alpha, y, X, n, q,s2);

      do niter = 1 to maxiter while (max(abs(delta)) > epsilon);
         n_modi = 0;

         SINGFLAG=0;
         ALPFLAG=0;
         do until(n_modi > max_modi | rangeflag = 0);
            ustarold = ustar;

            NPSDFLAG=0;
            NPSDADJFLAG=0;
            run SCORE(U, UUtran, ustar, ustarold, beta, alpha, y, X, Z, n, p, q,s2, 0, 
                      rangeflag,VEEFLAG,NPSDFLAG,NPSDADJFLAG,clusterid);

            if VEEFLAG=1 then do;
                       print 'VEEFLAG=1 indicates program terminated due to division by zero in variance';
                       goto theend;
            end;
            
            if rangeflag = 1 then do;
                  %if %upcase(&shrink) = THETA %then %do;
                     if niter = 1 then alpha = J(q, 1, 0);
                     else do;
                        theta = theta - (0.5)**(n_modi+1) * delta;
                        beta = theta[1:p];
                        alpha = theta[p+1:p+q];
                     end;
                  %end;
                  %else %if %upcase(&shrink) = ALPHA %then %do;
                     if niter = 1 then alpha = J(q, 1, 0);
                     else do;
                        alpha = 0.95 * alpha;
                     end;
                  %end;
                  n_modi = n_modi + 1;
                  %if %upcase(&printrange) = YES %then %do; 
                     print "Iteration and Shrink Number", niter n_modi;
                  %end;
            end;  /* if rangeflag=1 then do */
           
         end; /* do until... */

         if n_modi > max_modi then do;
            %if %upcase(&printrange) = YES %then %do; 
               ALPFLAG=1;
               print "n_modi too great, more than 20 shrinks";
            %end;
               
*               goto theend;
         end;

         theta = beta // alpha;
         *print theta beta delta s2;
		 *print theta s2 ;	 
		 *print ustar u s2 niter; *all ustar to be non psd;
         psdustar = is_pos_def(ustar);
         mineig = min(eigval(ustar));
 *Allow the Ustar to be non psd;   
*         if is_pos_def(ustar)=1 then do;
*   if mineig < 1E-12 then          print ustar u psdustar mineig;
              delta = ginv(ustar)*U;
              *print delta;
              theta = theta + delta;
              beta = theta[1:p];
              alpha = theta[p+1:p+q];
              %if %upcase(&makephione) = YES %then %do; s2 = 1;%end;
              %else  %do;RUN PHI(Ustar, beta, alpha, y, X, Z, n, p, q,s2);%end;
              
              converge = (max(abs(delta)) <= epsilon); 

              *print theta beta delta s2;
*           end;
           if is_pos_def(ustar)=0 then do; 
              SINGFLAG=1;  
              *print ustar u psdustar mineig;
*              goto theend;
           end;   

      end;  /* do niter = 1 ... */
         
      niter = niter - 1;

      ustarold = ustar;
      
      if converge = 1 then run MAKEVAR_WT_DIAG(robust, naive, varMD, varKC,varFG, beta,
         alpha, ustarold, y, X, Z, n,  p, q,s2, VEEFLAG, ROBFLAG, NPSDADJFLAG, nobs,clusterid);
      if ROBFLAG = 1 then print "The variance estimators of parameters are not positve positive";
* ROBFLAG =1 indicates negative robust variance estimates;
      
      theend:finish FITPRENTICE;




/**************************************************************************************

MODULE: RESULTS
Creates printed output to screen of parameters and other information

INPUT
  beta: Vector of marginal mean parameters
 alpha: Vector of marginal correlation Parameters  
robust: Robust covariance matrix for beta and alpha
 naive: Naive (Model-Based) covariance matrix for beta 
 niter: Number of iterations until convergence
     n: Vector of cluster sample sizes

OUTPUT
To Screen

**************************************************************************************/

 start RESULTS(beta, alpha, robust, naive, varMD, varKC,varFG, niter, n,s2, CIC); 
      p = nrow(beta);
      q = nrow(alpha);
      printfinal = &printfinal;
 
      I=nrow(n);
 if   printfinal =1 then do;   
	 /* 2: Put variables into a table and use the TABLEPRINT statement. */
	Tbl1 = TableCreate('Characters','Type');
	call TableAddVar(Tbl1, {'Number of Clusters' 'Maximum Cluster Size' 'Minimum Cluster Size' 'Number of Iterations' },
	                        (nrow(n)) || (max(n))    || (min(n))  ||(niter )    );
	call TableSetVarFormat(Tbl1, {'Number of Clusters' 'Maximum Cluster Size' 'Minimum Cluster Size'  'Number of Iterations'},
	                            {'3.0' '3.0'   '3.0' '3.0'    });
	call TablePrint(Tbl1) ID='Characters' label='Summary of Clusters';
	
	Tbl2 = TableCreate('Characters',' Type');
	call TableAddVar(Tbl2, { 'Outcome Variable' 'Outcome Type' 'Outcome Link' 'Correlation Link' 'Adjustment'},
	                         {&yvar}   || {&ytype}    || {&link} || {&corrlink} ||{&alpadj} );
	call TablePrint(Tbl2) ID='Characters' label='GEEMAEE Specifications';

	criterion = {'CIC (smaller is better)'};
	values = CIC ; 
	Tb3 = TableCreate('Criterion', criterion);
	call TableAddVar(Tb3, {'Values'}, values);
	call TableSetVarFormat(Tb3, {'Values'}, {'D11.3'});
	call TablePrint(Tb3) ID='Criterion' label='Correlation Selection Criteria';
	
	 
	 
	      bsebc0 = sqrt(vecdiag(robust[1:p, 1:p]));
	      Tbeta_bc0=beta / bsebc0;
	      bpval0 = 2#(1 - probt(abs(Tbeta_bc0), I-p));
	
	      bsebc1 = sqrt(vecdiag(varKC[1:p, 1:p]));
	      Tbeta_bc1 = beta / bsebc1;
	      bpval1 =2#(1 - probt(abs(Tbeta_bc1), I-p));
	
	      bsebc2 = sqrt(vecdiag(varMD[1:p, 1:p]));
	      Tbeta_bc2 = beta / bsebc2;
	      bpval2 =2#(1 - probt(abs(Tbeta_bc2), I-p));
	      
	      bsebc3 = sqrt(vecdiag(varFG[1:p, 1:p]));
	      Tbeta_bc3 = beta / bsebc3;
	      bpval3 =2#(1 - probt(abs(Tbeta_bc3), I-p));
	    
	      stderr = sqrt(vecdiag(naive));
	      Tnai = beta / stderr;
	      pvalue =2#(1 - probt(abs(Tnai), I-p));
	
	Tbl3 = TableCreate('Variable', {&xvar}`);
	call TableAddVar(Tbl3, { 'Estimate' 'StdErr-MB' 'tValue-MB' 'P value MB' 'StdErr-BC0' 'tValue-BC0' 'P value BC0'
	'StdErr-BC1' 'tValue-BC1' 'P value BC1' 'StdErr-BC2' 'tValue-BC2' 'P value BC2' 'StdErr-BC3' 'tValue-BC3' 'P value BC3' },
	                        beta    || stderr    || Tnai     || pvalue ||bsebc0 || Tbeta_bc0 ||bpval0 ||
	                        bsebc1 || Tbeta_bc1 ||bpval1 ||bsebc2|| Tbeta_bc2 ||bpval2 ||bsebc3 || Tbeta_bc3 ||bpval3 );
	call TableSetVarFormat(Tbl3, { 'Estimate' 'StdErr-MB' 'tValue-MB' 'P value MB' 'StdErr-BC0' 'tValue-BC0' 'P value BC0'
	'StdErr-BC1' 'tValue-BC1' 'P value BC1' 'StdErr-BC2' 'tValue-BC2' 'P value BC2' 'StdErr-BC3' 'tValue-BC3' 'P value BC3' },
	                            { 'D11.3'   'D11.3'  '7.2'    'PVALUE6.4' 'D11.3'  '7.2'    'PVALUE6.4' 'D11.3'  '7.2'    'PVALUE6.4'
	                            'D11.3'  '7.2'    'PVALUE6.4' 'D11.3'  '7.2'    'PVALUE6.4' });
	call TablePrint(Tbl3) ID='Variable' label='Marginal Mean Parameter Estimates with Model-Based and Robust and Bias-corrected Standard Errors';
	
	Tblphi = TableCreate('Variable', 'PHI');
	call TableAddVar(Tblphi, { 'Estimate'  },
	                        s2 );
	call TableSetVarFormat(Tblphi, { 'Estimate'},
	                            { 'D11.3' });
	call TablePrint(Tblphi) ID='Variable' label='Dispersion Parameter Estimation';
	
	
	
	
	      asebc0 = sqrt(vecdiag(robust[p+1:p+q, p+1:p+q]));
	      Talpha_bc0 = alpha / asebc0;
	      apval0 = 2#(1 - probt(abs(Talpha_bc0), I-q));
	
	      asebc1 = sqrt(vecdiag(varKC[p+1:p+q, p+1:p+q]));
	      Talpha_bc1 = alpha / asebc1;
	      apval1 =2#(1 - probt(abs(Talpha_bc1), I-q));
	
	      asebc2 = sqrt(vecdiag(varMD[p+1:p+q, p+1:p+q]));
	      Talpha_bc2 = alpha / asebc2;
	      apval2 =2#(1 - probt(abs(Talpha_bc2), I-q));
	      
	      asebc3 = sqrt(vecdiag(varFG[p+1:p+q, p+1:p+q]));
	      Talpha_bc3 = alpha / asebc3;
	      apval3 =2#(1 - probt(abs(Talpha_bc3), I-q));
	      
	
	      %if &alpadj = UEE %then %do ;print 'SAS macro MAEE v2, 2022: Uncorrected Estimating Equations for Rho (Prentice, 1988)';   %end  ;    
	      %if &alpadj = SAEE %then %do; print 'SAS macro MAEE v2, 2022: Scalar-based Corrected Estimating Equations for Rho (Sharples and Breslow, 1992)'; %end  ;   
	      %if &alpadj = MAEE %then %do; print 'SAS macro MAEE v2, 2022: Matrix-adjusted Estimating Equations for Rho (Preisser, Lu and Qaqish, 2008)'; %end  ;   
	
	
	%if %sysevalf(%superq(zvar)^=,boolean) %then %do;
		Tbl4 = TableCreate('Variable', {&zvar}`);
	%end;
	%else %do;
	
	    if q=1 then zvar = {'Z1'};
		else if q=2 then zvar = {'Z1', 'Z2'};
		else if q=3 then zvar = {'Z1', 'Z2', 'Z3'};
		else if q=4 then zvar = {'Z1', 'Z2', 'Z3', 'Z4'};
		Tbl4 = TableCreate('Variable', zvar);
	%end;
	call TableAddVar(Tbl4, { 'Estimate' 'StdErr-BC0' 'tValue-BC0' 'P value BC0'
	'StdErr-BC1' 'tValue-BC1' 'P value BC1' 'StdErr-BC2' 'tValue-BC2' 'P value BC2' 'StdErr-BC3' 'tValue-BC3' 'P value BC3'},
	                        alpha     ||asebc0 || Talpha_bc0 ||apval0 ||
	                        asebc1 || Talpha_bc1 ||apval1 ||asebc2|| Talpha_bc2 ||apval2 ||asebc3 ||Talpha_bc3 ||apval3 );
	call TableSetVarFormat(Tbl4, { 'Estimate'  'StdErr-BC0' 'tValue-BC0' 'P value BC0'
	'StdErr-BC1' 'tValue-BC1' 'P value BC1' 'StdErr-BC2' 'tValue-BC2' 'P value BC2' 'StdErr-BC3' 'tValue-BC3' 'P value BC3'},
	                            { 'D11.3'   'D11.3'  '7.2'    'PVALUE6.4' 'D11.3'  '7.2'    'PVALUE6.4'
	                            'D11.3'  '7.2'    'PVALUE6.4' 'D11.3'  '7.2'    'PVALUE6.4'});
	call TablePrint(Tbl4) ID='Variable' label='Marginal Correlation Parameter Estimates with Robust and Bias-corrected Standard Errors';
end; /*print final*/	

/******************************************
write out estimated parameters into a SAS table
******************************************/
       %IF (%length(&ESTOUT) ^= 0) %THEN %DO; 
       *generate theta estimation output;
       	if q=1 then theta1 = {&XVAR "ICC_parm1"};
		else  if q=2 then theta1 = {&XVAR "ICC_parm1" "ICC_parm2"};
		else  if q=3 then theta1 = {&XVAR "ICC_parm1" "ICC_parm2" "ICC_parm3"};
		else  if q=4 then theta1 = {&XVAR "ICC_parm1" "ICC_parm2" "ICC_parm3" "ICC_parm4"};;
       theta = (beta//alpha)`;
       /*Generate the parameter estimations*/
       create &ESTOUT from theta [colname=theta1];
       append from theta ;
       close &ESTOUT;
       %end;
/******************************************
write out estimated covariance matrix into a SAS table
******************************************/      
       
       %IF (%length(&VAR_BETA) ^= 0) %THEN %DO; 
       *generate the Robust covariance estimation output;
       Covname = {&xvar};
       %IF %upcase(&VAR_BETA)=MB %THEN %DO; Varbeta = naive[1:p, 1:p]; %END;
       %IF %upcase(&VAR_BETA)=BC0 %THEN %DO; Varbeta = robust[1:p, 1:p]; %END;
       %IF %upcase(&VAR_BETA)=BC1 %THEN %DO; Varbeta = varKC[1:p, 1:p]; %END;
       %IF %upcase(&VAR_BETA)=BC2 %THEN %DO; Varbeta = varMD[1:p, 1:p]; %END;
       %IF %upcase(&VAR_BETA)=BC3 %THEN %DO; Varbeta = varFG[1:p, 1:p]; %END;
       create VAR_BETA from Varbeta [colname=Covname];
       append from Varbeta ;
       close VAR_BETA;
       %end;
 
       %IF (%length(&VAR_ALPHA) ^= 0) %THEN %DO; 
       *generate the Robust covariance estimation output;
	        if q=1 then alphaname = {"ICC_parm1"};
			else  if q=2 then alphaname = {"ICC_parm1" "ICC_parm2"};
			else  if q=3 then alphaname = {"ICC_parm1" "ICC_parm2" "ICC_parm3"};
			else  if q=4 then alphaname = {"ICC_parm1" "ICC_parm2" "ICC_parm3" "ICC_parm4"};;
       %IF %upcase(&VAR_ALPHA)=MB %THEN %DO; Varalpha = naive[p+1:p+q, p+1:p+q]; %END;
       %IF %upcase(&VAR_ALPHA)=BC0 %THEN %DO; Varalpha = robust[p+1:p+q, p+1:p+q]; %END;
       %IF %upcase(&VAR_ALPHA)=BC1 %THEN %DO; Varalpha = varKC[p+1:p+q, p+1:p+q]; %END;
       %IF %upcase(&VAR_ALPHA)=BC2 %THEN %DO; Varalpha = varMD[p+1:p+q, p+1:p+q]; %END;
       %IF %upcase(&VAR_ALPHA)=BC3 %THEN %DO; Varalpha = varFG[p+1:p+q, p+1:p+q]; %END;
       create VAR_ALPHA from Varalpha [colname=alphaname];
       append from Varalpha ;
       close VAR_ALPHA;
       %end;
           

finish RESULTS;

/**************************************************************************************/

* reasons for non-results are identified and tallied;
   VEEFLAG=0;
   SINGFLAG=0;
   ROBFLAG=0;
   ALPFLAG=0;
   NPSDFLAG=0;
   NPSDADJFLAG=0;
   Ypair_error=0;
   run GRABDATA(y, X, n, Z,  maxiter, epsilon ,nobs,clusterid,Ypair_error);  

   
  

   run FITPRENTICE(beta, alpha, robust, naive, varMD, varKC,varFG, niter, converge,
   y, X, Z, n,s2, maxiter, epsilon, VEEFLAG, SINGFLAG, ROBFLAG, ALPFLAG, NPSDFLAG, NPSDADJFLAG, nobs,clusterid);
 
   if Ypair_error = 1 then print  'THE ALGORITHM terminated due to incorrect Ypair specification';
   else if VEEFLAG = 1 then print  'THE ALGORITHM terminated due to nonpositive variance in weights';
   else if SINGFLAG = 1 then print  'THE ALGORITHM terminated due to singular MB
   covariance matrix';
   else if NPSDADJFLAG=1 then print 'THE ALGORITHM terminated due to singular finite sample adjusted matrix';
   else if converge = 0 then print 'THE ALGORITHM DID NOT CONVERGE';

   else if (converge=1 ) then do;
    run CIC_Criteria(CIC, n, y, X, beta, varKC,s2);
    run RESULTS(beta, alpha, robust, naive, varMD, varKC,varFG, niter, n,s2,CIC);
	end;


quit;




/*************************************************************************************/
%mend GEEMAEE;
/*************************************************************************************/
