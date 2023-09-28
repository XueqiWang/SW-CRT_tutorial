

* The updated SAS codes for the Tutorial paper with covclusters intercept periods effects;
* Aim: To see whether there is an increasing trend of the intervention effect estimation under three intervention effects model;


FILENAME REFFILE9 '/home/u41116673/Tutorial Paper/sample_data_PCGS_tutorial_paper_ICC01.csv';
%include '/home/u41116673/sasuser.v94/MAEEV2.02.sas';
libname dat "/home/u41116673/CH";
*read in data for Connect Home study;
PROC IMPORT DATAFILE=REFFILE9
	DBMS=csv
	OUT=CH_data;
	GETNAMES=YES;
RUN;

proc sort data =CH_data;
by clusters periods;
run;


* For the caregiver preparedness score;
data pcg_comp_CH;
     set CH_data;
     where missing(y)=0;
run;


proc sort data =pcg_comp_CH;
by clusters periods;
run;

* generate the extended intervention effect model;
data pcg_comp_CH_2;
   set pcg_comp_CH;
   			if clusters=1  then do;
   			           if periods <=5 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-7)/10;
						   ex_incre_trt = min((periods-7)/5,1);
						   end; * The first site has intervention from the 8th periods;

   			end;
   			
   			else if clusters=2 then do;
   			           if periods <=7 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-9)/10;
						   ex_incre_trt = min((periods-9)/5,1);
						   end; * The first site has intervention from the 10th periods;

   			end;
   			
   			else if clusters=3 then do;
   			           if periods <=9 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-11)/10;
						   ex_incre_trt = min((periods-11)/5,1);
						   end; * The first site has intervention from the 12th periods;

   			end;
   			
   			else if clusters=4 then do;
   			           if periods <=11 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-13)/10;
						   ex_incre_trt = min((periods-13)/5,1);
						   end; * The first site has intervention from the 14th periods;

   			end;
   			
   			else if clusters=5 then do;
   			           if periods <=13 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-15)/10;
						   ex_incre_trt = min((periods-15)/5,1);
						   end; * The first site has intervention from the 16th periods;

   			end;
   			
   			else if clusters=6 then do;
   			           if periods <=15 then 
   			           		do;
   			           		incre_trt = intervention;
   			           		ex_incre_trt = intervention;
   			           		end;
   			           else do;
   			               incre_trt = (periods-17)/10;
						   ex_incre_trt = min((periods-17)/5,1);
						   end; * The first site has intervention from the 18th periods;

   			end;
   			
   	pcgs = pcgs_integer	;
   	period = periods-1; * The real period used in the marginal mean model;
 run;


proc sort data = pcg_comp_CH_2;
  by  clusters period;
 run;
 



title "The GEE/MAEE analysis of caregivers' preparedness with linear periods effect and average intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar=  int period intervention_binary, clusterID= clusters, 
periodID = period, corr  = NE , corrlink = Fishersz, makevone=NO ,makephione=NO, alpadj=MAEE,maxiter=300, epsilon=0.00001 );

 title "The GEE/MAEE analysis of caregivers' preparedness with linear periods effect and incremental intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar=  int period intervention, clusterID= clusters, 
periodID = period, corr  = NE , corrlink = Fishersz,makevone=NO ,makephione=NO, alpadj=MAEE,maxiter=100, epsilon=0.00001 );
 
 title "The GEE/MAEE analysis of caregivers' preparedness with linear periods effect and extended incremental intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar= int period ex_incre_trt, clusterID= clusters, 
periodID = period, corr  = NE ,corrlink = Fishersz, makevone=NO ,makephione=NO, alpadj=MAEE,maxiter=100, epsilon=0.00001 );


 title "The GEE/UEE analysis of caregivers' preparedness with linear periods effect and average intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar=  int period intervention_binary, clusterID= clusters, 
periodID = period, corr  = NE , corrlink = Fishersz, makevone=NO ,makephione=NO, alpadj=UEE,maxiter=100, epsilon=0.00001 );

 title "The GEE/UEE analysis of caregivers' preparedness with linear periods effect and incremental intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar=  int period incre_trt, clusterID= clusters, 
periodID = period, corr  = NE , corrlink = Fishersz, makevone=NO ,makephione=NO, alpadj=UEE,maxiter=100, epsilon=0.00001 );
 
 title "The GEE/UEE analysis of caregivers' preparedness with linear periods effect and extended incremental intervention effect under NE correlation structure";
%GEEMAEE(xydata= pcg_comp_CH_2,  yvar= pcgs ,ytype = normal, link =identity, xvar= int period ex_incre_trt, clusterID= clusters, 
periodID = period, corr  = NE , corrlink = Fishersz, makevone=NO ,makephione=NO, alpadj=UEE,maxiter=100, epsilon=0.00001 );

