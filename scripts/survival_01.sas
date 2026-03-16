libname in "/home/indrajitsg0/work";

/* LOAD DATA */
data in.prost;
    infile '/home/indrajitsg0/raw/prostate_cancer.txt' dlm=' ' dsd missover firstobs=1;
    input patient treatment time status age sh size index;
run;

proc contents data=in.prost;
run;

proc print data=in.prost;
run;


/* Kaplan-Meier Method */
proc lifetest data=in.prost plots=survival(cl);
	time time*status(0);
run;

/* Generate separate survival curves for each strata*/
proc lifetest data=in.prost plots=survival(cl);
	time time*status(0);
	strata treatment;
run;

/* Perform Log Rank test */
proc lifetest data=in.prost plots=survival(cl);
	time time*status(0);
	strata treatment / test=all;
run;

/* Load motorette data */
data in.motorette;
    infile '/home/indrajitsg0/raw/motorette.csv' dlm=',' dsd missover firstobs=2;
    input x cens y;
run;

proc contents data=in.motorette varnum;
run;

proc print data=in.motorette;
run;

/* Create a scatterplot */
proc sgplot data=in.motorette;
    scatter x=x y=y / markerattrs=(symbol=circlefilled size=8);
    title "Motorette Scatter Plot";
run;

proc lifetest data=in.motorette plots=survival(cl);
	time y*cens(0);
run;

/* Create a binary variable based on a threshold of 180 */
data in.motorette;
	set in.motorette;
	format temp_binary $8.;
	temp_binary = "low";
	if x >= 180 then temp_binary = "high";
run;

proc lifetest data=in.motorette plots=survival(cl);
	time y*cens(0);
	strata temp_binary / test=logrank;
run;

/* Fit an accelerated failure time model  */
proc lifereg data=in.prost;
    model time * status(0) = treatment age sh size index/ dist=lnormal;
run;
/* AIC: 57.956 */

proc lifereg data=in.prost;
    model time * status(0) = treatment age sh size index/ dist=weibull;
run;
/* AIC: 52.956 */

proc lifereg data=in.prost;
    model time * status(0) = treatment age sh size index/ dist=llogistic;
run;

