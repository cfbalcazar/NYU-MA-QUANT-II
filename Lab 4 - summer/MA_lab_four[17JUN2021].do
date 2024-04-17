/*
=========================================================================================
Description: This is the do-file for the 4th lab of Quant II for MA students
in politics at NYU.

Authors: Felipe Balcazar
Ali T. Ahmed
Pedro L. Rodriguez 

Date begun: 02/16/16

Date last modified: 04/07/21

Purpose: Model specification
=========================================================================================
*/

clear all 
set more off , permanently 
capture log close 

* Set the path
loc dir=subinstr("`c(pwd)'","\do","",.)
capture log us


* Set the seed
set seed 12345


*=========================================================================================
/* Suppose we have a population and we want to estimate a given mean (ex. mean age)
 Ideally we would ask everyone their age and take the average, this unfortunately is not always possible
 Next best alternative: take a sample of this population, compute its average and use that as an estimate
 Estimator = formula we use to compute the average
 Estimate = the value produced by the formula for a given sample
 Bias is a characteristic of estimators not estimates
 Is this estimator biased?
 Formally: B(mu_hat) = E(mu_hat) - mu 
 Note: (E = expectation = average if we were to replicate the sampling process over and over).
 Intuition: an estimator is biased if the estimates it produces are systematically different from the population parameter.

 Write a program (allows us to repeat a set of commands easily)
*/
*=========================================================================================
set seed 12345 // Recall this allow us to replicate results when there're draws from random variables 
* =========================================================================================
* BIAS AND CONSISTENCY
* =========================================================================================
cap prog drop sim_simple
program sim_simple
		syntax [, nobs(integer 10) mean(real 0) sd(real 1) ]
        drop _all
		set obs `nobs'
		gen y = rnormal(`mean',`sd')
	    sum y
		gen meany=r(mean)
end

/*
The previous syntaxis shows a very easy way to construct a program in stata. The purpose
of this program is to simulate a normal distribution with a given mean, sd and observations.
*/

* Very imprecise estimator
sim_simple, nobs(100) mean(5) sd(10)
ttest y==5

* What happens if we draw many samples?
simulate meany, reps(10000) nodots: sim_simple, nobs(100) mean(5) sd(10)

* Now let's plot it
histogram _sim_1, normal xtitle(Mean Y) title(Sampling Distribution (Replications = 10000))

* More precise
sum _sim_1
ttest _sim_1==5


* What happens if we increase the number of observations?
foreach N in 10 100 1000 10000 {
	simulate meany, reps(10000) nodots: sim_simple, nobs(`N') mean(5) sd(10)

	* Now let's plot it
	histogram _sim_1, normal xtitle(Mean Y) title("Sampling Distribution (N = `N')") ///
	saving(g`N', replace)
}
graph combine g10.gph g100.gph, xcommon ycommon

graph combine g100.gph g10000.gph, xcommon ycommon

graph combine g10.gph g100.gph g1000.gph g10000.gph, xcommon ycommon


/* 
 Another characteristics of estimators is consistency
 Is our estimator consistent?
 Formally: plim(mu_hat) -> mu as n -> infinity
 Our estimator is consistent if the estimates it produces converge to the true population parameter as we collect more obs (increase N)
 How do you think the sampling distribution changes as we collect more observations? Increase the number of obs in the program?
*/


* Can we generate a biased estimator? Sure
cap prog drop sim_simple
program sim_simple
		syntax [, nobs(integer 10) mean(real 0) sd(real 1) ]
        drop _all
		set obs `nobs'
		gen y = rnormal(`mean',`sd')
		drop if y<0 // Our sample is no longer random!!!
	    sum y
		gen meany=r(mean)
end

* Things seem to fail now!
simulate meany, reps(10000) nodots: sim_simple, nobs(100) mean(5) sd(10)
sum _sim_1
ttest _sim_1==5

* Can a larger N help? Certainly not!
simulate meany, reps(10000) nodots: sim_simple, nobs(10000) mean(5) sd(10)
sum _sim_1
ttest _sim_1==5

*=========================================================================================
/*
In today's recitation we are going to focus on making inference about X_1. Even if we have
more variables these are not of our main interest. Thus X_2,... X_N will be controls.

Important to keep in mind:
- Define Y=f(X)+e. All elements in X is the set of controls.
- Let's focus on X_1 as our main independent regressor, and we are interested in knowing 
wheter X_1->Y (i.e., X_1 is related to Y).
- If X_2->X_1 and X_2->Y we will call the control X_2 a confounder.
- If X_2->Y but it is not the case that X_2->X_1, X_2 is a control but not a confounder.
- The set of confounders is contained within the set of all possible controls.
*/
*=========================================================================================
* Well-specified model
*=========================================================================================
clear all
program sim_regression_correct
        drop _all
        set obs 1000
		* generate independent variables
		gen x1 = abs(rnormal())
		* generate error term
		gen u = rnormal()
		* generate dependent variable
			gen y = x1 + u
		* run regression
		reg y x1
		* store coefficients
		scalar b_x1 = _b[x1]
		* store standard errors
		scalar se_x1 = _se[x1]
end

* sample once
sim_regression_correct
disp b_x1

* use simulate to run the program many times and store its output
simulate b_x1, reps(1000) nodots: sim_regression_correct
*simulate se_x1, reps(1000) nodots: sim_regression_correct

* take average of coefficient estimates
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Well-specified (Replications = 1000)) addplot(pci 0 1 8 1) legend(off) ///
saving(g_well_spec, replace)

*=========================================================================================
* Misspecified model? variable uncorrelated with x1
*=========================================================================================
clear all
program sim_regression_correct
	drop _all
	set obs 1000
	* generate independent variables
	gen x1 = abs(rnormal())
	gen x2 = abs(rnormal())*100
	* generate error term
	gen u = rnormal()
	* generate dependent variable
	gen y = x1 + x2 + u
	* run regression
	reg y x1 // Try ading x2 here!
	* store coefficients
	scalar b_x1 = _b[x1]
	* store standard errors
	scalar se_x1 = _se[x1]
end

* sample once
sim_regression_correct
disp b_x1

* use simulate to run the program many times and store its output
simulate b_x1, reps(1000) nodots: sim_regression_correct
*simulate se_x1, reps(1000) nodots: sim_regression_correct

* take average of coefficient estimates
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Omitted variable (Replications = 1000)) addplot(pci 0 1 0.15 1) legend(off) ///
saving(g_omitted_variable, replace)
/*
We omitted a control! What are the consequences? 
- The estimator will still be unbiased.
- The standard error will be less efficient (i.e., larger). Why? Because the more variables
 we include the SSE in this formula: SSE/[(SSx_1)*(1-R2_x1)]^(0.5), becomes smaller.
- Does that mean we need to include all variables we can in a model? No! Spurious correlations
between X_1, say, and other variables that are not part of the data generating processs for Y
can increase the standard error, so if you want to improve efficiency, you can include 
controls that are correlated with Y as long as they are completely uncorrelated to X_1. Of 
course, you can and should incluse controls that are confounders... we will see why below. 
*/

*=========================================================================================
* Misspecified model 1: omitted variable
*=========================================================================================
clear all
program sim_regression_correct
	drop _all
	set obs 1000
	* generate independent variables
	gen x2 = abs(rnormal())*100
	gen x1 = abs(rnormal()) + x2
	* generate error term
	gen u = rnormal()
	* generate dependent variable
			gen y = x1 + x2 + u
	* run regression
	reg y x1 // We assume is part of the error term E(u|X_2)!=0; try adding x2!
	* store coefficients
	scalar b_x1 = _b[x1]
	* store standard errors
	scalar se_x1 = _se[x1]
end

* sample once
sim_regression_correct
disp b_x1

* use simulate to run the program many times and store its output
simulate b_x1, reps(1000) nodots: sim_regression_correct
*simulate se_x1, reps(1000) nodots: sim_regression_correct

* take average of coefficient estimates
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Omitted variable bias (Replications = 1000)) addplot(pci 0 1 0.03 1) legend(off) ///
saving(g_omitted_var_bias, replace)


/*
We omitted a confounder! What are the consequences? 
- The estimator will be biased! This is a huge problem!
- Most of what practicioners do in most empirical papers is to address this issue.
- This means we need to control for confounders! Always!
*/

*=========================================================================================
* Misspecified model 2: irrelevant variable
*=========================================================================================
clear all
program sim_regression_correct
	drop _all
	set obs 1000
	* generate independent variables
	gen x1 = abs(rnormal())
	gen x2 = abs(rnormal())*100
	* generate error term
	gen u = rnormal()
	* generate dependent variable
	gen y = x1 + u
	* run regression
	reg y x1 x2
	* store coefficients
	scalar b_x1 = _b[x1]
	* store standard errors
	scalar se_x1 = _se[x1]
end

* sample once
sim_regression_correct
disp b_x1

* use simulate to run the program many times and store its output
simulate b_x1, reps(1000) nodots: sim_regression_correct
*simulate se_x1, reps(1000) nodots: sim_regression_correct

* take average of coefficient estimates
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Irrelevant variable (Replications = 1000)) addplot(pci 0 1 8 1) legend(off) ///
saving(g_irrelevant_variable, replace)


/*
What happens if we include an irrelevant variable:

- Recall: SSE/[(SSx_1)*(1-R2_x1)]^(0.5). Then if there are spurious correlations between the irrelevant 
variable and the outcome (X_1) or the independent variable (Y), we may poorly estimate the stanard
errors because SSE and R2_x1 will be affect as a consequence. Likely, the standard error will inflate. 

- This basically means that we need to be careful about what variables we choose as controls in
our regressions. We cannot just drop into the regression anything we feel like.


This is additional: In some cases, the estimator may be incorrect if the irrelevant variable 
is spuriously correlated to either Y and X_1 in this case. (If you ever get curious about 
why is this the case, you can use the Frisch-Waugh-Lovell theorem. )
*/


*=========================================================================================
* SUMMARY OF THE AFORMENTIONED
*=========================================================================================
* We will use our running example to go through the previous cases

drop _all
set obs 1000 
gen parent_schooling = 0+int((16-0+1)*runiform()) // Notice this is a confounder 
gen schooling = 0+int((16-0+1)*runiform()) ///
	+parent_schooling 
gen parent_income = rnormal(100,10) // Notice this is a control but not a confounder
gen randn=runiform()
egen gene_health = cut(randn), group(2) 
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 1.5
scalar beta_3 = 0.001
scalar beta_4 = -10
gen income = alpha + beta_1*schooling + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps

* If we do not include the confounder we will obtain a biased estimator
reg income schooling // The coefficient is 0.59! not 0.5
reg income schooling parent_schooling // By adding the confounder we solve the bias
reg income schooling parent_schooling parent_income gene_health // By adding controls we increase efficiency/precission (SE smaller)
* Lets add an irrelevant variable
gen irrelevant = int((16-0+1)*runiform()) 
reg income schooling irrelevant parent_schooling parent_income gene_health // The standard error became somewhat bigger

*=========================================================================================
* Misspecified model 3: non-linearity
*=========================================================================================
clear all
program sim_regression_correct
	drop _all
	set obs 1000
	* generate independent variables
	gen x1 = abs(rnormal())
	* generate error term
	gen u = rnormal()
	* generate dependent variable
			gen y = x1 - x1^2 + u
	* run regression
	reg y x1 c.x1#c.x1 // This command creates an interaction; similar to have a x1sq=x1^2
	* store coefficients
	scalar b_x1 = _b[x1]
	* store standard errors
	scalar se_x1 = _se[x1]
end

* sample once
sim_regression_correct
disp b_x1

* use simulate to run the program many times and store its output
simulate b_x1, reps(1000) nodots: sim_regression_correct
*simulate se_x1, reps(1000) nodots: sim_regression_correct

* take average of coefficient estimates
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Distribution (Replications = 1000)) addplot(pci 0 1 5 1) legend(off)

/*
A wise thing to do sometimes is to check the scatter plot of X_1 vs Y to determine whether
there are nonlinearites. 

This is additional: A more advanced, and correct, way to do this analysis is to use the 
Frisch-Waugh-Lovell theorem. That is to residualize Y and X_1 on all other covariates and 
then look at the scatter plot of residualized X_1 vs. residualized Y. See the end of this
do-file.
*/

*=========================================================================================
* LET'S TAKE A LOOK AT THIS IN OUR CLASSIC EXAMPLE
*=========================================================================================
drop _all
set obs 1000 
gen parent_schooling = 0+int((16-0+1)*runiform())  
gen schooling = 0+int((16-0+1)*runiform()) ///
+0.5*runiform()*parent_schooling 
gen parent_income = parent_schooling*100 + rnormal(10,10)
gen randn=runiform()
egen gene_health = cut(randn), group(2) 
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 0.2
scalar beta_3 = 0.001
scalar beta_4 = -10
gen income2 = alpha + beta_1*schooling-0.02*schooling^2 + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps

* Linear regression
reg income2 schooling  parent_schooling parent_income gene_health

* Quadratic fit
gen schooling_sq=schooling^2
reg income2 schooling schooling_sq parent_schooling parent_income gene_health
* Similarly
reg income2 schooling c.schooling#c.schooling parent_schooling parent_income gene_health

/*
We observe a worrisome misspecification problem here; not adding the square fo income 
schooling when the data generating process depends on it leads to very bad inference
about the effect of schooling on income2.
*/

*=========================================================================================
* Functional forms
*=========================================================================================

* logarithmic transformation (non-linear relationship)
clear all
* set number of obser
set obs 250
* generate independent variable
gen x1 = abs(rnormal())*10
* generate error term
gen u = rnormal()
* generate dependent variable
gen y = 10*log(x1) + u 
* plot relationship
twoway (scatter y x1) (lfit y x1)
* run regression
reg y x1
* look at residuals
predict e, resid
twoway (scatter e x1)
* take logs
gen logx1 = log(x1)
* plot
twoway (scatter y logx1) (lfit y logx1)
* re-run regression
reg y logx1
* look at residuals
predict e2, resid
twoway (scatter e2 logx1)

* quadratic relationship (non-linear relationship)
clear all
* set number of obser
set obs 250
* generate independent variable
gen x1 = abs(rnormal())*100
* generate error term
gen u = rnormal()
* generate dependent variable
gen y = 10*x1^2 + u
* plot relationship
twoway (scatter y x1) (lfit y x1)
* run regression
reg y x1
* look at residuals
predict e, resid
twoway (scatter e x1)
* include quadratic
gen x1sqr = x1^2  
* re-run regression
reg y x1 x1sqr
* look at residuals
predict e2, resid
twoway (scatter e2 x1)

*=========================================================================================
* RAMSEY RESET TEST
*=========================================================================================

* No non-linearities present
clear all
* set number of observations
set obs 250
* generate independent variable
gen x1 = abs(rnormal())*10
* generate error term
gen u = rnormal()
* generate dependent variable
gen y = 10*x1 + u
* run regression
reg y x1
* predict fitted values
predict yfit
* construct polynomials of fitted values
gen yfit2 = yfit^2
gen yfit3 = yfit^3
* re-run regression including polynomials of fitted values
reg y x1 yfit2 yfit3
* perform F-test
test yfit2 yfit3

* Non-linearities present
clear all
* set number of observations
set obs 250
* generate independent variable
gen x1 = abs(rnormal())*10
* generate error term
gen u = rnormal()
* generate dependent variable
gen y = 10*log(x1) + u
* run regression
reg y x1
* predict fitted values
predict yfit
* construct polynomials of fitted values
gen yfit2 = yfit^2
gen yfit3 = yfit^3
* re-run regression including polynomials of fitted values
reg y x1 yfit2 yfit3
* perform F-test
test yfit2 yfit3

/*
A reset test it tests whether non-linear combinations of the fitted values help 
explain the response variable. The intuition behind the test is that if non-linear
combinations of the explanatory variables have any power in explaining the 
response variable, the model is misspecified in the sense that the data generating 
process might be better approximated by a polynomial or another non-linear 
functional form.
*/
*=========================================================================================
* F-TESTS FOR NON-LINEAR SPECIFICATIONS (POLYNOMIALS)
*=========================================================================================
* No non-linearities present
clear all
* set number of observations
set obs 250
* generate independent variable
gen x1 = abs(rnormal())*10
* generate error term
gen u = rnormal()
* generate dependent variable
gen y = 10*x1 + x1^2 + 5*x1^3 + u
* generate polinomial terms
forvalues i=2/4 {
	gen x`i'=x1^`i'
}
* run regression and test compare the models using an F on the added terms
reg y x1 // We check our standard linear regression
test x1 // The coefficient is significant, but look at the coefficient!

* 2nd degree polynomial
reg y x1 x2 
test x2 // The coefficient is significant, but look at the coefficients!

* 3nd degree polynomial
reg y x1 x2 x3  
test x3 // The coefficient are significant, but the look right!

* 4th degree polynomial
reg y x1 x2 x3 x4   
test x4 // We reject accept the null that x4=0, so we do need 4th degree polynomial!

*=========================================================================================
* LET'S TAKE A LOOK AT THIS IN OUR CLASSIC EXAMPLE
*=========================================================================================
drop _all
set obs 1000 
gen parent_schooling = 0+int((16-0+1)*runiform())  
gen schooling = 0+int((16-0+1)*runiform()) ///
+0.5*runiform()*parent_schooling 
gen parent_income = parent_schooling*100 + rnormal(10,10)
gen randn=runiform()
egen gene_health = cut(randn), group(2) 
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 0.2
scalar beta_3 = 0.001
scalar beta_4 = -10
gen income2 = alpha + beta_1*schooling-0.02*schooling^2 + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps

* Ramsey test:
reg income2 schooling  parent_schooling parent_income gene_health
predict yfit
gen yfit2 = yfit^2
gen yfit3 = yfit^3
reg income2 schooling  parent_schooling parent_income gene_health yfit2 yfit3
test yfit2 yfit3 // It suggests that there is a polynomial

* Now let's use the F to decide the degree
forvalues i=2/4 {
	gen schooling`i'=schooling^`i'
}
* Notice that the second degree polynomial is preferred to just the linear model
reg income2 schooling schooling2 parent_schooling parent_income gene_health
test schooling2

* But the third and fourth degree polynomials are not preferred to the 2nd degree one  
reg income2 schooling schooling2 schooling3 parent_schooling parent_income gene_health
test schooling3
reg income2 schooling schooling2 schooling3 schooling4 parent_schooling parent_income gene_health
test schooling4

* =========================================================================================
* =========================================================================================
* THIS IS EXTRA AND ONLY FOR THOSE WHO ARE INTERESTED IN GOING BEYOND THE RECITATION
* =========================================================================================
* =========================================================================================
* The objective here is to obtain beta_1, that is the coefficient for scholing

drop _all
set obs 1000 
gen parent_schooling = 0+int((16-0+1)*runiform())  
gen schooling = 0+int((16-0+1)*runiform()) ///
+0.5*runiform()*parent_schooling 
gen parent_income = parent_schooling*100 + rnormal(10,10)
gen randn=runiform()
egen gene_health = cut(randn), group(2) 
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 0.2
scalar beta_3 = 0.001
scalar beta_4 = -10
gen income = alpha + beta_1*schooling + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps

* First we do the standard way
reg income schooling parent_schooling parent_income gene_health
scalar beta_1_hat=_b[schooling]

* Then we apply Frisch-Waugh-Lovell or, what is the same, FWL
reg income parent_schooling parent_income gene_health
predict res_Y, res
reg schooling parent_schooling parent_income gene_health
predict res_X, res

reg res_Y res_X
scalar beta_1_hat_FWL=_b[res_X]

disp  beta_1_hat-beta_1_hat_FWL // Awesome, right?!

* Why is this useful to investigate non-linear relationships visually?
gen income2 = alpha + beta_1*schooling-0.02*schooling^2 + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps
reg income2 schooling c.schooling#c.schooling parent_schooling parent_income gene_health

* Checking at the non-residualized scatter plot
twoway (scatter income2 schooling) (lfit income2 schooling), ///
title(Linear fit - not residualized) legend(off) saving(g_linear, replace) 

* Checking at the residualized scatter plot
reg income2 parent_schooling parent_income gene_health
predict res_Y2, res
reg schooling parent_schooling parent_income gene_health
predict res_X2, res

twoway (scatter res_Y2 res_X2) (qfit res_Y2 res_X2), ///
	title(Quadratic fit - residualized via FWL) ///
	saving(g_quadratic, replace) legend(off) // Huge difference!!!!

* Huge difference here!
graph combine g_linear.gph g_quadratic.gph, xcommon ycommon

* Finally, we conclude our regression should include the square of schooling
reg income2 schooling c.schooling#c.schooling parent_schooling parent_income gene_health



