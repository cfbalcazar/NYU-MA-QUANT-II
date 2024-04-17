/*
=========================================================================================
Description: This is the do-file for the 8th lab of Quant II for MA students
			in politics at NYU.

Authors: 	Felipe Balcazar
			Ali T. Ahmed
			Pedro L. Rodriguez 

Date begun: 02/16/16

Date last modified: 04/14/21

Purpose: 
=========================================================================================
*/

clear all 						
set more off , permanently 
capture log close 

* Set the path
loc dir=subinstr("`c(pwd)'","\do","",.)			
capture log us

/*
Consider our usual regression equation:

Y=a+beta_1*X+e

If the variability of the random disturbance (e) is different across elements of X, 
if Var(e)!=constant, we say then that error term is heteroskedastic.
e might be heteroskedastic if for instance var(e)=g(X), where g(.) is some function.  


The presence of heteroskedasticity implies the absence of homoskedasticity, violating
the assumptions for the Gauss-Markov theorem. This means that OLS estimators are not
BLUE because their variance is not the lowest of all other unbiased estimators. 

Heteroscedasticity does not cause OLS estimates to be biased, although it can cause 
biased standard errors. Therefore inferences obtained from data analysis are 
suspect---and likely wrong. 
*/

*=========================================================================================
* WHAT WOULD HAPPEND IF ERRORS ARE NOT HOMOSKEDASTIC?
*=========================================================================================
clear
set obs 1000 
* Independent variables
	gen schooling = 0+int((16-0+1)*runiform()) 	// Individuals' years of schooling
		
	gen parent_income = rnormal(10,10)	// Parent's income
	
	gen ability=  rnormal(10,10)

* Error term
	gen eps = rnormal() 

* Coefficients

	scalar alpha = 5                           
	scalar beta_1 = 0.5
	scalar beta_2 = 0.2
	scalar beta_3 = 0.001

* Dependent variables (the data generating process)
	gen income = alpha + beta_1*schooling + beta_2*ability + beta_3*parent_income + eps
	
	
* Let us introduce heteroskedasticity to the data generating process 
	gen s = exp(0.16*schooling)
	gen eps_het = s*eps
	sum eps_het
	replace eps_het=eps_het-r(mean) // Centering at zero

	gen income_het = alpha + beta_1*schooling + beta_2*ability + beta_3*parent_income + eps_het
	
*=========================================================================================	
* First let us establish a baseline
 reg income schooling
 predict y_hat
 predict resid, resid
  
* Then we run the regression with the heteroskedasticity issue
 reg income_het schooling 
 predict y_hat_het
 predict resid_het, resid
 
 /*
 Note that the SE are bigger; precission falls and therefore our estimate is larger.
 This is a result of specification error.
 */

 
 * Let's take a look visually
 scatter resid y_hat, title(Homoskedastic errors) saving(hom.gph, replace)
 scatter resid_het y_hat_het, , title(Heteroskedastic errors) saving(het.gph, replace)

 graph combine hom.gph het.gph

* We can also carry out diagnostic tests: in both cases we reject the null of constant variance
 reg income_het schooling 
 estat imtest // White’s test 
	
 estat hettest // Breusch-Pagan test


*=========================================================================================
* WEIGHTED LEAST SQUARES
*=========================================================================================
/* 
We incorporate knowledge of the variance of observations into the regression. That is, 
we weight the observations by their sigma value.

Think of a variable (Z) and suppose that sigma_i=lambda*Z_i, where lamabda is the 
scaling factor. If we divide the model by Z:

Y/Z=a*1/Z+beta_1*X/Z+u/Z

Y_i/Z_i

then it follows

E((u/Z)^2)=lambda^2

thus the variance term is constant again. Furthermore we dont need to know lambda.
*/

gen w_ones = 1
replace w_ones = w_ones/s
gen w_y=income_het/s
gen w_schooling=schooling/s
reg w_y w_ones w_schooling, nocons 
predict yhat_wls
predict resid_wls, resid
scatter resid_wls yhat_wls 
	
/* 
What happens if we don't know the structure of the error term? In general we don't kow this! 

Sometimes Z can be a population count, sometimes X itself, sometimes the error
in a two-stage regression.In practice, we use our best guess.
*/

gen wls_ones = 1/schooling
gen y_wls=income_het/schooling
reg y_wls wls_ones // Herein the constant is the coefficient for schooling!!!

* Notice, however, that this did not do a good job, why? schooling is not a good proxy for s.


* Using a two stage regression we may overcome the previous problem
reg income_het schooling 
predict y_hat_het
predict resid_het, resid
gen abs_resid_het=abs(resid_het)
reg  abs_resid_het y_hat_het 
predict wls
pwcorr wls s // The correlations is very high

* Let's test this
reg income_het schooling [w=1/wls^2] // It seems to work

* Alternatives include (wls0 performs weighted least squares)
ssc install wls0

wls0  income_het schooling, wv(schooling) type(abse)


*=========================================================================================
* HUBER-WHITE STANDARD ERRORS
*=========================================================================================
/*
Heteroscedasticity-consistent (HC) standard errors. In regressions analysis these are known
as Huber-White standard errors. When the type of heteroskedasticity is unknown, that it is
it's nature HC standard errors provide consistent estimates of the standard error in large
samples but they may be biased in small samples!
*/

reg income_het schooling, robust // This can be done by simpling adding the robust option. 

*=========================================================================================
* MESUREMENT ERROR
*=========================================================================================
/*
Measurement error in the dependent variable: Random measurement error in the depedenent
variable means that X_1 explains less variation of Y. This in turn means that the SSE is larger
and therefore the standard error of beta_1 will be less efficient, and thus bigger.

Measurement error in the independent variable: Random measurement error in the indepedenent
variable means that the measurement error in X_1, for instance, becomes part of the error 
term in the regression equation thus creating an endogeneity bias. OLS estimation of beta_1
will lead the coefficient to be biased towards zero. This bias is therefore
called attenuation bias.

The previous two types of measurement error are called classical measurement error.

Measurement error as a source of omitted variable bias: If there is a source for measurement error
and this element affects both the dependent and independent variable, then it becomes a confounder.
And therefore it can induce omitted variable bias.
*/

* First, let's check what happens when we introduce randome measruemen error
gen schooling_er=schooling+rnormal(10,1)
gen income_er=income+rnormal(10,1)

* This estimate is our benchmark
reg income schooling

* Measureent error in X, leads to attenuation bias. Beta is closer to zero.
reg income schooling_er 

* Measureent error in Y leads to bigger standard errors!.
reg income_er schooling

* Uncorrelated measurement error in Y and X leads to attenuation bias with larger SE
reg income_er schooling_er

* What if the measurement error is correlated?
gen merror= rnormal(10,1)
gen schooling_cer=schooling+merror
gen income_cer=alpha + beta_1*schooling + beta_2*ability + beta_3*parent_income + 5*merror+eps

* Note that correlated measurement error leads to bias, it is like omitted vairiable bias
* because merror->schooling_cer and merror->income_cer
reg income_cer schooling_cer 

* Therefore once we control for it we solve the problem
reg income_er schooling_er merror

/*
The big problem with non-classical measurement error is that we do not observe it, we
do not know what shape it takes! 

To address this problem we use often instrumental variables!!!
*/



	
