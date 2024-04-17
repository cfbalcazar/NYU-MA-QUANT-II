/*
=========================================================================================
Description: This is the do-file for the 3rd lab of Quant II for MA students
			in politics at NYU.

Authors: 	Felipe Balcazar
			Ali T. Ahmed
			Pedro L. Rodriguez 

Purpose: (1) Show how to use a do-file (2) Show basic commands.

Date begun: 02/16/16

Date last modified: 03/03/21

Purpose: 	(1) Multiple regression
		(2) F-tests
		(3) the problem of collinearity
=========================================================================================
*/

clear all 						
set more off , permanently 
capture log close 

* Set the path
loc dir=subinstr("`c(pwd)'","\do","",.)			
capture log using "`dir'\Log_`c(current_date)'", replace text
cd "`dir'"

* =========================================================================================
* MULTIPLE REGRESSION
* =========================================================================================
clear
set obs 1000 
* Independent variables
	gen parent_schooling = 0+int((16-0+1)*runiform())  		// Parents' years of schooling
	
	gen schooling = 0+int((16-0+1)*runiform()) ///
		+0.1*parent_schooling 
		
	gen parent_income = parent_schooling*100 +rnormal(10,10)	// Parent's income
	
	gen randn=runiform()
	egen gene_health = cut(randn), group(2) 			// Individuals' genetic predisposition to illness 

* Error term
	gen eps = rnormal() 

* Coefficients

	scalar alpha = 5                           
	scalar beta_1 = 0.5
	scalar beta_2 = 0.2
	scalar beta_3 = 0.001
	scalar beta_4 = -10

** Dependent variables (the data generating process)
	gen income = alpha + beta_1*schooling + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps
	
* Let us compute the regression
	reg income schooling parent_schooling parent_income gene_health

* Do you notice something strange? Is there collinearity? Why?
	
log close


* =========================================================================================
* F TESTS
* =========================================================================================

/*
F-statistic: used in testing the null hypothesis that all of the model coefficients are 0 (combined significance test).
for a quick review of how to interpret the stata regression table see: https://bit.ly/2JhBLcy
*/

* Load the data
webuse census3, clear
 
* Run regression
 
regress brate medage medagesq i.region // i. indicates that the variable is a factor and creates dummies for each value of region
 
* Test coefficient on 3.region is 0
test 3.region=0

* Shorthand for the previous test command
test 3.region

* Estimate linear combination of coefficients

lincom 2.region + 3.region
	
* Two ways of running an F-test of joint significance (i.e., multiple restricitons)
/*
HO: 2.region=0
	3.region=0
Ha: i.region!= with i={2,3}
*/
testparm 2.region 3.region
test (2.region=0) (3.region=0)

* Test coefficient on 2.region=coefficient on 4.region
test 2.region=4.region
lincomest 2.region-4.region
/* Yes we can use either F or t statistics. The reason is that for one single linear
	 restriciton the squared of the t-statisct equals the F statistict. 
	So we should obtain roughly the same p-values.*/
	
	
	
	
* =========================================================================================
* Multicollinearity
* =========================================================================================
/*
The Gauss-Markov assumes the regressors being estimated aren't perfectly correlated with 
each other. This basically entails that the model is correctly specificed. Why?

This question is answered much better using linear algebra.

If the model is not correctly specified then the matrix of coefficients X does not have full rank.
In other words rank(X)!=n, where n is the number of regressors. This means that there are some variables
that are linear combinations of other; in other words the are linear dependencies in the X matrix. 

Why is this problematic?

We can write the OLS coefficient in matrix form as follows

B=(X'X)^(-1)X'Y

Here X' means X transposed, and (X'X)^(-1) is the inverse of X'X.

The problem with multicollinearity is that if there is perfect multicollinearity, X is not 
full rank, and thus X'X does not have an inverse. As a consequence we cannot estimate the coefficients B.
Stata drops one of the perfect multicolinear regressors by default to avoid this issue. 

If there's a high degree of multicollinearity, but it is not perfect, then we can estimate the coefficients.
Multicollinearity inflates the standard errors of the coefficients. For instance, if x_1 and x_2 are related,
there is less unique information we know about each with their respective relationships with Y. It becomes 
harder to understand how x_1 or x_2 relates to Y, so we have less certainty in the estimated values of beta_1 
or beta_2. This higher level of uncertainty can give us the wrong coefficients.

To see this he standard error for beta_1 is equivalent to SSE/[(SSx_1)*(1-R2_x1)]^(0.5), where SSE is the sum of squared errors,
SSx_1 is the sum of squared deviations for x_1 about its mean, and R2_x1 is the R2 for predicting x_1 with x_2 
(and any other variables in X if there are more than 2 involved) including an intercept. This tells us the proportion 
of variation in x_1 explained (shared) by x_2. Subtracting this from 1 gives us the unique (unshared) variation in x_1. 
Multicollinearity therefore makes the denominator of the standard error smaller, so the quotient gets bigger.

Another point to note is that 1/(1-R2_x1) is called the Variance Inflation Factor (VIF). It tells us the factor 
by which the variance of beta_1 is inflated due to X’s relationship with the other independent variables in 
the model. The square root VIF tells us how much the standard error of beta_1 (for instance) is inflated.
*/
* =========================================================================================
* We are going to use our example from lab 3, wherein we simulated income using education
* We will focus on the coeffcient for parent_schooling.

* First we create a simultion for a good model. 
clear
cap prog drop sim_regression_correct
program sim_regression_correct
	drop _all
	set obs 1000 
	gen parent_schooling = 0+int((16-0+1)*runiform())  
	gen schooling = 0+int((16-0+1)*runiform())
	gen parent_income = rnormal(10,10)
	gen randn=runiform()
	egen gene_health = cut(randn), group(2) 
	gen eps = rnormal() 
	scalar alpha = 5   
	scalar beta_1 = 0.5
	scalar beta_2 = 0.2
	scalar beta_3 = 0.001
	scalar beta_4 = -10
	gen income = alpha + beta_1*schooling + beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps
	reg income parent_schooling parent_income gene_health
	scalar b_x=_b[parent_schooling]
end

* use simulate to run the program many times and store its output
simulate b_x, reps(1000) nodots: sim_regression_correct
sum _sim_1

* plot distribution of output
histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Distribution (Replications = 1000))  ///
	saving(gwellspec, replace)
 
* Then we use our model with collinear regressors on parent schooling and parent income
clear
cap prog drop sim_regression_collinear
program sim_regression_collinear
	drop _all
	set obs 1000 
	gen parent_schooling = 0+int((16-0+1)*runiform())  
	gen schooling = 0+int((16-0+1)*runiform()) 
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
	reg income parent_schooling parent_income gene_health
	scalar b_x=_b[parent_schooling]
end

* use simulate to run the program many times and store its output
	simulate b_x, reps(1000) nodots: sim_regression_collinear
	sum _sim_1

* plot distribution of output
	histogram _sim_1, normal xtitle(OLS Coefficient of X1) title(Distribution (Replications = 1000))  ///
	saving(gcollinear, replace)

* Notice the albeit collinearity does not generate bias, it does generate much more uncertainty 
* about the parameter
graph combine gwellspec.gph gcollinear.gph

* =========================================================================================
* Diagnosis
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
reg income parent_schooling parent_income gene_health

* First we can take a look at the correlations
pwcorr schooling parent_schooling parent_income gene_health

/* 
Indeed schooling is correlated with parent income and schooling but the correlation
is not too high. However notice that the correlation between parent income and parent
schooling is quite high! This is a red flag! 
*/
* Indeed, in the regression
reg income schooling parent_schooling gene_health
* schooling was not facing issues.

* Using the Variance Inflation Factor (VIF)
reg income schooling parent_income parent_schooling gene_health
vif 

/* 
The VIF and tolerance (1/VIF) values for parent_schooling parent_income are worrisome.  
The very high VIF values indicate that these variables are possibly redundant.  
For example, after you know parent_schooling, you probably can predict parent_income very well.
In this example, multicollinearity arises because we have put in too many variables 
that measure the same thing, parent education.
*/

* Let's compute it by hand for parent_income
reg parent_income schooling parent_schooling gene_health
scalar R2x=e(r2)
disp 1/(1-R2x)


* What happens if we exclude parent_income?
reg income schooling parent_schooling gene_health
vif 

* This means that we can also check at the stability of parameters to adding and removing covariates.


* =========================================================================================
* BRIEF NOTES FOR HOME
* =========================================================================================

* =========================================================================================
* PROPERTIES OF LOG FUNCTION AND EXPONENTIAL FUNCTION

/*
In state you can log-transform a variable using log(x) or ln(x)

Recall as well:

log(a) if a<0 is not defined

log(0) is not defined

log(1) = 0

log(a*b) = log(a)+log(b)

log(a/b) = log(a)-log(b)

log(a^b) = b*log(a)

exp(0) = 1

exp(1) = 2.718 The value of the exponential constant

exp(a*b) = exp(a)^b for all b rational

exp(a+b) = exp(a)*exp(b)

exp(a-b) = exp(a)/exp(b)

log(exp(a)) = a = exp(log(a))
*/

disp log(100)

disp ln(100)

disp ln(100*50)-(ln(100)+ln(50))

disp ln(100/50)-(ln(100)-ln(50))

disp ln(100^5)-5*ln(100)

disp exp(0)

disp exp(1) // The value of the exponential constant

disp exp(10*(1/3))-(exp(10))^(1/3)

disp exp(10+5)-(exp(10)*exp(5))

disp exp(10-5)-(exp(10)/exp(5))

disp ln(exp(10)) - exp(ln(10))

disp ln(exp(10))


* =========================================================================================
* INTERPRETATION

/* 
If we want to log transform a variable (x) that contains zeroes, we can log transform
it as log(x+1). This is a very common adjustment. 

If somebody wants to dig really deep into this you can use the inverse hyperbolic sine
when you want to log-transform a variable. Since this can get tricky, I recommend it
just for those that would like to read econometric theory. 
*/

* For instance the schooling variable has zeroes
gen zeroes_schooling=(schooling==0) if schooling!=.
tab zeroes
* If we log transform without adding 1 those zeroes become missing
gen log_schooling_wzeroes=ln(schooling)
gen log_schooling_wozeroues=ln(schooling+1)
sum log_schooling*

* Let's use the simulated data
gen ln_schooling=ln(schooling+1)
gen ln_parent_schooling=ln(parent_schooling+1)

* Univariate regression: lin-lin regression
reg schooling parent_schooling

disp "One more year of parent schooling increases children schooling by " _b[parent_schooling]

* Univariate regression: log-lin regression
reg ln_schooling parent_schooling

disp "One more year of parent schooling increases children schooling by " _b[parent_schooling] "%"

* Univariate regression: lin-log regression 
reg schooling ln_parent_schooling

disp "One % increase in parent schooling increases children schooling by " _b[ln_parent_schooling]

* Univariate regression: log-log regression
reg ln_schooling ln_parent_schooling

disp "One % increase in parent schooling increases children schooling by " _b[ln_parent_schooling] "%"

log close




 
