/*
=========================================================================================
Description: This is the do-file for the lab review of Quant II for MA students
			in politics at NYU.

Authors: 	Felipe Balcazar
			Ali T. Ahmed
			Pedro L. Rodriguez 

Purpose: (1) Show how to use a do-file (2) Show basic commands.

Date begun: 02/16/16

Date last modified: 03/31/21

Purpose: Review of the previous labs
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

/*
What we are trying to do: 

- y = f(x) + u -> f(x) = systematic component, u = random error term (irreducible error).
- Want to estimate f(x).
- Aim: (1) Prediction: use observed values of X to predict Y (2) Inference: understand how Y is affected by X (preferably causally) (3) A combination of both.
- In general we want to: min E(Y - Yhat)^2 = min E(f(X) + u - fhat(X)]^2 = min [f(X) - fhat(X)]^2 + Var(u) = min [reducible error] + irreducible error.
- There are different methods to estimate f(x), Ordinary Least Squares (OLS) is one way. It is the one we use most of the time for inference.

Ordinary Least Squares (OLS):

- Advantages: (1) easy to interpret (2) reduces reality to a simple parametric form (linear) -> easy to estimate. 
- Disadvantages: (2) heavy reliance on assumptions. 
- Ass 1. f(X) is linear in X -> y ~ c + b_1*x_1 + b_2*x_2 + ... + u. u is the error term; X={x_1,x_2,...,x_n}.
- OLS is one way of estimatin f(X) assuming linearity (but it is not the only way).
- Intuition of OLS: minimize the distance between the observed value and the predicted value -> (y - yhat).
- OLS proposes: choose the line that minimizes the sum of the square of these distances.
- Why the square? We treat negative and positive differences equally.
- OLS makes certain assumptions (which if they do not hold then OLS will do a poor job of estimating f(X). We will later be more clear what "poor job" means).
- Ass 2. E(u) = 0 (not strictly necessary) -> holds if we include a constant in our regression.
- Ass 3. the error u is mean independent of the Xs -> E(u|x_1,x_2,...,x_n) = 0 (often referred to as the strict exogeneity assumption).
- Ass 2. + Ass 3. = zero conditional mean assumption. 
- Ass 4. The variance of the error term is independent of the Xs (homoskedastic) -> Var(u|X) = sigma.
- Ass 5. The errors are independent of each other -> E(u_i, u_j) = 0 for all i!=j (recall: ! = not equal to).
- Ass 6. Errors are distributed normally (not strictly necessary) -> u ~ N(0, sigma). For instance you can use bootstrap to compute the standard errors.
- Often we assume (y, x) are iid, that is (i)ndependent and (i)dentically distributed -> implies Ass 3, 4 and 5. 

Interpretation of coefficients:

Supose y=a+b*x+u
- b is the change in y when x increases in one unit: b*(x' - x) -> b*(101 - 100) = b

Supose log(y)=a+b*x+u
- b is the % change in y when x increases in one unit: b*(x' - x) -> exp(b*(101 - 100)) = exp(log(y') - log(y))=exp(log(y'))/exp(log(y)) 
														-> log(exp(log(y'))/exp(log(y)))= log(exp(log(y')))-log(exp(log(y)))=log(y')-log(y) ~ (y'-y)/y

Supose y=a+b*log(x)+u
- b is the change in y when x increases in one %: b*(log(x')-log(x)) = b*log(x'/x) -> b*log(101/100) ~ b*(101-100)/100 = b*1%

Supose log(y)=a+b*log(x)+u
- b is the % change in y when x increases in one %: b*(log(x')-log(x)) -> b*log(101/100) ~ b*1% = log(y') - log(y) ~ (y'-y)/y

*/
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
	
* Univariate regression
	reg income schooling
	
* Obtain predicted values and residuals
	predict yhat
	predict e, resid

* Do it manually 
	gen yhat_manual = _b[_cons] + _b[schooling]*schooling
	gen e_manual = income - yhat_manual
	gen diff_yhat=yhat-yhat_manual // They should be the same
	sum diff_yhat // Yes, they are!

* Checks
	corr yhat e                  // = 0 by Ass 3
	corr schooling e             // = 0 by Ass 3
	egen sume = sum(e)       
	sum sume                     // = 0 by Ass 2

* Averages
	sum schooling
	scalar mschooling = r(mean)
	sum income
	scalar mincome = r(mean)
	disp mincome-_b[_cons]-_b[schooling]*mschooling // Formally: E(Y) = E(c + b*X + u) = c + b*E(x)

*=========================================================================================
* What would happen if assumption 2 does not hold: assume that the constant is part of the error term so E(u)!=0
 
reg income schooling, nocons // Notice that the coefficient is biased!
predict e2, resid
sum e2
* In other words, when there's a constant and we do not include it we run into issues.

* Conversely, if there is no constant and we include it there should be no problem
gen income2 = beta_1*schooling + beta_2*ability + beta_3*parent_income + eps
reg income2 schooling
predict e3, resid
sum e3

*=========================================================================================
* What would happen if assumption 3 does not hold:
gen schooling2=schooling+ability+rnormal()
 
reg income schooling2 // Notice that the coefficient is biased!
predict yhat2
predict e4, resid       
corr schooling e4  // There is correlation between the regressor and the error term

*This is often referred to as omitted variable bias. Thus if we include the omitted variable (this is an extra)

reg income schooling2 ability // Voila!

*=========================================================================================
* What would happen if the errors are not homoskedastic
gen s = exp(0.16*schooling)
gen eps2 = s*eps
sum eps2
replace eps2=eps2-r(mean)

gen income3 = alpha + beta_1*schooling + beta_2*ability + beta_3*parent_income + eps2
reg income3 schooling 
predict yhat3
predict e5, resid
scatter e5 yhat3 

* We can fix this by using weighted least squares (this is an extra)
gen w_ones = 1
replace w_ones = w_ones/s
gen w_y=income3/s
gen w_schooling=schooling/s
reg w_y w_ones w_schooling, nocons 
drop yhat3 e5
predict yhat3
predict e5, resid
scatter e5 yhat3 

*=========================================================================================
* What would happen if the errors are correlated
gen n=_n
xtile ng=n, n(10)
gen eps3 = rnormal(ng,1) // We generate correlated errors
sum eps3
replace eps3=eps3-r(mean)

gen income4 = alpha + beta_1*schooling + beta_2*ability + beta_3*parent_income + eps3
reg income schooling
reg income4 schooling // Notice that the standard error is wrong!

* We can fix this with Huber-White standard errors (this is an extra).
reg income4 schooling, r


*=========================================================================================
* Multivariate linear regression

* First, notice what happen to the coefficients and the R squared when we add more variables
reg income schooling 
reg income schooling ability
reg income schooling ability parent_income

* Interpretation of coefficient: schooling. Income increase by one unit when schooling goes up by one unit
disp _b[schooling]*11-_b[schooling]*10

* Lets compute the t-statistic,  the p-value, confidence interval for 95% confidence (5% significance), and R squared

* t-stat
disp (_b[schooling] - 0)/_se[schooling]

* p-value (uses the model degrees of freedom and the t-stat
disp  ttail(e(df_r), (_b[schooling] - 0)/_se[schooling]) 
/* 
significance level: ex. it is significant at alpha%
alpha <- pvalue
interpretation: probability of rejecting the null hypothesis given that it is true (= type 1 error).
*/

* Confidence interval
disp "[`=round(_b[schooling]-1.96*_se[schooling],0.001)',`=round(_b[schooling]+1.96*_se[schooling],0.001)']"

* R-squared	
disp 1 - e(rss)/(e(mss)+e(rss))
disp e(mss)/(e(mss)+e(rss))

* F-test for model's significance

* F-statistic
disp e(mss)/e(rss)*e(df_r)/e(df_m) 

* Compute p-value for F-test
disp Ftail(e(df_m),e(df_r),e(mss)/e(rss)*e(df_r)/e(df_m))

