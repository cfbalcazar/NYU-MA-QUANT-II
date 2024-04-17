/*
=========================================================================================
Description: This is the do-file for the 2nd lab of Quant II for MA students
in politics at NYU.

Authors: Felipe Balcazar
		 Ali T. Ahmed
		 Pedro L. Rodriguez 

Purpose: Run a regression,interpret the output and hypothesis testing.

Date begun: 02/16/16

Date last modified: 02/24/21

Notes: it is not necessarily expected of you to be able to run simulations, 
but hopefully it is pedagogically useful
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
* CONSTRUCT DATA
* =========================================================================================
* In this part we will construct our own data prior to exploring it.

* Define parameters
scalar alpha = 5                                  // defines a parameter rather than a variable
scalar beta = 0.5

* Generate data
set obs 1000                                      // fixes number of observations at 1000
set seed 021616                                   // allows us to replicate the draw by the pseudo-random number generator

** Independent variables
gen schooling = 0+int((16-0+1)*runiform())        // we draw the years of schooling from a uniform distribution on the interval [0,16]

** Error term
gen eps = rnormal()                               // we assume the error term is standard normal with mean 0 and std deviation of 1

** Dependent variables (the data generating process)
gen income = alpha + beta*schooling + eps         // generates our dependent variable

* Take a sub-sample of data
sum income
preserve
	sample 90
	sum income
restore


/*
Note: preserve and restore are useful insofar as the preserve the data that you loaded as a
temporary file that lasts the entire session, until you use restore. 

preserve/restore can be used in creative ways to merge data (we are not going to work on this).
*/


* =========================================================================================
* Clean data and transform if necessary

* Recode if necessary (not our case)
recode schooling 0/12=0 13/16=1, generate(highschool) label("High School")

* =========================================================================================
** Regression interpretation

* Univariate regression: using the dummy
reg income highschool

disp "High-school increases income by " _b[highschool]

twoway (scatter income highschool, sort jitter(3)) (lfit income highschool), ///
	title(Effect of Schooling on Income) ///
	ytitle(Income (thousands of dollars)) ///
	xtitle(Schooling)

/*
Notes: 
A. SS – These are the Sum of Squares associated with the three sources of variance, Total, Model and Residual. 
	SSTotal: The total variability around the mean.  Sum(Y – Ybar)^2. 
	SSResidual:  The sum of squared errors in prediction. Sum(Y – Yhat)^2. 
	SSModel: The improvement in prediction by using the predicted value of Y over just using the mean of Y. 
	Sum(Ypredicted – Ybar)^2. 

Another way to think of this is the SSModel is SSTotal – SSResidual. Note that the SSTotal = SSModel + SSResidual.  

Note that SSModel / SSTotal is equal to the value of R-Square. This is because R-Square is the proportion of
	the variance explained by the independent variables.

B. df – These are the degrees of freedom associated with the sources of variance.  
	The total variance has N-1 degrees of freedom, because one degree of freedom is lost to the constant. 
	So there are the total of 999 df.  One more df is lost for every indepdendent variable.
	Thus residuals degress of freedome are N-1-K where K is the number of regressors.

C. MS – These are the Mean Squares, the Sum of Squares divided by their respective DF.

D. Number of obs – This is the number of observations used in the regression analysis.

E. F and Prob > F – The F-value is the Mean Square Model divided by the Mean Square Residual 
	, yielding.  
	The p-value associated with this F value is very small (0.0000). 
	These values are used to answer the question 
	“Do the independent variables reliably predict the dependent variable?”. 

F. R-squared – R-Squared is the proportion of variance in the dependent variable (income) 
	which can be predicted from the independent variables (highschool). 

G. Adj R-squared – Adjusted R-square.  As predictors are added to the model, each predictor will 
	explain some of the variance in the dependent variable simply due to chance. 
	One could continue to add predictors to the model which would continue to improve the ability 
	of the predictors to explain the dependent variable, although some of this increase in R-square 
	would be simply due to chance variation in that particular sample.  
*/

* Now let's run the regression of interest
reg income schooling

* Now let's plot the fitted line
twoway (scatter income schooling, sort jitter(3)) (lfit income schooling), ///
	title(Effect of Schooling on Income) ///
	ytitle(Income (thousands of dollars)) ///
	xtitle(Schooling (years of schooling))
	
	
/*
_cons is the intercept on the y-axis. The coefficient for highschool is the slope of the line.
*/

* Access the stored results
ereturn list

* Total sum of squares
disp e(mss)  + e(rss)

* Degrees of freedom
disp e(df_m) // The number of regressors
disp e(df_r) - e(rank) // The number of observations minus the number of parameters 

* R-squared	
disp 1 - e(rss)/(e(mss)+e(rss))
disp e(mss)/(e(mss)+e(rss))

* F-test for model's significance

* F-statistic
disp e(mss)/e(rss)*e(df_r)/e(df_m) 

* Compute p-value for F-test
disp Ftail(e(df_m),e(df_r),e(mss)/e(rss)*e(df_r)/e(df_m))

* =========================================================================================
* Hypothesis testing
reg income schooling

/* Simple hypothesis:

H0: _b[schooling]=0

Ha: _b[schooling]!=0
*/

lincom schooling // Significance
testparm schooling // F-test

/* Yes we can use either F or t statistics. The reason is that for one single linear
	 restriciton the squared of the t-statistic equals the F statistic. 
	So we should obtain roughly the same p-values.*/ 



* Interpretation of on coefficient: schooling. Income increase by one unit when schooling goes up byu one unit
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

/*
For one-sided tests you can follow: 
https://www.stata.com/support/faqs/statistics/one-sided-tests-for-coefficients/
*/

* Let's add some complexity
gen noise=runiform() // this new variable should not be related to income

reg income schooling noise // Indeed noise is not significant; practice interpretation here

lincom schooling - noise// Equality of coefficients
lincom schooling + noise// Linear combination hypothesis

* =========================================================================================
* Now let's do postestimation

reg income schooling

* Obtain fitted values
predict yhat

* obtain residuals
predict e, resid

* check to see if y = yhat + resid
gen ypred = yhat + e
corr ypred income

* Explore the residuals
twoway (histogram e), title(Residuals histogram) ///
	ytitle(Income (thousands of dollars)) xtitle(Value of residuals)

* Let's check they are almost white noise
twoway (scatter  e yhat, sort jitter(3)) (lfit e yhat), title(Residuals vs. fitted) ///
	ytitle(Residuals) xtitle(Fitted values)

* Remember to close the log file
log close
