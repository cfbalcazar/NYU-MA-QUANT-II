/*
=========================================================================================
Description: This is the do-file for the 5th lab of Quant II for MA students
			in politics at NYU.

Authors: 	Felipe Balcazar
			Ali T. Ahmed
			Pedro L. Rodriguez 

Date begun: 02/16/16

Date last modified: 04/14/21

Purpose: Heteroskedasticity, dummies and interactions effects
=========================================================================================
*/

clear all 						
set more off , permanently 
capture log close 

* Set the path
loc dir=subinstr("`c(pwd)'","\do","",.)			
capture log us


	
*=========================================================================================
*=========================================================================================
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
* MEASUREMENT ERROR
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
reg income_cer schooling_cer merror

/*
The big problem with non-classical measurement error is that we do not observe it, we
do not know what shape it takes! 

To address this problem we use often instrumental variables!!!
*/


*=========================================================================================
* DUMMY VARIABLES
*=========================================================================================

/*
Dummy variables are those variables that only take two values: 1 or 0. Usually we thinik of
them as follows: 
- The dummy variable X_1 evaluates whether a conditions is met. If the condition is met it 
takes the value of 1, zero otherwise.
- An example is whether you have high-school education or not, thus for instance X_1=1 if
the observation has high-school diploma or zero if she doesn't. 
*/


*=========================================================================================
* Dummy variables in our running example
*=========================================================================================
drop _all
set obs 1000 
gen parent_schooling = 0+int((16-0+1)*runiform()) // Notice this is a confounder 
gen schooling = 0+int((16-0+1)*runiform()) ///
	+parent_schooling 
gen parent_income = rnormal(100,10) // Notice this is a control but not a confounder
gen randn=runiform()
* Notice that gene_health is a dummy variable, we say that it takes the value of 1 if poor health
egen gene_health = cut(randn), group(2)

* Let's define our data generating process
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 1.5
scalar beta_3 = 0.001
scalar beta_4 = -100 // I inflated this term
gen income = alpha + beta_1*schooling+ beta_2*parent_schooling + beta_3*parent_income + beta_4*gene_health + eps

* Checking at the scatter plot, the two clouds are a result of the gene_health dummy
* therefore it is wise to control for it to reduce the standard error
twoway (scatter income schooling) (lfit income schooling), ///
	title(Linear fit) legend(off)
	
* How do we know this?! Let's use a bit of FWL to demonstrate this
reg income gene_health
predict res_Y2, res
twoway (scatter res_Y2 schooling) (qfit res_Y2  schooling), ///
	title(Linear fit with residualized income) legend(off) 

* We can also generate a dummy variable here for schooling, assuming that having a high school diploma
* is equivalent to exhibiting at least 6 years of schooling (this may happend if we cannot measure years)
gen high_school=(schooling>12) if  schooling!=. // Adding the if statement here is very important in stata 

* In this regression the coefficient for high_school is interpreted as the effect of meeting the condition
* that is, the effect of having high-school education
reg income high_school parent_schooling parent_income gene_health

* Let's assume now a more complicated case wherein we have a categorical variable
gen degree_schooling=""
replace degree_schooling="No education" if schooling==0
replace degree_schooling="Primary" if schooling>0
replace degree_schooling="Secondary" if schooling>6
replace degree_schooling="Tertiary" if schooling>12

* Let's convert this to a factor variable for using it in the regression
encode degree_schooling, gen(dsch)
tab dsch
tab dsch, nolab // Notice it transformed it to numbers

* Now let's run the regression, it is important to use the prefix i. to tell 
* stata this is a factor (or a collection of dummies)
reg income i.dsch parent_schooling parent_income gene_health // the base category is no education

*=========================================================================================
* INTERACTION EFFECTS; MODERATORS; HETEROGENOUS EFFECTS
*=========================================================================================
/*
An interaction effect explores the heterogenous effect of a regressor, say X_1, on our
outcome Y. The variable that is used for the interaction, say X_2, is called a moderator.
Therefore when people use the terms interaction effects, moderators or heterogenous effects 
they are basically referring to the same thing!

We say that X_2 is a moderator because it moderates the effect of X_1 on Y. Importantly X_2
must always occur before X_1, that is: X_2->X_1 and not the other way around!

In a regression setting this looks like

	Y=alpha+beta_1*X_1+beta_2*X_2+beta_3*X_1*X_2+e

beta_3 is called the interaction effect. X_1 and X_2 are called the component terms of
the interaction X_1*X_2, and they should always be included.

If X_2 is a dummy variable (the easiest way to get this), then we have the following (assuming 
X_1 is also a dummy):

beta_1=E(Y|X_1=1,X_2=0)-E(Y|X_1=0,X_2=0)

beta_1 + beta_2=E(Y|X_1=1,X_2=1)-E(Y|X_1=0,X_2=1) This is called the full interaction effect and 
	usually what we interpret when we talk about an interaction effect.

beta_2 = E(Y|X_1=1,X_2=1)-E(Y|X_1=0,X_2=1)-[E(Y|X_1=1,X_2=0)-E(Y|X_1=0,X_2=0)] 
This is the difference in the effect of X_1 between the groups in your sample when they take 
different values for X_2.
 
Interaction effects are tricky at first. I have seen people getting the interpretation
wrong. Therefore keep always in mind the difference between beta_2 (the interaction effect)
and beta_1+beta_2 (the full interaction effect).
*/
clear all

* Let's simplify our running example to address this issue
drop _all
set obs 1000 
gen schooling = 0+int((16-0+1)*runiform())
* Let's also create a variable for gender, which will act as our moderator
gen randn=uniform()
egen gender = cut(randn), group(2)

* Let's define our data generating process
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 1.5
scalar beta_3 = 1.2 // I'm making this quite big to make everything visible

* Notice that in this data generating process gender is a moderator!
gen income = alpha + beta_1*schooling+ beta_2*gender + beta_3*gender*schooling + eps

* explore graphically
* plot without interaction
	twoway (scatter income schooling) (lfit income schooling), legend(off)
	
* plot accounting for interaction (using jitter for a nicer plot)
	twoway (scatter income schooling if gender==0, jitter(15)) ///
	(lfit income schooling if gender==0) ///
	(scatter income schooling if gender==1, jitter(15)) ///
	(lfit income schooling if gender==1), ///
	legend(on order(2 "gender = 0" 4 "gender = 1"))

* regression analysis
* only schooling
	reg income schooling
* schooling and gender without interaction
	reg income schooling gender 
* schooling and gender with interaction
	gen inter = schooling*gender
	reg income schooling gender inter
	
* The effect of one unit increas on schooling given different values of gender:
	
	lincom schooling // This is the effect of schooling conditional on gender=0
	
	lincom schooling + inter // This is the effect of schooling conditional on gender=1	
	
* we cannot also run this regression as follows, without using inter
	reg income schooling gender c.schooling#i.gender // c is telling stata we deal with a continuous variable
	reg income c.schooling##i.gender // ## is telling stata we deal with a full interaction

* The hypothesis testing is very similar herein
		
	lincom schooling // This is the effect of schooling conditional on gender=0
	
	lincom schooling + c.schooling#1.gender // This is the effect of schooling conditional on gender=1	
	
	
* Let's use margins to plot this
	reg income c.schooling##i.gender
	margins, dydx(schooling) at(gender=(0(1)1)) post
	marginsplot, recast(scatter) yline(0) xscale(range(-.5 1.5)) ///
	title("Interaction effect") ytitle("Effect of schooling") ///
	xtitle("Gender")


	
	
*=========================================================================================
* CONTINUOUS MODERATORS (TO TAKE HOME)
*=========================================================================================
/*
With continous moderators, interactions effects are interpreted in the same way, however,
given that we have a continuous variable, we may want to vary the value of the moderator
X_2 accross all its values (it's support).
*/
clear all

* Let's simplify our running example to address this issue
drop _all
set obs 1000 
gen schooling = 0+int((16-0+1)*runiform())
* Let's also create a variable for gender, which will act as our moderator
gen randn=uniform()
gen health_index = randn*100

* Let's define our data generating process
gen eps = rnormal() 
scalar alpha = 5                           
scalar beta_1 = 0.5
scalar beta_2 = 1.5
scalar beta_3 = 0.001 

* Notice that in this data generating process gender is a moderator!
gen income = alpha + beta_1*schooling+ beta_2*health_index + beta_3*health_index*schooling + eps

* I'm going to use this method because is the one I find most useful for coding properly
reg income c.schooling##c.health_index // Notice i. changed to c. here because the moderator is continuous

* But what value should we use to analize the full effect?
lincom schooling + c.schooling#c.health_index // 1? But it goes until 99!
lincom schooling + 50*c.schooling#c.health_index // 50? But it goes until 99!
lincom schooling + 70*c.schooling#c.health_index // 70? But it goes until 99!
 
* Let's use a margins plot to address this conundrum!
reg income c.schooling##c.health_index 
margins, dydx(schooling) at(health_index=(0(1)100)) post cont
marginsplot, recast(line) plot1opts(lcolor(gs0)) ///
	ciopts(recast(rarea) color(gs9)) graphregion(color(white)) ///
	title("Full interaction effect") ///
	ytitle("Effect of schooling on income") ///
	xtitle("Value of the helath index") level(95) ///
	legend(off) 
	
	

