/*
=========================================================================================
Description: This is the do-file for the 7th lab of Quant II for MA students
			in politics at NYU.

Authors: 	Felipe Balcazar
			Ali T. Ahmed
			Pedro L. Rodriguez 

Date begun: 02/16/16

Date last modified: 04/14/21

Purpose: Dummies and interactions effects
=========================================================================================
*/

clear all 						
set more off , permanently 
capture log close 

* Set the path
loc dir=subinstr("`c(pwd)'","\do","",.)			
capture log us

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
