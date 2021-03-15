//The following code was used to calculate the true values.
// A large sample of 100,000,000 observations is generated following the data-generating mechanisms used to generate the exposure and confounder data of the simulation. Then, the 1-year and 5-year survival of each individual is obtained using the Weibull formula. The true values is then derived by obtaining the mean over the whole study population.

set seed 1805

forvalues cor=1/3 {

preserve
clear

set obs 100000000

local cor=`cor'

global hrz1 1.02
global hrz2 1.3
global hrz3 1.5
global lambda 0.2
global gamma 0.5
global hrdep 1.2
global maxt 5
global prevexp 0.5



tempname cmatrix
if `cor'==1 {
matrix `cmatrix' = (1, 0.8, 1, 0.8, 0.8, 1)
}
if `cor'==2 {
matrix `cmatrix' = (1, 0.5, 1, 0.5, 0.5, 1)
}
if `cor'==3 {
matrix `cmatrix' = (1, 0, 1, 0, 0, 1)
}


drawnorm z1 z2 z3,  means(60, 0, 0) sds(13, 1, 1) corr(`cmatrix') cstorage(lower) n(`obs')
replace z1=floor(z1)
replace z1 = 18 if z1<18
replace z1 = 99 if z1>99

foreach time in 1 5 {
	capture drop xb0 xb1
	gen xb0 = ln(${lambda}) + ln(${hrz1})*(z1-60) + ln(${hrz2})*z2 + ln(${hrz3})*z3
	gen xb1 = xb0 + ln(${hrdep})
	gen S`time'_dep0 = exp(-exp(xb0)*`time'^(${gamma})) 
	gen S`time'_dep1 = exp(-exp(xb1)*`time'^(${gamma})) 
}

foreach time in 1 5 {
	summ S`time'_dep0
	global truth_S0_T`time' `r(mean)'
	summ S`time'_dep1
	global truth_S1_T`time' `r(mean)'
	global truth_diff_T`time' = ${truth_S1_T`time'} - ${truth_S0_T`time'}
}
restore

macro dir
}


