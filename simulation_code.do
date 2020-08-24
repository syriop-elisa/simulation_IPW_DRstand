



local nsim 1000



//////////////////////////////////////////	
//Define program to run one repetition
/////////////////////////////////////////	
capture program drop CIsim
program define CIsim, rclass
syntax, [ COR(integer 1) ///	  
		]

local obs 2000
local orz1 1.01 
local orz2 1.3 
local orz3 1.4 
local hrz1 1.02
local hrz2 1.3
local hrz3 1.5
local lambda 0.2
local gamma 0.5
local hrdep 1.2
local maxt 5
local prevexp 0.5


local df 3
local dftvc 2
	



local time 1


//covariance matrix
tempname Cmat
if "`cor'"=="1" {
matrix `Cmat' = (1, 0.8, 1, 0.8, 0.8, 1)
}
if "`cor'"=="2" {
matrix `Cmat' = (1, 0.5, 1, 0.5, 0.5, 1)
}
if "`cor'"=="3" {
matrix `Cmat' = (1, 0, 1, 0, 0, 1)
}



// generate covariates
capture drop z1 z2 z3
drawnorm z1 z2 z3,  means(60, 0, 0) sds(13, 1, 1) corr(`Cmat') cstorage(lower) n(`obs')
replace z1=floor(z1)
replace z1 = 18 if z1<18
replace z1 = 99 if z1>99



// generate treatment
capture drop logoddsdep dep
di "`=log(`orz1')'*z1 +  `=log(`orz2')'*z2 + `=log(`orz3')'*z3"
gen logoddsdep = logit(`prevexp') + `=log(`orz1')'*(z1-60) +  `=log(`orz2')'*z2 + `=log(`orz3')'*z3
gen dep = rbinomial(1,invlogit(logoddsdep))
tab dep

capture drop z1c
gen z1c=z1-60

// Simulate Relative Survival 
capture drop deadcanc
survsim timerel deadcanc, dist(weibull) lambda(`lambda') gamma(`gamma') ///
	cov(dep `=log(`hrdep')' z1c `=log(`hrz1')' z2 `=log(`hrz2')' z3 `=log(`hrz3')') maxt(`maxt')



//Simulate the other-cause time
capture drop agediag
rename z1 agediag
capture drop yydx sex age year  yeardiag dx
gen yydx=2009
gen sex=1
gen yeardiag=2009
gen age=agediag
replace age=99 if age>=100
gen year=yeardiag
gen dx = mdy(1, 1, 2009)
capture drop patid
gen patid=_n

forvalues i=0/4 {
               capture drop _merge survprob rate
               replace age=min(agediag+`i',99)
               replace age=99 if age>=100
               replace year=yydx+`i'
               sort sex year age dep
               quietly merge m:1 sex year age dep using popmort, keep(matched master)
			   sort patid
               gen timeback`i'=(-log(runiform()))/(-ln(survprob))
               replace timeback`i'=1 if timeback`i'>=1
}
               
gen timeback=timeback0 if timeback0<1

forvalues i=1/4 {
               quietly replace timeback=timeback`i'+`i' if timeback`i'<1 & timeback==.
}

replace timeback=5 if timeback==. 

drop timeback? 
 


//Create the all-cause time and stset 
capture drop t
gen t=min(timeback,timerel)	


capture drop d id
gen d=1 if timeback>=timerel
replace d = 2 if timeback<timerel
replace d= 0 if t>=5
drop  year age 
drop _merge 
drop rate
gen id=_n
capture drop timeback timerel
stset t, failure(d=1,2) id(id)



capture drop tt
//gen tt=`time' in 1
range tt  1 5 2

//perform analysis

gen age = min(int(agediag + _t),99)
gen year = int(yydx + _t)

quietly merge m:1 sex year age dep using popmort, keep(matched master)

gen z1=agediag 
capture drop tt
//gen tt=`time' in 1
range tt  1 5 2



///////// Regression standardisation
foreach covs in "z1 z2 z3" "z1 z2" "z1"  ""   {
	capture drop _at* _contrast*
	stpm2 dep `covs', df(`df') scale(hazard)  bhaz(rate)
	standsurv, at1(dep 0) at2(dep 1) contrast(difference) timevar(tt) se  transform(none)
	
	
	forvalues t = 1/2 {
		return scalar S_RA`model'_dep0_T`t'= _at1[`t']
		return scalar S_RA`model'_dep1_T`t' = _at2[`t']
		return scalar S_RA`model'_dep0_T`t'_se = _at1_se[`t']
		return scalar S_RA`model'_dep1_T`t'_se = _at2_se[`t']
		return scalar DIFF_RA`model'_T`t' = _contrast2_1[`t']
		return scalar DIFF_RA`model'_T`t'_se = _contrast2_1_se[`t']
	}
	capture drop _at* _contrast*
	standsurv, at1(dep 0) at2(dep 1) contrast(difference) timevar(tt) se mestimation transform(none)
	forvalues t = 1/2 {
		return scalar S_RA`model'_dep0_T`t'_m = _at1[`t']
		return scalar S_RA`model'_dep1_T`t'_m = _at2[`t']
		return scalar S_RA`model'_dep0_T`t'_m_se = _at1_se[`t']
		return scalar S_RA`model'_dep1_T`t'_m_se = _at2_se[`t']
		return scalar DIFF_RA`model'_T`t'_m = _contrast2_1[`t']
		return scalar DIFF_RA`model'_T`t'_m_se = _contrast2_1_se[`t']
	}
	
	
	local model=`model'+1
	
}



////////// IPW

local model 1
foreach covs in "z1 z2 z3" "z1 z2" "z1"  ""  {
    
	capture drop _at* _contrast*
	capture drop wt pr 
	capture drop wt_logit

	logit dep `covs'
	predict pr
	summ dep
	local Pdep `r(mean)'
	gen wt_logit = cond(dep==1,`Pdep',1-`Pdep')/cond(dep==1,pr,1-pr)
	
	
	mrsprep using popmort.dta, 	pmother(sex dep) ///
							agediag(agediag) ///
							datediag(dx) ///
							pmmaxyear(2009) pmyear(year) pmage(age) ///
							breaks(0(0.2)5) ///
							indweights(wt_logit) ///
							newframe(,replace) /// 
							by(dep)
							
	

	
						   
    stset tstop [iw=wt], enter(tstart) failure(event==1)	
	capture drop tt
	
	range tt  1 5 2					


	
	stpm2 dep, df(`df') scale(hazard)  bhaz(meanhazard_wt) tvc(dep) dftvc(3) vce(cluster id) 
	
	capture drop _at* _contrast*
	standsurv if _n==1, at1(dep 0) at2(dep 1) contrast(difference) timevar(tt) se  transform(none)
	forvalues t = 1/2 {
		return scalar S_IPW`model'_dep0_T`t' = _at1[`t']
		return scalar S_IPW`model'_dep1_T`t' = _at2[`t']
		return scalar S_IPW`model'_dep0_T`t'_se = _at1_se[`t']
		return scalar S_IPW`model'_dep1_T`t'_se = _at2_se[`t']
		return scalar DIFF_IPW`model'_T`t' = _contrast2_1[`t']
		return scalar DIFF_IPW`model'_T`t'_se = _contrast2_1_se[`t']
	}
	local model=`model'+1
	frame change default
}



////////// Doubly Robust 
// True Model RA - WRONG IPW
capture drop tt
range tt  1 5 2
local model 1
foreach covs in "z1 z2 z3" "z1 z2" "z1"  ""  {
	capture drop _at* _contrast*
	capture drop wt pr
	logit dep `covs'
	predict pr
	gen wt = cond(dep==1,`Pdep',1-`Pdep')/cond(dep==1,pr,1-pr)
	//stset t d [iw=wt]
	stset t [iw=wt], failure(d=1,2)   
	stpm2 dep z1 z2 z3, df(`df') scale(hazard)    bhaz(rate) vce(cluster id) 
	standsurv , at1(dep 0) at2(dep 1) contrast(difference) timevar(tt) se 
	forvalues t = 1/2 {
		return scalar S_DRa`model'_dep0_T`t' = _at1[`t']
		return scalar S_DRa`model'_dep1_T`t' = _at2[`t']
		return scalar S_DRa`model'_dep0_T`t'_se = _at1_se[`t']
		return scalar S_DRa`model'_dep1_T`t'_se = _at2_se[`t']
		return scalar DIFF_DRa`model'_T`t' = _contrast2_1[`t']
		return scalar DIFF_DRa`model'_T`t'_se = _contrast2_1_se[`t']
	}	
	
	local model=`model'+1
	
}

///////// Doubly Robust 
// Wrong Model RA - True IPW
local model 1
foreach covs in "z1 z2 z3" "z1 z2" "z1"  ""  {
	capture drop _at* _contrast*
	capture drop wt pr
	logit dep z1 z2 z3
	predict pr
	gen wt = cond(dep==1,`Pdep',1-`Pdep')/cond(dep==1,pr,1-pr)
	//stset t d [iw=wt]
	stset t [iw=wt], failure(d=1,2)   
	stpm2 dep `covs', df(`df') scale(hazard)  bhaz(rate) vce(cluster id) 
	standsurv, at1(dep 0) at2(dep 1) contrast(difference) timevar(tt) se 
	
	forvalues t = 1/2 {
		return scalar S_DRb`model'_dep0_T`t' = _at1[`t']
		return scalar S_DRb`model'_dep1_T`t' = _at2[`t']
		return scalar S_DRb`model'_dep0_T`t'_se = _at1_se[`t']
		return scalar S_DRb`model'_dep1_T`t'_se = _at2_se[`t']
		return scalar DIFF_DRb`model'_T`t' = _contrast2_1[`t']
		return scalar DIFF_DRb`model'_T`t'_se = _contrast2_1_se[`t']
	}	
	
	local model=`model'+1
	
}


	
ereturn clear
end

tempname estimates states
postfile `estimates' i method cor covs tvc time  unexp unexpse exp expse dif difse using RightEstimates.dta, replace
postfile `states' i cor  str2000 s1 str2000 s2 str1100 s3 using RightStates.dta, replace

quietly{
	set seed 1805
noi _dots 0, title("Simulation running...")
profiler on
	forval i = 1/`nsim' {
			forval cor = 1/3 {
					
				post `states' (`i') (`cor')  (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
				CIsim, cor(`cor')
				//local t 1
				local tvc 0
				
				forvalues covs=1/4 {
				 forvalues t=1/2{
					post `estimates' (`i') (1) (`cor') (`covs') (`tvc') (`t') (r(S_RA`covs'_dep0_T`t'))  (r(S_RA`covs'_dep0_T`t'_se)) (r(S_RA`covs'_dep1_T`t'))   (r(S_RA`covs'_dep1_T`t'_se)) (r(DIFF_RA`covs'_T`t'))  (r(DIFF_RA`covs'_T`t'_se))
					post `estimates' (`i') (2) (`cor') (`covs') (`tvc') (`t') (r(S_RA`covs'_dep0_T`t'_m)) (r(S_RA`covs'_dep0_T`t'_m_se))  (r(S_RA`covs'_dep1_T`t'_m))  (r(S_RA`covs'_dep1_T`t'_m_se)) (r(DIFF_RA`covs'_T`t'_m))  (r(DIFF_RA`covs'_T`t'_m_se))
					post `estimates' (`i') (3) (`cor') (`covs') (`tvc') (`t') (r(S_IPW`covs'_dep0_T`t')) (r(S_IPW`covs'_dep0_T`t'_se)) (r(S_IPW`covs'_dep1_T`t'))    (r(S_IPW`covs'_dep1_T`t'_se)) (r(DIFF_IPW`covs'_T`t'))  (r(DIFF_IPW`covs'_T`t'_se))
					post `estimates' (`i') (4) (`cor') (`covs') (`tvc') (`t') (r(S_DRa`covs'_dep0_T`t'))  (r(S_DRa`covs'_dep0_T`t'_se))  (r(S_DRa`covs'_dep1_T`t'))   (r(S_DRa`covs'_dep1_T`t'_se)) (r(DIFF_DRa`covs'_T`t'))  (r(DIFF_DRa`covs'_T`t'_se))
					post `estimates' (`i') (5) (`cor') (`covs') (`tvc') (`t') (r(S_DRb`covs'_dep0_T`t')) (r(S_DRb`covs'_dep0_T`t'_se))  (r(S_DRb`covs'_dep1_T`t'))    (r(S_DRb`covs'_dep1_T`t'_se)) (r(DIFF_DRb`covs'_T`t'))  (r(DIFF_DRb`covs'_T`t'_se))
					
					}
				}
				
				

			}
    noi _dots `i' 0
profiler off

	}

post `states' (`=`nsim'+1') (.)  (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
postclose `estimates'
postclose `states'

}

