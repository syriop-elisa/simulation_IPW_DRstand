# Inverse Probability Weighting and Doubly Robust Standardisation in the Relative Survival Framework

This repository contains the code required for the simulation study of the manuscript titled _Inverse probability weighting and doubly robust standardisation
in the relative survival framework_ by Syriopoulou et al. (2020).

The code for the simulation study is included in the file [`simulation_code.do`](https://github.com/syriop-elisa/simulation_IPW_DRstand/blob/master/simulation_code.do).

The main aim of this simulation study is to assess how sensitive point estimates of relative survival obtained from regression standardisation,
inverse probability weighting (IPW) and doubly robust standardisation are to model misspecification. 
As a secondary aim, different ways for obtaining standard errors for the point estimates are explored. 
For regression standardisation, standard errors are obtained by using either the delta method or M-estimation. 
For IPW and doubly robust standardisation, standard errors are obtained with the delta method while using robust clustered standard errors.  
