/* _____________________________________________________________________________________

 	Do-file title:		cr_SimData
	
	The purpose of this program is to create simulated data sets based on the 
	LSAC data.
	
	Inputs:
	miss_scen: missingness scenario 
				1 - CDM-LSAC
				2 - CDM-inflated
				3 - CDM-extreme
	simno: number of simulated data sets to create
	seed: starting seed
			 
	Datasets used:	LSAC.dta
	Datasets created: 	SimData_miss_scen`miss_scen'_`s'.dta
						                   	
	* Written by R Mainzer, April 2020
	* Based on code written by JA for Masters project
	
	UPDATES: 
	02-10-2020 - Added missingness scenario 3 (CDM-extreme)
_________________________________________________________________________________ */

capture program drop cr_SimData

program define cr_SimData

	syntax, miss_scen(integer) simno(integer) seed(integer)

	set seed `seed'
	
	*** Set up parameters for generating data sets ***

	* Obtain variables for simulation study
	use PedsQL_a1_w1 - PedsQL_b3c_w1 PedsQL_a1_w2 - PedsQL_b3c_w2 ///
		PedsQL_a1_w3 - PedsQL_b3c_w3 PedsQL_a1_w4 - PedsQL_b3c_w4 ///
		bmi_w1 sdq_w4 ///
	using "analysis_datafiles\LSAC.dta", clear
	order PedsQL*

	* Obtain correlation matrix of data
	matpwcorr PedsQL_a1_w1 - sdq_w4
	mat R = corr

	* Obtain mean and sd of bmi_w1 and sdq_w4
	foreach var of varlist bmi_w1 sdq_w4{
	qui summarize `var'
	scalar m_`var' = r(mean)
	scalar sd_`var' = r(sd)
	}

	* Create a 92x1 vector of means with the first 90 elements = 0 and the last two
	* elements = mean of bmi_w1 and sdq_w4
	matrix a = J(90,1,0)		
	matrix b = (m_bmi_w1 \ m_sdq_w4)	
	matrix m = a \ b				

	* Create a 92x1 vector of sd's with the first 90 elements = 1 and the the last 
	* two elements = sd of bmi_w1 and sdq_w4
	matrix c = J(90,1,1)		
	matrix d = (sd_bmi_w1 \ sd_sdq_w4)	
	matrix sd = c \ d				

	* Drop variables, matrices retained in memory
	drop _all

	* Specify cut points
	scalar a = invnormal(0.4)
	scalar b = invnormal(0.4 + 0.3)
	scalar c = invnormal(0.4 + 0.3 + 0.2)
	scalar d = invnormal(0.4 + 0.3 + 0.2 + 0.07)
	
	*** Generate data sets ***
	
	forvalues s = 1(1)`simno'{

		* Draw n = 1000 observations on 92 variables from a MVN distribution
		set obs 1000
		drawnorm y1-y92, means(m) sds(sd) corr(R)
		gen id = _n

		* Categorise items to create categorical variables
		foreach v of varlist y1-y90 {
			gen `v'_cat = 1
			replace `v'_cat = 2 if `v'>= a & `v'< b
			replace `v'_cat = 3 if `v'>= b & `v'< c
			replace `v'_cat = 4 if `v'>= c & `v'< d
			replace `v'_cat = 5 if `v'>= d
		}

		rename (y1_cat-y21_cat) (y_1_#), addnumber
		rename (y22_cat-y44_cat) (y_2_#), addnumber
		rename (y45_cat-y67_cat) (y_3_#), addnumber
		rename (y68_cat-y90_cat) (y_4_#), addnumber
	
		gen x = y91
		gen z = y92

		drop y1-y92

		*** Set data to missing ***

		*Rescale x and z
		egen sd_x = sd(x)
		gen x_scale = x/sd_x

		egen sd_z = sd(z)
		gen z_scale = z/sd_z

		*** Set data to missing according to missingness scenario ***
		 
		if `miss_scen' == 1 {
			gen logit_p_cases = -3.74 + ln(1.2) * x_scale + ln(1.2) * z_scale
			gen logit_p_items = -7.31 + ln(1.2) * x_scale + ln(1.2) * z_scale
		}

		if `miss_scen' == 2 {
			gen logit_p_cases = -2.84 + ln(1.2) * x_scale + ln(1.2) * z_scale
			gen logit_p_items = -6.58 + ln(1.2) * x_scale + ln(1.2) * z_scale
		}
		
		if `miss_scen' == 3 {
			gen logit_p_cases = -8.57 + ln(2) * x_scale + ln(2) * z_scale
			gen logit_p_items = -12.31 + ln(2) * x_scale + ln(2) * z_scale  				
		}
		
		* Generate probabilities used to create missing cases and items
		gen p_cases = (1/(1 + exp(-(logit_p_cases))))
		gen p_items = (1/(1 + exp(-(logit_p_items))))
	
		* Set data missing
		forvalues i = 1/4 {
			
			* Set up indicators for missing cases at each wave
			gen miss_cases_w`i' = 0
			qui replace miss_cases_w`i' = 1 if runiform() < p_cases 
		
			* Set data missing at wave 1
			if `i' == 1 {
				foreach v of varlist y_`i'_1 - y_`i'_21 {
					qui replace `v' = . if miss_cases_w`i' == 1
					qui replace `v' = . if runiform() < p_items & miss_cases_w`i' == 0	
				}
			}
			
			* Set data missing at waves 2 - 4
			if `i' != 1 {
				foreach v of varlist y_`i'_1-y_`i'_23 {
					qui replace `v' = . if miss_cases_w`i' == 1
					qui replace `v' = . if runiform() < p_items & miss_cases_w`i' == 0	
				}
			}
		}
		
		* Remove extra variables from data
		drop sd_x - miss_cases_w4

		save "analysis_datafiles\sim\SimData_miss_scen`miss_scen'_`s'.dta", replace
		clear
	}
	
end
 