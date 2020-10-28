/* _____________________________________________________________________________________

Program: an_SimData

Purpose: Analyse data using different multiple imputation strategies

Input:
	miss_scen: 
		1 - realistic missingness scenario, 
		2 - extreme missingness scenario
	strat: 
		0 - complete case analysis
		1 - item level imputation using items from other scales as auxiliary 
			variables
		2 - scale level imputation using scale scores from other waves as 
			auxiliary variables
		3 - item level imputation using scale scores from other waves as 
			auxiliary variables
		4 - item level imputation using the same items from other 
		    waves and total scores from all waves as auxiliary variables in 
		    the imputation model
		5 - item level imputation using principal components, derived from
			the items of scales measured at other waves, as auxiliary 
			variables
	method: 
		"cc" - complete case analysis
		"mvn" - multivariate normal imputation
		"fcs" - fully conditional specification where items are imputed using (regress)
	simno_start: number of simulated data set to start with
	simno_end: 	 number of simulated data set to end with
	rule: 		 
		0 - skip computation of scale score (default)
		1 - compute scale score if no items are missing (strat 2 - 4) 
		2 - compute scale score if <50% of items are missing (strat 2 - 4)

Output: result file of estimated quantities
 
Written by R Mainzer, April 2020
Based on code written by JA for Masters project 

_________________________________________________________________________________ */

capture program drop an_SimData

program define an_SimData

	syntax, miss_scen(integer) strat(integer) method(string) simno_start(integer) ///
		simno_end(integer) [rule(real 0)]
	
	* Start timer
	timer on 1
	
	* Set up postfile for results 
	capture postfile results
	else {
		postfile results replication rc b se_b lower_lim_b upper_lim_b mean ///
		se_mean lower_lim_mean upper_lim_mean median se_median ///
		lower_lim_median upper_lim_median ///
		using an_SimData_miss_scen`miss_scen'_strat`strat'_rule`rule'_method`method'_datasets`simno_start'to`simno_end', replace
	}
	
	forvalues s = `simno_start'(1)`simno_end' {
		noi di "Iteration number = `s'"
		noi di "*********************************************************"
		
		* Load data
		use "analysis_datafiles/sim/SimData_miss_scen`miss_scen'_`s'.dta", clear
	
		* Complete case analysis
		if `strat' == 0 {
			
			qui{
				
				* Calculate scale score at wave 4 according to rule
				calc_total_score_sim, wave(4) rule(`rule')
			
				* Regression of y_4 on x
				regress y_4 x
				matrix table = r(table)
				scalar b = table[1,1]
				scalar se_b = table[2,1]
				scalar ll_b = table[5,1]
				scalar ul_b = table[6,1]
		
				* Mean of y_4
				mean y_4
				matrix table = r(table)
				scalar mean = table[1,1]
				scalar se_mean = table[2,1]
				scalar ll_mean = table[5,1]
				scalar ul_mean = table[6,1]
		
				* Median of y_4
				qreg y_4
				matrix table=r(table)
				scalar median=table[1,1]
				scalar se_median = table[2,1]
				scalar ll_median = table[5,1]
				scalar ul_median = table[6,1]		
			}
		
			post results (`s') (.) (b) (se_b) (ll_b) (ul_b) (mean) (se_mean) ///
			(ll_mean) (ul_mean) (median) (se_median) (ll_median) (ul_median)
			
		}
	
		* MI analysis
		if `strat' != 0 {
		
			* MI set and register regular variables
			mi set flong
			mi register regular x z
					
			* Strategy 1: item-level imputation using items from other scales as 
			* auxiliary variables
			if `strat' == 1 {
					
				* MVN imputation
				if "`method'" == "mvn" {
					mi register imputed y_1_1 - y_4_23
					cap mi impute mvn y_1_1 - y_4_23 = x z, noisily add(40)
				}
				
				if "`method'" == "fcs" {
					mi register imputed y_1_1 - y_4_23
					cap mi impute chained (pmm, knn(10)) y_1_1 - y_4_23 = x z, ///
					noisily augment add(40)
				}
				
			}
		
		
			* Strategy 2: scale-level imputation using scale scores from other waves 
			* as auxiliary variables
			if `strat' == 2 {
				
				* Calculate scale scores at each wave
				calc_total_score_sim, wave(1) rule(`rule')
				calc_total_score_sim, wave(2) rule(`rule')
				calc_total_score_sim, wave(3) rule(`rule')				
				calc_total_score_sim, wave(4) rule(`rule')
			
				* Register scale scores as auxiliary variables
				mi register imputed y_1 y_2 y_3 y_4
			
				* MVN imputation
				if "`method'" == "mvn" {
					cap mi impute mvn y_1 y_2 y_3 y_4 = x z, noisily add(40)
				}
			
				* FCS imputation
				if "`method'" == "fcs" {
					cap mi impute chained (pmm, knn(10)) y_1 y_2 y_3 y_4 = x z, ///
					augment add(40)
				}
			
			}
	
	
			* Strategy 3: item-level imputation using scale scores from other waves as
			* auxiliary variables
			if `strat' == 3 {
				
				* Calculate scale scores at waves 1 - 3
				calc_total_score_sim, wave(1) rule(`rule')
				calc_total_score_sim, wave(2) rule(`rule')
				calc_total_score_sim, wave(3) rule(`rule')	
			
				* MVN imputation
				if "`method'" == "mvn" {
					mi register imputed y_4_1 - y_4_23 y_1 y_2 y_3
					cap mi impute mvn y_4_1- y_4_23 y_1 y_2 y_3 ///
					= x z, add(40)
				}
				
				* FCS imputation
				if "`method'" == "fcs" {
					mi register imputed y_4_1 - y_4_23 y_1 y_2 y_3
					cap mi impute chained (pmm, knn(10)) y_4_1 - y_4_23 ///
										y_1 y_2 y_3 = x z, ///
										augment add(40)
				}
			
			}
			
			
			* Strategy 4: item-level imputation using the same items from other 
			* waves and total scores from all waves as auxiliary variables in 
			* the imputation model
			if `strat' == 4 {
			
				* Calculate total scale score at each wave
				calc_total_score_sim, wave(1) rule(`rule')
				calc_total_score_sim, wave(2) rule(`rule')
				calc_total_score_sim, wave(3) rule(`rule')
				calc_total_score_sim, wave(4) rule(`rule')

				* Drop reverse scored items (will be re-computed after mi)
				drop rs*

				* Re-label items that were not measured at all four waves to 
				* distinguish them from each other
				rename (y_1_19 y_1_20 y_1_21) (y_1_19a y_1_22 y_1_23)
				rename (y_2_19 y_3_19 y_4_19)(y_2_19b y_3_19b y_4_19b)

				* Register scale scores at all waves for imputation
				mi register imputed y_1 y_2 y_3 y_4

				* Register items for imputation
				mi register imputed y_1_1 - y_4_23 
				
				* FCS imputation
				cap mi impute chained ///
				  (regress, omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) /// 
								 y_1 y_2 y_3 y_4 ///
				  (pmm, knn(10) omit( y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_1 ///
				  (pmm, knn(10) omit(y_*_1  y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_2 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2  y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_3 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3  y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_4 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4  y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_5 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5  y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_6 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6  ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_7 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								  y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_8 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8  y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_9 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9  y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_10 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10  y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_11 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11  y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_12 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12  ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_13 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								  y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_14 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14  y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_15 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15  y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_16 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16  y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_17 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17  y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_18 ///
				 (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19b ///
								  y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_19a ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								  y_*_20 y_*_21 y_*_22 y_*_23)) ///
								 y_*_19b ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b  y_*_21 y_*_22 y_*_23)) ///
								 y_*_20 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20  y_*_22 y_*_23)) ///
								 y_*_21 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21  y_*_23)) ///
								 y_*_22 ///
				  (pmm, knn(10) omit(y_*_1 y_*_2 y_*_3 y_*_4 y_*_5 y_*_6 y_*_7 ///
								 y_*_8 y_*_9 y_*_10 y_*_11 y_*_12 y_*_13 ///
								 y_*_14 y_*_15 y_*_16 y_*_17 y_*_18 y_*_19a ///
								 y_*_19b y_*_20 y_*_21 y_*_22 )) ///
								 y_*_23 ///								 
				  = x z, add(40) augment
				
				* Drop total score at wave 6 (will be computed from items)
				drop y_4
			
			}
			
			
			* Strategy 5: item-level imputation using principal components, derived 
			* from the items of scales measured at other waves, as auxiliary variables
			if `strat' == 5 {
			
				*** Calculate principal components *** 
				
				* Fill in missing items at waves 1 - 3 using one fcs imputation
				mi register imputed y_1_1 - y_4_23
				cap mi impute chained (pmm, knn(10)) y_1_1 - y_4_23 = x z, ///
					augment add(1)
				
				if _rc == 0 {
						
					* Extract principal components 
					qui{
						pca y_1_1 - y_3_23 if _mi_m == 1, com(7) 
						mat eigenvalues = e(Ev)
						mat eigenvectors_components = e(L)
						gen ExpVar = e(rho)
		
						* predict the retained component scores
						predict pc1 - pc7 if _mi_m == 1, score
					}
				
					* Drop imputed data. Retain principal components and original data. 
					gsort _mi_id _mi_m 
					foreach v of varlist pc1 - pc7 {
						replace `v' = `v'[_n+1] if missing(`v')
					}
					mi extract 0, clear
				
					*** Do multiple imputation ***
				
					* mi register
					mi set flong
					mi register regular x z pc1 - pc7
				
					* MVN imputation
					if "`method'" == "mvn" {
						mi register imputed y_4_1 - y_4_23
						cap mi impute mvn y_4_1 - y_4_23 ///
						= pc1 - pc7 x z, add(40)
					}
					
					* FCS imputation
					if "`method'" == "fcs" {
						mi register imputed y_4_1 - y_4_23
						cap mi impute chained (pmm, knn(10)) y_4_1 - y_4_23 ///
						= pc1 - pc7 x z, add(40) augment
					}
				}			
			}
			
			display "return code = " _rc
		
			if _rc == 0 {	
				
				* Calculate total scale score for strategies 1, 3 & 5 
				* (rs_`var' and y_4 not yet defined for these strategies)
				if `strat' != 2 {
				
					* Create reverse scored items
					foreach var of varlist y_4_1 - y_4_23 {
						gen rs_`var' = 125 - 25 * `var'
					}
					
					* Create total scale score
					egen y_4 = rowmean(rs_y_4_1 - rs_y_4_23) if _mi_id != 0
					
				}
			
				qui{ 
					* Regression coeff
					mi estimate: regress y_4 x
					matrix table=r(table)
					scalar b=table[1,1]
					scalar se_b = table[2,1]
					scalar ll_b = table[5,1]
					scalar ul_b = table[6,1]
				
					* Mean
					mi estimate: mean y_4
					matrix table=r(table)
					scalar mean=table[1,1]
					scalar se_mean = table[2,1]
					scalar ll_mean = table[5,1]
					scalar ul_mean = table[6,1]
				
					* Median
					mi estimate: qreg y_4
					matrix table=r(table)
					scalar median=table[1,1]
					scalar se_median = table[2,1]
					scalar ll_median = table[5,1]
					scalar ul_median = table[6,1]
				}
				
				post results (`s') (.) (b) (se_b) (ll_b) (ul_b) (mean) (se_mean) ///
				(ll_mean) (ul_mean) (median) (se_median) (ll_median) (ul_median)
			}
		
			if _rc > 0 {
				post results (`s') (_rc) (.) (.) (.) (.) (.) (.) (.) (.) ///
				(.) (.) (.) (.)
			}
		
		}
		
	}
	
	* Close postfile
	postclose results
	
	* Stop timer
	timer off 1
	qui timer list 1
	display "Time taken = " r(t1)
	timer clear 1
	
end
