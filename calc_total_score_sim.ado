/* _____________________________________________________________________________________

Program: calc_total_score_sim

Purpose: Calculate the composite score of the scale at a specified 
		 wave using one of two rules.
		 Rule 1: Calculate row means if all items are observed
		 Rule 2: Calculate row means according to PedsQL guidelines
	
Written by R Mainzer, June 2020

_________________________________________________________________________________ */


capture program drop calc_total_score_sim

program define calc_total_score_sim

	syntax, wave(integer) rule(integer)
	
	* Count number of missing PedsQL items at wave for each individual
	egen w`wave'_miss = rmiss(y_`wave'_*)
	
	* Reverse score items
	foreach var of varlist y_`wave'_* {
		qui recode `var' (1=100) (2=75) (3=50) (4=25) (5=0), gen(rs_`var')
	}
	
	* Rule 1: Calculate row means if all items are observed at given wave
	if `rule' == 1 {
		qui egen y_`wave' = rowmean(rs_y_`wave'_*) if w`wave'_miss == 0
	}
	
	* Rule 2: Calculate row means if at least half the items are observed
	if `rule' == 2 {
		if `wave' == 1 {
			qui egen y_1 = rowmean(rs_y_1_*) if w1_miss <= 10
		}
		
		if `wave' != 1 {
			qui egen y_`wave' = rowmean(rs_y_`wave'_*) if w`wave'_miss <= 11
		}
	}
	
	drop w`wave'_miss 
		
end
