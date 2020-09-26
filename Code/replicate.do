* tell something about data
* tell something about order of running things
* tell something about output

* data: CPS output, QCEW output (or DLR's data)
* separate MW series? maybe not

* need to install estout package, synth
* code in R for county SC, code in Stata for state SC
* cluster2areg

global folder "C:\Education\University of California Irvine\Winter 2018\ECON 190BW\Data"

clear all

cd "C:\Education\University of California Irvine\Winter 2018\ECON 190BW\Data"
use CPSdata, clear
fvset base 140 date
rename epteen lnepteen
rename minwage lnminwage

* Create state time trends
gen trend1 = (date-119)/(10^1)
gen trend2 = (trend1^2)/(10^2)
gen trend3 = (trend1^3)/(10^3)
gen trend4 = (trend1^4)/(10^4)
gen trend5 = (trend1^5)/(10^5)
gen trend6 = (trend1^6)/(10^6)
gen trend7 = (trend1^7)/(10^7)

* Table 1
eststo tab1_1: qui reg lnepteen lnminwage urall rpteen i.state i.date                                     [w=popteen], cluster(state)
eststo tab1_2: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.trend1                    [w=popteen], cluster(state)
eststo tab1_3: qui reg lnepteen lnminwage urall rpteen i.state i.date                  i.censusdiv#i.date [w=popteen], cluster(state)
eststo tab1_4: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.trend1 i.censusdiv#i.date [w=popteen], cluster(state)
esttab tab1_1 tab1_2 tab1_3 tab1_4 using $folder\output\Table1a.csv, keep(lnminwage urall rpteen) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) r2 obslast

eststo tab1_t2: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend2) [w=popteen], cluster(state)
eststo tab1_t3: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend3) [w=popteen], cluster(state)
eststo tab1_t4: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend4) [w=popteen], cluster(state)
eststo tab1_t5: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend5) [w=popteen], cluster(state)
eststo tab1_t6: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend6) [w=popteen], cluster(state)
eststo tab1_t7: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.(trend1-trend7) [w=popteen], cluster(state)
esttab tab1_t2 tab1_t3 tab1_t4 tab1_t5 tab1_t6 tab1_t7 using $folder\output\Table1b.csv, keep(lnminwage urall rpteen) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) obslast

* verify Table 1 results same as NW 2011
eststo tab1_1x: qui reg lnepteen lnminwage urall rpteen i.state i.date                  if date>=yq(1994,1) & date<=yq(2007,2) [w=popteen], cluster(state)
eststo tab1_2x: qui reg lnepteen lnminwage urall rpteen i.state i.date i.state#c.trend1 if date>=yq(1994,1) & date<=yq(2007,2) [w=popteen], cluster(state)
esttab tab1_1 tab1_1x tab1_2 tab1_2x, keep(lnminwage urall rpteen) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) r2 obslast

* Table 2
forval i = 1/9 {
    eststo tab2_`i'a: qui reg lnepteen lnminwage urall rpteen i.state i.date if censusdiv==`i'                    [w=popteen], cluster(state)
    eststo tab2_`i'b: qui reg lnepteen lnminwage urall rpteen i.state i.date if censusdiv==`i' & date<=yq(2009,4) [w=popteen], cluster(state)
}
forval i = 1/9 {
    esttab tab2_`i'a tab2_`i'b using $folder\output\Table2.csv, append keep(lnminwage) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) r2 obslast
}

* Table 3, show results that take out invalid obs in 1q and 4q difference results
* also show residuals estimated with or without MW variable














use QCEWindustry_minwage_contig.dta, clear

rename lnMW lnmw
rename lnemp_rest_both lnempboth
rename lnemp_TOT lnemptotal
rename pair_id pairid
rename nonmissing_rest_both nonmiss

gen yearqtr = yq(year, quarter)
gen sampleperiod = yearqtr>=yq(1990,1) & yearqtr<=yq(2006,2)
gen sample = sampleperiod==1 & nonmiss==66
egen pairidc = group(pairid countyreal)

foreach j in lnempboth lnmw lnpop lnemptotal {
    qui areg `j' if sampleperiod, absorb(countyreal)
    predict `j'_r, residual
}

gen  spec5 = yearqtr
egen spec6 = group(pairid yearqtr)
* drop if pair_id_county==.
* egen pairidc = group(pairid countyreal)

gen state1 = real(substr(pairid, 1, 2))
gen state2 = real(substr(pairid, 7, 2))
gen st_min = min(state1, state2)
gen st_max = max(state1, state2)
egen bordersegment = group(st_min st_max)
* xtset pairidc yearqtr

* weight
gen one = 1

* Table 4 (DLR replication)
foreach i in 5 6 {
    eststo spec`i'emp : qui cluster2areg lnempboth_r lnmw_r lnpop_r              if sample==1, a(spec`i') fcluster(state_fips) tcluster(bordersegment) w(one)
    eststo spec`i'empt: qui cluster2areg lnempboth_r lnmw_r lnpop_r lnemptotal_r if sample==1, a(spec`i') fcluster(state_fips) tcluster(bordersegment) w(one)
}
esttab spec5emp spec5empt spec6emp spec6empt using $folder\output\Table4.csv, keep(lnmw_r lnpop_r lnemptotal_r) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) r2 obslast



* save sample corresponding to spec1
use QCEWindustry_minwage_all, clear

rename lnMW lnmw
rename lnemp_rest_both lnempboth
rename lnemp_TOT lnemptotal
rename nonmissing_rest_both nonmiss
rename state_fips state

gen yearqtr = yq(year, quarter)
gen sampleperiod = yearqtr>=yq(1990,1) & yearqtr<=yq(2006,2)
gen sample = sampleperiod==1 & nonmiss==66

keep if sample==1
keep lnempboth lnmw lnpop lnemptotal countyreal yearqtr pop
gen lneprest = ln(exp(lnempboth)/exp(lnpop))
gen lneppriv = ln(exp(lnemptotal)/exp(lnpop))
rename yearqtr date
rename countyreal county
tab county if lneppriv==.
gen priv = (county~=9013 & county~=17001)
save qcew2, replace









* Table 6, also add column 2 with earlier date 1997,4

use QCEWindustry_minwage_contig.dta, clear
gen yearqtr = yq(year,quarter)
format yearqtr %tq
keep if yearqtr>=yq(1990,1) & yearqtr<=yq(2006,2)

rename lnMW lnmw
rename lnemp_rest_both lnempboth
rename lnemp_TOT lnemptotal
rename state_fips state
fvset base 185 yearqtr

* setup for balanced sample identification
gen empboth = emp7221 + emp7222
bysort pair_id countyreal: egen nonmissA = count(empboth)
bysort pair_id countyreal: egen nonmissB = count(empboth) if yearqtr>=yq(1998,3)

* assign paired county's data
sort yearqtr pair_id
foreach x in "minwage" "lnmw" "lnemptotal" "lnempboth" "cntyarea" "lnpop" "nonmissA" "nonmissB" "state" "countyreal" {
    by yearqtr pair_id: gen     `x'_pairother = `x'[_n-1]
    by yearqtr pair_id: replace `x'_pairother = `x'[_n+1] if `x'_pairother==.
}

** create vars for later use
* same as federal
gen mweqfedmw = float(minwage)==float(federalmin)
* not same as neighbor
gen mwdiff = float(minwage)~=float(minwage_pairother)
egen mwdiffc = max(mwdiff), by(pair_id)

*** define conditions for full sample and subsample 98q3 onwards ***************************************************************
* specify sample conditions
local condA1 "yearqtr==yearqtr"
local condA2 "mwdiffsumA>=1 & mwdiffsumA<."
local condB1 "yearqtr>=yq(1998,3)"
local condB2 "mwdiffsumB>=1 & mwdiffsumB<."
* specify conditions for estimation
local condA   "counterfA==1"
local condA_p "counterfA==1 & cntyarea          <2000 & nonmissA          ==66"
local condA_a "counterfA==1 & cntyarea_pairother<2000 & nonmissA_pairother==66"
local condB   "counterfB==1"                                   
local condB_p "counterfB==1 & cntyarea          <2000 & nonmissB          ==32"
local condB_a "counterfB==1 & cntyarea_pairother<2000 & nonmissB_pairother==32"

foreach j in A B {

    * identify counties with state MW never above federal level
    egen counterf`j'  = min(mweqfedmw) if `cond`j'1', by(pair_id_county)

    * identify cases
    egen mwdiffsum`j' = total(mwdiff) if `cond`j'1', by(pair_id_county)
    
    * show proportion of cases where placebo MW is indeed a placebo
    forval i = 1/2 {
        tab mwdiff if `cond`j'`i'' & `cond`j'_p'
    }
    forval i = 1/2 {
        tab mwdiffc if `cond`j'`i'' & `cond`j'_p'
    }

}


foreach j in A B {
    forval i = 1/2 {
        * estimate DLR's actual MW sample and then placebo MW sample
        eststo est`j'`i'_a: qui areg lnempboth_pairother lnmw_pairother lnemptotal_pairother lnpop_pairother i.yearqtr if `cond`j'`i'' & `cond`j'_a', cluster(state) absorb(countyreal)
        eststo est`j'`i'_p: qui areg lnempboth           lnmw_pairother lnemptotal           lnpop           i.yearqtr if `cond`j'`i'' & `cond`j'_p', cluster(state) absorb(countyreal)
    }
}



esttab estA1_a estB1_a estB2_a using $folder\output\Table6a.csv, plain compress keep(lnmw_pairother) b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) parentheses
esttab estA1_p estB1_p estB2_p using $folder\output\Table6b.csv, plain compress keep(lnmw_pairother) b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) parentheses


* Table 3 (Stata runs)

* Table 5 (R runs)

* Table 7 (states)

* Table 8 (counties)

* generate census division key for later use
use CPSdata, clear
keep censusdiv state
duplicates drop
rename state state_fips
save censusdivkey, replace




