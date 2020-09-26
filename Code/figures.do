
global folder "C:\JMIS\UCI\RA-work\MW\ILRR\data"

cd $folder

use trtctl4, clear
expand 3
sort trtdate trtstate ctlstate
by trtdate trtstate ctlstate: gen n = _n
gen depvar = ""
replace depvar = "resid1"  if n==1
replace depvar = "resid2"  if n==2
replace depvar = "lnepteen" if n==3
drop n
* change to 5 if want to expand selection window without MW changes
append using trtctl4
replace depvar = "d1lnepteen" if depvar==""
drop if trtdate<yq(1990,1)+5 & depvar=="d1lnepteen"
* change to 8 if want to expand selection window without MW changes
append using trtctl4
replace depvar = "d4lnepteen" if depvar==""
drop if trtdate<yq(1990,1)+8 & depvar=="d4lnepteen"
save trtctl_state, replace

* ranked prediction error approach for states
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    use trtctl_state, clear
    keep if depvar=="`j'"
    expand 4
    bysort trtdate trtstate ctlstate: gen date = trtdate - _n
    rename trtstate state
    merge m:1 state date using cps_depvar, nogen keep(match) keepusing(`j')
    rename `j' `j'trt
    rename state trtstate
    rename ctlstate state
    merge m:1 state date using cps_depvar, nogen keep(match) keepusing(`j')
    rename `j' `j'ctl
    rename state ctlstate
    gen spe = (`j'trt - `j'ctl)^2
    collapse (mean) spe, by(depvar trtdate trtstate ctlstate)
    rename spe mspe
    bysort trtdate trtstate: egen rank = rank(-mspe)
    * Weibull variant
    gen one = 1
    egen N = total(one), by(depvar trtdate trtstate)
    drop one
    gen pcrank = 100*(rank)/(N+1)
    save statepctile_`j', replace
}

clear
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    append using statepctile_`j'
    erase statepctile_`j'.dta
}
save statepctile, replace

* percentile graphs
use statepctile, clear
rename ctlstate state_fips
merge m:1 state_fips using censusdivkey, nogen keepusing(censusdiv) keep(master match)
rename censusdiv ctldiv
rename state_fips ctlstate
rename trtstate state_fips
merge m:1 state_fips using censusdivkey, nogen keepusing(censusdiv) keep(master match)
rename censusdiv trtdiv
rename state_fips trtstate
gen samediv = ctldiv==trtdiv

foreach j in resid1 lnepteen d1lnepteen d4lnepteen {
    summ pcrank if samediv==1 & depvar=="`j'", detail
    local `j'_mdn = round(r(p50))
}

local depvar resid1 lnepteen d1lnepteen d4lnepteen
local title_resid1 `"`"Regression"' `"residuals"'"'
local title_lnepteen `"`"Log teen employment-to-"' `"population ratio"'"'
local title_d1lnepteen `"`"1-quarter diff. in log teen"' `"employment-to-population ratio"'"'
local title_d4lnepteen `"`"4-quarter diff. in log teen"' `"employment-to-population ratio"'"'
* "
foreach j in resid1 lnepteen d1lnepteen d4lnepteen {
    hist pcrank if samediv==1 & depvar=="`j'", width(1) start(0) lcolor(sand) title(`title_`j'', size(medlarge)) xtitle("Percentile ranking (Median = ``j'_mdn')")
    graph save statepct_`j', replace
}

graph combine statepct_resid1.gph statepct_lnepteen.gph statepct_d1lnepteen.gph statepct_d4lnepteen.gph, rows(2) ysize(6.5) xsize(9) ycommon
graph save statepct, replace
graph export $folder\output\Figure2.png, replace height(1200)







gen rmspe = sqrt(mspe)
label var rmspe "Pre-treatment RMSPE"
local trt 6 6  6 
local con 2 15 41
local dte 2001q1 2001q1 2001q1
local x: word count `trt'
forval i = 1/`x' {
    local A: word `i' of `trt'
    local B: word `i' of `con'
    local C: word `i' of `dte'
    foreach k in resid1 lnepteen d1lnepteen d4lnepteen {
        summ rmspe if depvar=="`k'" & trtdate==tq(`C') & trtstate==`A' & ctlstate==`B' 
        local `k'_`C'_`A'_`B'_mdn = r(mean)
    }
}
* 53 not included in 6_2001q1 since also had coincident MW increase

* wsize is about a tenth of the SD
foreach k in resid1 lnepteen d1lnepteen d4lnepteen {
    summ rmspe if depvar=="`k'" & trtdate==tq(2001q1) & trtstate==6
}
local depvar resid1 lnepteen d1lnepteen d4lnepteen
local wsize  .004  .011      .005       .005
local x: word count `depvar'
forval i = 1/`x' {
    local K: word `i' of `depvar'
    local B: word `i' of `wsize'
    hist rmspe if trtdate==yq(2001,1) & trtstate==6  & depvar=="`K'", width(`B') start(0) lcolor(sand) freq legend(off) title(`title_`K'', size(medlarge)) xline(``K'_2001q1_6_2_mdn', lcolor(blue) lpattern(solid)) xline(``K'_2001q1_6_15_mdn', lcolor(blue) lpattern(longdash)) xline(``K'_2001q1_6_41_mdn', lcolor(blue) lpattern(shortdash))
    graph save `K'_2001q1_6, replace
}

graph combine resid1_2001q1_6.gph lnepteen_2001q1_6.gph d1lnepteen_2001q1_6.gph d4lnepteen_2001q1_6.gph, rows(2) ysize(6.5) xsize(9) ycommon
graph save 2001q1_6, replace
graph export $folder\output\Figure1.png, replace height(1200)
* edit with AK, HI, OR lines







use trtctl4_cty2, clear
expand 4
sort trtdate trtcounty ctlcounty
by trtdate trtcounty ctlcounty: gen date = trtdate - _n
rename trtcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(resid1 resid2 lneprest)
rename resid? resid?trt
rename lneprest lnepresttrt
rename county trtcounty
rename ctlcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(resid1 resid2 lneprest)
rename resid? resid?ctl
rename lneprest lneprestctl
rename county ctlcounty
foreach j in resid1 resid2 lneprest {
    gen spe_`j' = (`j'trt - `j'ctl)^2
}
collapse (mean) spe_* contig, by(trtdate trtcounty ctlcounty)
rename spe_* mspe_*
foreach j in resid1 resid2 lneprest {
preserve
    sort trtdate trtcounty mspe_`j'
    keep trtdate trtcounty ctlcounty contig mspe_`j'
    gen depvar = "`j'"
    save trtctl_`j'_all, replace
restore
}

* change to 5 if want to expand selection window without MW changes
use trtctl4_cty2, clear
drop if trtdate<yq(1990,1)+5
expand 4
sort trtdate trtcounty ctlcounty
by trtdate trtcounty ctlcounty: gen date = trtdate - _n
rename trtcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(d1lneprest)
rename d1lneprest d1lnepresttrt
rename county trtcounty
rename ctlcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(d1lneprest)
rename d1lneprest d1lneprestctl
rename county ctlcounty
foreach j in d1lneprest {
    gen spe_`j' = (`j'trt - `j'ctl)^2
}
collapse (mean) spe_* contig, by(trtdate trtcounty ctlcounty)
rename spe_* mspe_*
foreach j in d1lneprest {
    sort trtdate trtcounty mspe_`j'
    gen depvar = "`j'"
    save trtctl_`j'_all, replace
}

* change to 8 if want to expand selection window without MW changes
use trtctl4_cty2, clear
drop if trtdate<yq(1990,1)+8
expand 4
sort trtdate trtcounty ctlcounty
by trtdate trtcounty ctlcounty: gen date = trtdate - _n
rename trtcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(d4lneprest)
rename d4lneprest d4lnepresttrt
rename county trtcounty
rename ctlcounty county
merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(d4lneprest)
rename d4lneprest d4lneprestctl
rename county ctlcounty
foreach j in d4lneprest {
    gen spe_`j' = (`j'trt - `j'ctl)^2
}
collapse (mean) spe_* contig, by(trtdate trtcounty ctlcounty)
rename spe_* mspe_*
foreach j in d4lneprest {
    sort trtdate trtcounty mspe_`j'
    gen depvar = "`j'"
    save trtctl_`j'_all, replace
}


foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    use trtctl_`j'_all, clear
    * take out 9013 and 17001 only for resid
    drop if mspe_`j'==.
    keep trtdate trtcounty ctlcounty depvar
    save trtctl_`j'_all, replace
}

clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    append using trtctl_`j'_all
    erase trtctl_`j'_all.dta
}
save trtctl_all, replace

* ranked prediction error approach for counties
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    use trtctl_all, clear
    keep if depvar=="`j'"
    expand 4
    bysort trtdate trtcounty ctlcounty: gen date = trtdate - _n
    rename trtcounty county
    merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(`j')
    rename `j' `j'trt
    rename county trtcounty
    rename ctlcounty county
    merge m:1 county date using qcew_depvar, nogen keep(match) keepusing(`j')
    rename `j' `j'ctl
    rename county ctlcounty
    gen spe = (`j'trt - `j'ctl)^2
    collapse (mean) spe, by(depvar trtdate trtcounty ctlcounty)
    rename spe mspe
    bysort trtdate trtcounty: egen rank = rank(-mspe)
    * Weibull variant
    gen one = 1
    egen N = total(one), by(depvar trtdate trtcounty)
    drop one
    gen pcrank = 100*(rank)/(N+1)
    save countypctile_`j', replace
}

clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    append using countypctile_`j'
    erase countypctile_`j'.dta
}
save countypctile, replace

* percentile graphs
use countypctile, clear
merge m:1 trtcounty ctlcounty using borderpairs, keep(master match)
gen contig = _merge==3
drop _merge

foreach j in resid1 lneprest d1lneprest d4lneprest {
    summ pcrank if contig==1 & depvar=="`j'", detail
    local `j'_mdn = round(r(p50))
}

local depvar resid1 lneprest d1lneprest d4lneprest
local title_resid1 `"`"Regression"' `"residuals"'"'
local title_lneprest `"`"Log restaurant employment-to-"' `"county population ratio"'"'
local title_d1lneprest `"`"1-quarter diff. in log restaurant"' `"employment-to-county pop. ratio"'"'
local title_d4lneprest `"`"4-quarter diff. in log restaurant"' `"employment-to-county pop. ratio"'"'
* "
foreach j in resid1 lneprest d1lneprest d4lneprest {
    hist pcrank if contig==1 & depvar=="`j'", width(1) start(0) lcolor(sand) title(`title_`j'', size(medlarge)) xtitle("Percentile ranking (Median = ``j'_mdn')")
    graph save countypct_`j', replace
}

graph combine countypct_resid1.gph countypct_lneprest.gph countypct_d1lneprest.gph countypct_d4lneprest.gph, rows(2) ysize(6.5) xsize(9) ycommon
graph save countypct, replace
graph export $folder\output\Figure3.png, replace height(1200)
