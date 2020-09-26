
global folder "C:\JMIS\UCI\RA-work\MW\ILRR\data"

cd $folder

** **************************************************** COUNTY VERSION
* difference is period only up to 2006q2

use qcew2, clear
keep county
duplicates drop
sort county
save countylist, replace

gen state = floor(county/1000)
sort state county
save countystate, replace

use borderpairs, clear
keep ctlcounty
duplicates drop
rename ctlcounty county
save a, replace

* set of border counties
use qcew2, clear
gen state = floor(county/1000)
merge m:1 county using a

eststo e1_A: qui areg lneprest lnmw lneppriv i.date if priv==1        , absorb(county) cluster(state)
eststo e1_B: qui areg lneprest lnmw lneppriv i.date if priv==1 [w=pop], absorb(county) cluster(state)

keep if _merge==3
drop _merge
save qcew2_border, replace

use qcew2_border, clear
keep county
duplicates drop
save bordercty_1, replace
drop if county==9013
drop if county==17001
save bordercty_2, replace

use bordercty_1, clear
save border, replace

use borderpairs, clear
rename ctlcounty county
merge m:1 county using border, gen(border1)
rename county cty
rename trtcounty county
merge m:1 county using border, gen(border2)
keep if border1==3 & border2==3
drop border*
rename county trtcounty
rename cty county
save borderpair, replace
drop trtcounty
duplicates drop
save borderpaircty, replace

gen state = floor(county/1000)
sort state county
save borderpairstate, replace

* set of border counties that are part of a border county pair
use qcew2, clear
gen state = floor(county/1000)
merge m:1 county using borderpaircty
eststo e2_A: qui areg lneprest lnmw lneppriv i.date if priv==1 & _merge==3        , absorb(county) cluster(state)
eststo e2_B: qui areg lneprest lnmw lneppriv i.date if priv==1 & _merge==3 [w=pop], absorb(county) cluster(state)

keep if _merge==3
drop _merge
save qcew2_borderpair, replace

use qcew2_borderpair, clear
keep county
duplicates drop
save borderpaircty_1, replace
drop if county==9013
drop if county==17001
save borderpaircty_2, replace


use borderpaircty, clear
levelsof county, local(bcountylist)
* ********************* contig donor pool (treated border counties with cross-state neighboring border counties)
foreach k in `bcountylist' {
    use if trtcounty==`k' using borderpair, clear
    merge 1:m county using qcew2, keep(match) nogen
    collapse (mean) lneprest lnmw lneppriv, by(trtcounty date)
    rename trtcounty county
    replace county = 100000 + county
    save qcew2_`k', replace
}
clear
foreach k in `bcountylist' {
    append using qcew2_`k'
    erase qcew2_`k'.dta
}
save qcew2_cr, replace
append using qcew2
gen county2 = county
replace county = . if county2>100000

gen id = county2
replace id = county2-100000 if county>100000
egen exp = group(date id)

gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen county3 = group(treated control)
save qcew2_cr, replace

* ********************* borderpair donor pool, has to be counties in other states
use qcew2_borderpair, clear
levelsof county, local(countylist)
foreach k in `countylist' {
    preserve
    collapse (mean) lneprest lnmw lneppriv if state~=floor(`k'/1000), by(date)
    gen county = 100000+`k'
    save qcew2_borderpair_`k', replace
    restore
}

clear
foreach k in `countylist' {
    append using qcew2_borderpair_`k'
    erase qcew2_borderpair_`k'.dta
}
save qcew2_borderpair_c, replace
append using qcew2_borderpair
gen county2 = county
replace county = . if county2>100000

gen id = county2
replace id = county2-100000 if county>100000
egen exp = group(date id)

gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

save qcew2_borderpair_c, replace

* ********************* full donor pool, has to be counties in other states
use qcew2, clear
gen state = floor(county/1000)
levelsof county, local(countylist)
foreach k in `countylist' {
    preserve
    collapse (mean) lneprest lnmw lneppriv if state~=floor(`k'/1000), by(date)
    gen county = 100000+`k'
    save qcew2_`k', replace
    restore
}

clear
foreach k in `countylist' {
    append using qcew2_`k'
    erase qcew2_`k'.dta
}
save qcew2_c, replace
append using qcew2
gen county2 = county
replace county = . if county2>100000

gen id = county2
replace id = county2-100000 if county>100000
egen exp = group(date id)

gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

save qcew2_c, replace



* ******************************************** all county sample
use alldmw_sample, clear
drop if state>100
rename id state
joinby state using countystate
drop state state2
rename county id
expand 2
sort trtdate id date
by trtdate id date: gen number = _n
gen county2 = id
replace county2 = 100000+id if number==2
drop number
save alldmw_samplecty, replace

* all MW experiments less 90q2
use alldmw_samplecty, clear
keep if trtdate<=yq(2005,3)
drop if trtdate==yq(1990,2)
merge m:1 date id county2 using qcew2_c, keep(match) keepusing(lneprest lnmw lneppriv pop)
egen exp = group(trtdate id)

gen county = county2 if county2<100000
gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen county3 = group(treated control)

egen group = group(exp county2)

tab county if lneppriv==.
gen priv = (id~=9013 & id~=17001)
eststo e0_A: qui areg lneprest lnmw lneppriv i.date if priv==1        , cluster(group) absorb(group)
eststo e0_B: qui areg lneprest lnmw lneppriv i.date if priv==1 [w=pop], cluster(group) absorb(group)


* ******************************************** border pair sample
use alldmw_sample, clear
drop if state>100
rename id state
joinby state using borderpairstate
drop state state2
rename county id
expand 2
sort trtdate id date
by trtdate id date: gen number = _n
gen county2 = id
replace county2 = 100000+id if number==2
drop number
save alldmw_samplecty, replace

* all MW experiments less 90q2, counties in other states
use alldmw_samplecty, clear
keep if trtdate<=yq(2005,3)
drop if trtdate==yq(1990,2)
merge m:1 date id county2 using qcew2_borderpair_c, keep(match) keepusing(lneprest lnmw lneppriv pop)
egen exp = group(trtdate id)

gen county = county2 if county2<100000
gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen county3 = group(treated control)

gen poptrt = pop if treated==county
sort id date county
by id date: replace poptrt = poptrt[_n-1] if poptrt==.
egen group = group(exp county2)

tab county if lneppriv==.
gen priv = (id~=9013 & id~=17001)
eststo e4_Ax: qui areg lneprest lnmw lneppriv i.date if priv==1           , cluster(group) absorb(group)
eststo e4_Bx: qui areg lneprest lnmw lneppriv i.date if priv==1 [w=poptrt], cluster(group) absorb(group)


* all MW experiments less 90q2, neighbor
use alldmw_samplecty, clear
keep if trtdate<=yq(2005,3)
drop if trtdate==yq(1990,2)
merge m:1 date id county2 using qcew2_cr, keep(match) keepusing(lneprest lnmw lneppriv pop)
egen exp = group(trtdate id)

gen county = county2 if county2<100000
gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen county3 = group(treated control)

gen poptrt = pop if treated==county
sort id date county
by id date: replace poptrt = poptrt[_n-1] if poptrt==.
egen group = group(exp county2)

tab county if lneppriv==.
gen priv = (id~=9013 & id~=17001)
eststo e4_A: qui areg lneprest lnmw lneppriv i.date if priv==1           , cluster(group) absorb(group)
eststo e4_B: qui areg lneprest lnmw lneppriv i.date if priv==1 [w=poptrt], cluster(group) absorb(group)


* ******************************************** border pair sample, paired with all border counties in other states
use alldmw_sample, clear
drop if state>100
rename id state
joinby state using borderpairstate
drop state state2
rename county id
expand 2
sort trtdate id date
by trtdate id date: gen number = _n
gen county2 = id
replace county2 = 100000+id if number==2
drop number
save alldmw_samplecty, replace

* all MW experiments less 90q2, counties in other states
use alldmw_samplecty, clear
keep if trtdate<=yq(2005,3)
drop if trtdate==yq(1990,2)
merge m:1 date id county2 using qcew2_c, keep(match) keepusing(lneprest lnmw lneppriv pop)
egen exp = group(trtdate id)

gen county = county2 if county2<100000
gen treated = county
gen control = exp if county==.
recode treated (.=0)
recode control (.=0)
egen county3 = group(treated control)

gen poptrt = pop if treated==county
sort id date county
by id date: replace poptrt = poptrt[_n-1] if poptrt==.
egen group = group(exp county2)

tab county if lneppriv==.
gen priv = (id~=9013 & id~=17001)
eststo e3_A: qui areg lneprest lnmw lneppriv i.date if priv==1           , cluster(group) absorb(group)
eststo e3_B: qui areg lneprest lnmw lneppriv i.date if priv==1 [w=poptrt], cluster(group) absorb(group)





* sets of sample for final analysis
use if date~=yq(1990,2) & date<=yq(2006,2)-3 using alldmw, clear
joinby state using borderpairstate
drop state
sort date county
keep if (county~=9013 & county~=17001)
save alldmw_borderpair_cty, replace

use if date~=yq(1990,2) & date<=yq(2006,2)-3 using alldmw, clear
joinby state using countystate
drop state
sort date county
keep if (county~=9013 & county~=17001)
save alldmw_allcounty_cty, replace



* UNWEIGHTED
* specify regression to obtain the residuals
local reg_5 "lneprest lnmw lneppriv i.date if priv==1, absorb(county) cluster(state)"
local reg_6 "lneprest      lneppriv i.date if priv==1, absorb(county) cluster(state)"

forval m = 5/6 {

    use qcew2, clear
//  use qcew2_border, clear
//  use qcew2_borderpair, clear
    gen state = floor(county/1000)
    areg `reg_`m''
    predict resid, resid
    keep date county resid
    xtset county date
    save resid_`m', replace

}

clear
input null
0
end
save temp_5, replace
save temp_6, replace

forval m = 5/6 {

    use alldmw_borderpair_cty, replace
//  use alldmw_allcounty_cty, clear
    summ
    local total = r(N)
    forval i = 1/`total' {
        use alldmw_borderpair_cty, clear
//      use alldmw_allcounty_cty, clear
        keep in `i'
        levelsof county, local(treated)
        levelsof date, local(trtdate)
//* for top50 border/pair counties
        rename date trtdate
        * need to only have previous 4 periods
        expand 4
        gen date = trtdate - _n
        merge 1:1 county date using resid_`m', keep(match) nogen
        rename county trtcounty
        rename resid residtrt
        merge 1:m date using resid_`m', keep(match) nogen
        gen spe = (resid-residtrt)^2
        collapse (mean) spe, by(trtdate trtcounty county)
        rename spe mspe
        sort mspe
        * take out other counties in same state
        gen state = floor(county/1000)
        drop if state==floor(`treated'/1000)
        * 51 if includes trtcounty itself, but not
        keep in 1/50
        save top50_`trtdate'_`treated', replace
        use alldmw_borderpair_cty, clear
//      use alldmw_allcounty_cty, clear
        keep in `i'
        merge 1:1 county using top50_`trtdate'_`treated', keep(using)
        erase top50_`trtdate'_`treated'.dta
//* for all border/pair counties
//      merge 1:1 county using bordercty, keep(using)
//      merge 1:1 county using borderpaircty, keep(using)
        levelsof county, local(control)
        use resid_`m', clear
        local trtdateA = `trtdate' - 4
        local trtdateB = `trtdate' - 3
        local trtdateC = `trtdate' - 2
        local trtdateD = `trtdate' - 1
        keep if date>=`trtdateA' & date<=`trtdate'
        synth resid resid(`trtdateA') resid(`trtdateB') resid(`trtdateC') resid(`trtdateD'), trperiod(`trtdate') trunit(`treated') counit(`control') xperiod(`trtdateA'/`trtdateD') mspeperiod(`trtdateA'/`trtdateD') customV(1 1 1 1) keep(sc_`trtdate'_`treated') replace
        use sc_`trtdate'_`treated', clear
        rename _Co_Number county
        rename _W_Weight weight
        drop _Y_treated _Y_synthetic _time
        gen trtdate = `trtdate'
        gen trtcounty = `treated'
        save, replace
        use temp_`m', clear
        append using sc_`trtdate'_`treated'
        save temp_`m', replace
        erase sc_`trtdate'_`treated'.dta
    }

    use temp_`m', clear
    drop in 1
    drop null
    save temp_`m', replace

}

    use alldmw_borderpair_cty, clear
//  use alldmw_allcounty_cty, clear
    rename  date trtdate
    rename county trtcounty
    cross using countylist
    expand 8
    sort trtdate trtcounty county
    by trtdate trtcounty county: gen time = _n-4
    gen date = trtdate - time
    drop time
    format date %tq
    save listmerge, replace


forval m = 5/6 {

    use listmerge, clear
    merge m:1 date county using qcew2, keep(match master) nogen
    merge m:1 trtdate trtcounty county using temp_`m', keep(match master) nogen
    replace weight = 1 if trtcounty==county
    save dup_`m', replace
    
    keep if trtcounty==county
    save dup1_`m', replace
    
    * all MW experiments, Abadie weights
    use dup_`m', clear
    drop if trtcounty==county
    foreach j in lneprest lnmw lneppriv {
        gen wt`j' = weight*`j'
    }
    collapse (sum) wt*, by(trtcounty trtdate date)
    rename wt* *
    * doesn't work, but should it?
    * collapse (sum) lneprest lnmw popteen [w=weight], by(trtcounty trtdate date)
    save dup0_`m', replace
    append using dup1_`m'
    
    * id is treatment county
    * county identifies treatment county, missing for control obs (51 FE)
    * county2 identifies each treatment county and each control county associated with each treatment county (102 FE)
    * county3 identifies each treatment county and each control county associated with each experiment (51 FE + # of experiments)
    
    egen exp = group(trtdate trtcounty)
    gen county2 = county
    replace county2 = trtcounty+100000 if county==.
    
    gen treated = county
    gen control = exp if county==.
    recode treated (.=0)
    recode control (.=0)
    egen county3 = group(treated control)
    
    gen poptrt = pop if treated==county
    sort exp date county
    by exp date: replace poptrt = poptrt[_n-1] if poptrt==.
    egen group = group(exp county2)
    
    drop priv
    gen priv = (trtcounty~=9013 & trtcounty~=17001)
    save dupresid_`m', replace

    eststo e`m'_A: qui areg lneprest lnmw lneppriv i.date if priv==1 & trtdate<=yq(2005,3)           , cluster(group) absorb(group)
}

* WEIGHTED
* specify regression to obtain the residuals
local reg_5 "lneprest lnmw lneppriv i.date if priv==1 [w=pop], absorb(county) cluster(state)"
local reg_6 "lneprest      lneppriv i.date if priv==1 [w=pop], absorb(county) cluster(state)"

forval m = 5/6 {

    use qcew2, clear
//  use qcew2_border, clear
//  use qcew2_borderpair, clear
    gen state = floor(county/1000)
    areg `reg_`m''
    predict resid, resid
    keep date county resid
    xtset county date
    save resid_`m', replace

}

clear
input null
0
end
save temp_5, replace
save temp_6, replace

forval m = 5/6 {

    use alldmw_borderpair_cty, replace
//  use alldmw_allcounty_cty, clear
    summ
    local total = r(N)
    forval i = 1/`total' {
        use alldmw_borderpair_cty, clear
//      use alldmw_allcounty_cty, clear
        keep in `i'
        levelsof county, local(treated)
        levelsof date, local(trtdate)
//* for top50 border/pair counties
        rename date trtdate
        * need to only have previous 4 periods
        expand 4
        gen date = trtdate - _n
        merge 1:1 county date using resid_`m', keep(match) nogen
        rename county trtcounty
        rename resid residtrt
        merge 1:m date using resid_`m', keep(match) nogen
        gen spe = (resid-residtrt)^2
        collapse (mean) spe, by(trtdate trtcounty county)
        rename spe mspe
        sort mspe
        * take out other counties in same state
        gen state = floor(county/1000)
        drop if state==floor(`treated'/1000)
        * 51 if includes trtcounty itself, but not
        keep in 1/50
        save top50_`trtdate'_`treated', replace
        use alldmw_borderpair_cty, clear
//      use alldmw_allcounty_cty, clear
        keep in `i'
        merge 1:1 county using top50_`trtdate'_`treated', keep(using)
        erase top50_`trtdate'_`treated'.dta
//* for all border/pair counties
//      merge 1:1 county using bordercty, keep(using)
//      merge 1:1 county using borderpaircty, keep(using)
        levelsof county, local(control)
        use resid_`m', clear
        local trtdateA = `trtdate' - 4
        local trtdateB = `trtdate' - 3
        local trtdateC = `trtdate' - 2
        local trtdateD = `trtdate' - 1
        keep if date>=`trtdateA' & date<=`trtdate'
        synth resid resid(`trtdateA') resid(`trtdateB') resid(`trtdateC') resid(`trtdateD'), trperiod(`trtdate') trunit(`treated') counit(`control') xperiod(`trtdateA'/`trtdateD') mspeperiod(`trtdateA'/`trtdateD') customV(1 1 1 1) keep(sc_`trtdate'_`treated') replace
        use sc_`trtdate'_`treated', clear
        rename _Co_Number county
        rename _W_Weight weight
        drop _Y_treated _Y_synthetic _time
        gen trtdate = `trtdate'
        gen trtcounty = `treated'
        save, replace
        use temp_`m', clear
        append using sc_`trtdate'_`treated'
        save temp_`m', replace
        erase sc_`trtdate'_`treated'.dta
    }

    use temp_`m', clear
    drop in 1
    drop null
    save temp_`m', replace

}

    use alldmw_borderpair_cty, clear
//  use alldmw_allcounty_cty, clear
    rename  date trtdate
    rename county trtcounty
    cross using countylist
    expand 8
    sort trtdate trtcounty county
    by trtdate trtcounty county: gen time = _n-4
    gen date = trtdate - time
    drop time
    format date %tq
    save listmerge, replace


forval m = 5/6 {

    use listmerge, clear
    merge m:1 date county using qcew2, keep(match master) nogen
    merge m:1 trtdate trtcounty county using temp_`m', keep(match master) nogen
    replace weight = 1 if trtcounty==county
    save dup_`m', replace
    
    keep if trtcounty==county
    save dup1_`m', replace
    
    * all MW experiments, Abadie weights
    use dup_`m', clear
    drop if trtcounty==county
    foreach j in lneprest lnmw lneppriv {
        gen wt`j' = weight*`j'
    }
    collapse (sum) wt*, by(trtcounty trtdate date)
    rename wt* *
    * doesn't work, but should it?
    * collapse (sum) lneprest lnmw popteen [w=weight], by(trtcounty trtdate date)
    save dup0_`m', replace
    append using dup1_`m'
    
    * id is treatment county
    * county identifies treatment county, missing for control obs (51 FE)
    * county2 identifies each treatment county and each control county associated with each treatment county (102 FE)
    * county3 identifies each treatment county and each control county associated with each experiment (51 FE + # of experiments)
    
    egen exp = group(trtdate trtcounty)
    gen county2 = county
    replace county2 = trtcounty+100000 if county==.
    
    gen treated = county
    gen control = exp if county==.
    recode treated (.=0)
    recode control (.=0)
    egen county3 = group(treated control)
    
    gen poptrt = pop if treated==county
    sort exp date county
    by exp date: replace poptrt = poptrt[_n-1] if poptrt==.
    egen group = group(exp county2)
    
    drop priv
    gen priv = (trtcounty~=9013 & trtcounty~=17001)
    save dupresid_`m', replace

    eststo e`m'_B: qui areg lneprest lnmw lneppriv i.date if priv==1 & trtdate<=yq(2005,3) [w=poptrt], cluster(group) absorb(group)
}

esttab e1_A e2_A e3_A e4_A e5_A e6_A using $folder\output\Table8a.csv, keep(lnmw) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ar2 r2 obslast
esttab e1_B e2_B e3_B e4_B e5_B e6_B using $folder\output\Table8a.csv, keep(lnmw) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ar2 r2 obslast

