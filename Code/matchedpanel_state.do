
global folder "C:\JMIS\UCI\RA-work\MW\ILRR\data"

cd $folder

* ********************* FULL DONOR POOL

use CPSdata, clear
levelsof state, local(statelist)
foreach k in `statelist' {
    preserve
    collapse (mean) lnepteen lnminwage urall rpteen popteen if state~=`k', by(date)
    gen state = 100+`k'
    save CPSdata_`k', replace
    restore
}

clear
foreach k in `statelist' {
    append using CPSdata_`k'
    erase CPSdata_`k'.dta
}
save CPSdata_c, replace
append using CPSdata
gen state2 = state
replace state = . if state2>100

gen id = state2
replace id = state2-100 if state>100
egen exp = group(date id)

gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

gen popteentrt = popteen if treated==state
sort id date state
by id date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
save CPSdata_c, replace

* same-div
levelsof censusdiv, local(divlist)
foreach j in `divlist' {
    use CPSdata, clear
    keep if censusdiv==`j'
    levelsof state, local(statelist)
    foreach k in `statelist' {
        preserve
        collapse (mean) lnepteen lnminwage urall rpteen popteen if state~=`k', by(date)
        gen state = 100+`k'
        save CPSdata_`j'_`k', replace
        restore
    }
    clear
    foreach k in `statelist' {
        append using CPSdata_`j'_`k'
        erase CPSdata_`j'_`k'.dta
    }
    save CPSdata_c_`j', replace
}

clear
forval i = 1/9 {
    append using CPSdata_c_`i'
    erase CPSdata_c_`i'.dta
}
save CPSdata_cr, replace
append using CPSdata
gen state2 = state
replace state = . if state2>100

gen id = state2
replace id = state2-100 if state>100
egen exp = group(date id)

gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)
save CPSdata_cr, replace

* state contig
use contigstatepairs, clear
keep state1 state2
levelsof state1, local(bstatelist)
foreach j in `bstatelist' {
    use state1 state2 if state1==`j' using contigstatepairs, clear
    drop state1
    rename state2 state
    merge 1:m state using CPSdata, keep(match) nogen
    collapse (mean) lnepteen lnminwage urall rpteen popteen, by(date)
    gen state = 100+`j'
    save CPSdata_`j', replace
}
clear
foreach j in `bstatelist' {
    append using CPSdata_`j'
    erase CPSdata_`j'.dta
}
save CPSdata_b, replace
append using CPSdata
gen state2 = state
replace state = . if state2>100

gen id = state2
replace id = state2-100 if state>100
egen exp = group(date id)

gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)
save CPSdata_b, replace





cd $folder

use CPSdata, clear
* take out AK and HI since ultimately want to see state border analysis analog
eststo e1: qui areg lnepteen lnminwage urall rpteen i.date [w=pop], absorb(state) cluster(state)

* 129 clean experiments
use trtctl4, clear
drop ctlstate
duplicates drop
rename trtdate date
rename trtstate state
save trt4, replace

use trt4, clear
rename date trtdate
rename state id
expand 8
sort trtdate id
egen a = group(trtdate id)
spell a
gen time = _seq-4
gen date = trtdate - time
format date %tq
keep trtdate id date
expand 2
sort trtdate id date
by trtdate id date: gen number = _n
gen state2 = id
replace state2 = 100+id if number==2
drop number
save trt4c, replace

merge m:1 date id state2 using CPSdata_c, keep(match) keepusing(lnepteen lnminwage urall rpteen popteen)
egen exp = group(trtdate id)

gen state = state2 if state2<100
gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

gen popteentrt = popteen if treated==state
sort id date state
by id date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
egen group = group(exp state2)

eststo e2: qui areg lnepteen lnminwage urall rpteen          i.date [w=popteentrt], cluster(group) absorb(group)


use CPSdata, clear
gen mw = exp(lnminwage)
tab mw
sort state date
xtset state date
by state: gen dmw = D.mw
tab dmw
di 4386-3791
keep if dmw>0 & dmw<.
keep date state
save alldmw, replace

rename date trtdate
rename state id
expand 8
sort trtdate id
egen a = group(trtdate id)
spell a
gen time = _seq-4
gen date = trtdate - time
format date %tq
keep trtdate id date
expand 2
sort trtdate id date
by trtdate id date: gen number = _n
gen state2 = id
replace state2 = 100+id if number==2
drop number
save alldmw_sample, replace


* all MW experiments less 90q2
use alldmw_sample, clear
drop if trtdate==yq(1990,2) | trtdate>yq(2010,3)
merge m:1 date id state2 using CPSdata_c, keep(match) keepusing(lnepteen lnminwage urall rpteen popteen)
egen exp = group(trtdate id)

gen state = state2 if state2<100
gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

gen popteentrt = popteen if treated==state
sort id date state
by id date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
egen group = group(exp state2)

eststo e3: qui areg lnepteen lnminwage urall rpteen          i.date [w=popteentrt], cluster(group) absorb(group)


* all MW experiments less 90q2, same-div
use alldmw_sample, clear
drop if trtdate==yq(1990,2) | trtdate>yq(2010,3)
merge m:1 date id state2 using CPSdata_cr, keep(match) keepusing(lnepteen lnminwage urall rpteen popteen)
egen exp = group(trtdate id)

gen state = state2 if state2<100
gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

gen popteentrt = popteen if treated==state
sort id date state
by id date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
egen group = group(exp state2)

eststo e4: qui areg lnepteen lnminwage urall rpteen          i.date [w=popteentrt], cluster(group) absorb(group)



* specify regression to obtain the residuals
use CPSdata, clear
areg lnepteen lnminwage urall rpteen i.date [w=popteen], absorb(state) cluster(state)
predict resida, resid
areg lnepteen           urall rpteen i.date [w=popteen], absorb(state) cluster(state)
predict residb, resid
keep date state resida residb
xtset state date
save resid, replace

use CPSdata, clear
keep state
duplicates drop
save statelist, replace

clear
input null
0
end
save tempa, replace
save tempb, replace

foreach x in a b {
    use if date~=yq(1990,2) & date<=yq(2010,3) using alldmw, clear
    summ
    local total = r(N)
    forval i = 1/`total' {
        use if date~=yq(1990,2) & date<=yq(2010,3) using alldmw, clear
        keep in `i'
        levelsof state, local(treated)
        levelsof date,  local(trtdate)
        local trtdateA = `trtdate' - 4
        local trtdateB = `trtdate' - 3
        local trtdateC = `trtdate' - 2
        local trtdateD = `trtdate' - 1
        merge 1:1 state using statelist, keep(using)
        levelsof state, local(control)
        use resid, clear
        keep if date>=`trtdateA' & date<=`trtdate'
        synth resid`x' resid`x'(`trtdateA') resid`x'(`trtdateB') resid`x'(`trtdateC') resid`x'(`trtdateD'), trperiod(`trtdate') trunit(`treated') counit(`control') xperiod(`trtdateA'/`trtdateD') mspeperiod(`trtdateA'/`trtdateD') customV(1 1 1 1) keep(sc_`trtdate'_`treated') replace
        use sc_`trtdate'_`treated', clear
        rename _Co_Number state
        rename _W_Weight weight
        drop _Y_treated _Y_synthetic _time
        gen trtdate = `trtdate'
        gen trtunit = `treated'
        save, replace
        use temp`x', clear
        append using sc_`trtdate'_`treated'
        save temp`x', replace
        erase sc_`trtdate'_`treated'.dta
    }
    
    use temp`x', clear
    drop in 1
    drop null
    save temp`x', replace
}


use if date~=yq(1990,2) & date<=yq(2010,3) using alldmw, clear
rename  date trtdate
rename state trtunit
cross using statelist
expand 8
sort trtdate trtunit state
by trtdate trtunit state: gen time = _n-4
gen date = trtdate - time
drop time
format date %tq
save listmerge, replace

foreach x in a b {
    use listmerge, clear
    merge m:1 date state using CPSdata, keep(match master) nogen
    merge m:1 trtdate trtunit state using temp`x', keep(match master) nogen
    replace weight = 1 if trtunit==state
    save dup, replace
    
    keep if trtunit==state
    save dup1, replace
    
    * all MW experiments, Abadie weights
    use dup, clear
    drop if trtunit==state
    foreach j in lnepteen lnminwage urall rpteen popteen {
        gen wt`j' = weight*`j'
    }
    collapse (sum) wt*, by(trtunit trtdate date)
    rename wt* *
    save dup0, replace
    append using dup1
    
    * id is treatment state
    * state identifies treatment state, missing for control obs (51 FE)
    * state2 identifies each treatment state and each control state associated with each treatment state (102 FE)
    * state3 identifies each treatment state and each control state associated with each experiment (51 FE + # of experiments)
    
    egen exp = group(trtdate trtunit)
    gen state2 = state
    replace state2 = trtunit+100 if state==.
    
    gen treated = state
    gen control = exp if state==.
    recode treated (.=0)
    recode control (.=0)
    egen state3 = group(treated control)
    
    gen popteentrt = popteen if treated==state
    sort exp date state
    by exp date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
    egen group = group(exp state2)
    
    eststo e5`x': qui areg lnepteen lnminwage urall rpteen          i.date [w=popteentrt], cluster(group) absorb(group)
}

esttab e1 e2 e3 e4 e5a e5b using $folder\output\Table7.csv, keep(lnminwage urall rpteen) compress b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ar2 r2 obslast

/*
* all MW experiments less 90q2, contig state
use alldmw_sample, clear
drop if trtdate==yq(1990,2) | trtdate>yq(2010,3)
merge m:1 date id state2 using CPSdata_b, keep(match) keepusing(lnepteen lnminwage urall rpteen popteen)
keep if state~=2 & state~=15
egen exp = group(trtdate id)

gen state = state2 if state2<100
gen treated = state
gen control = exp if state==.
recode treated (.=0)
recode control (.=0)
egen state3 = group(treated control)

gen popteentrt = popteen if treated==state
sort id date state
by id date: replace popteentrt = popteentrt[_n-1] if popteentrt==.
egen group = group(exp state2)

eststo e6: qui areg lnepteen lnminwage urall rpteen          i.date [w=popteentrt], cluster(group) absorb(group)
*/


