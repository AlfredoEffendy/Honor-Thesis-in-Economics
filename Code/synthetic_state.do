
*** PREP

global folder "C:\JMIS\UCI\RA-work\MW\ILRR\data"

cd $folder

* for treatment, check previous 4q, 5q, 8q
* for control, check previous 4q, 5q, 8q and future 3q

* code to produce Xq sample
* regardless of dMW in the future length, but should have no dMW in the past length
local y = 3
foreach x in 4 5 8 {

    use CPSdata, clear
    
    preserve
    keep date
    duplicates drop _all, force
    save dates, replace
    restore
    
    xtset state date
    by state: gen dminwage = minwage - minwage[_n-1]
    * only consider changes greater than 5 cents
    gen dmw = float(dminwage)>=.05
    replace dmw = . if dminwage==.
    spell dmw, by(state)
    
    preserve
        xtset state date
        by state: gen _seqminus1   = _seq[_n-1]
        by state: gen _seqminus`x' = _seq[_n-`x']
        * don't care about future MW changes
        keep if dmw==1 & _seq==1 & (_seqminus1 == _seqminus`x' + `x' - 1)
        keep state date
        save treatment`x', replace
    restore
    
    * list of states with continuous quarters of no MW changes
    xtset state date
    by state: gen _seqminus`x' = _seq[_n-`x']
    by state: gen _seqplus`y'  = _seq[_n+`y']
    keep if _seqplus`y'  == _seq + `y'
    keep if _seqminus`x' == _seq - `x'
    keep _end date state
    rename _end status
    save controls`x', replace

}


* specify regression to obtain the residuals
use CPSdata, clear
xtset state date
areg lnepteen lnminwage urall rpteen i.date [w=popteen], absorb(state) cluster(state)
predict resid1, resid
areg lnepteen           urall rpteen i.date [w=popteen], absorb(state) cluster(state)
predict resid2, resid
by state: gen d1lnepteen = lnepteen - lnepteen[_n-1]
by state: gen d4lnepteen = lnepteen - lnepteen[_n-4]
keep date state resid1 resid2 lnepteen d1lnepteen d4lnepteen
save cps_depvar, replace



* get treatments (date and state) and corresponding controls (donor states)
foreach x in 4 5 8 {
    use treatment`x', clear
    rename state trtstate
    joinby date using controls`x'
    rename state ctlstate
    rename date trtdate
    drop status
    sort trtdate trtstate ctlstate
    * drop if there's only one control state
    gen one=1
    by trtdate trtstate: egen total = total(one)
    drop if total==1
    drop one total
    save trtctl`x', replace
}

* get list of counties per state
use qcew2, clear
keep county
duplicates drop
gen state = floor(county/1000)
sort state county
save statecounties, replace

* plug in counties
foreach x in 4 5 8 {
    use trtctl`x', clear
    rename trtstate state
    joinby state using statecounties
    rename state trtstate
    rename county trtcounty
    rename ctlstate state
    joinby state using statecounties
    rename state ctlstate
    rename county ctlcounty
    save trtctl`x'_cty, replace
}

* borderpairs data has county pairings
* note that this file is 2x unique pairings (interchanged places), so no need to do twice

* winnow down to those with cross-border pairs
foreach x in 4 5 8 {
    use trtctl`x'_cty, clear
    merge m:1 trtcounty ctlcounty using borderpairs, keep(match) nogen
    keep trtdate trtcounty
    duplicates drop
    keep if trtdate<=yq(2006,2)-3
    * there should be 121
    count
    save trt`x'_cty, replace

    use trtctl`x'_cty, clear
    merge m:1 trtdate trtcounty using trt`x'_cty, keep(match) nogen
    merge m:1 trtcounty ctlcounty using borderpairs, keep(master match)
    gen contig = _merge==3
    drop _merge
    compress
    save trtctl`x'_cty2, replace
}

* specify regression to obtain the residuals
use qcew2, clear
gen state = floor(county/1000)
xtset county date
areg lneprest lnmw lneppriv i.date if priv==1, absorb(county) cluster(state)
predict resid1, resid
areg lneprest      lneppriv i.date if priv==1, absorb(county) cluster(state)
predict resid2, resid
by county: gen d1lneprest = lneprest - lneprest[_n-1]
by county: gen d4lneprest = lneprest - lneprest[_n-4]
keep date county resid1 resid2 lneprest d1lneprest d4lneprest
save qcew_depvar, replace
















foreach x in 4 5 8 {
    use treatment`x', clear
    append using controls`x'
    recode status (.=1)
    gsort date -status
    * take out dates with no treatment
    by date: egen total = total(status)
    drop if total==0
    drop total
    * take out dates with less than two controls
    replace status = 1-status
    by date: egen total = total(status)
    drop if total<2
    * get list of donors
    replace status = 1-status
    levelsof date, local(trtdates`x')
    foreach i in `trtdates`x'' {
        levelsof state if date==`i' & status==0, local(ctl`x'_`i')
        levelsof state if date==`i' & status==1, local(trt`x'_`i')
    }
}
macro list

clear
input null
0
end
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    save weights_`j', replace
}

* set limits for treatment dates for d1 and d4
* 1990q1 is 120
local start_resid1     124
local start_resid2     124
local start_lnepteen   124
local start_d1lnepteen 125
local start_d4lnepteen 128

* foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
foreach j in                         d1lnepteen d4lnepteen {
    foreach i in `trtdates4' {
        if `i' >= `start_`j'' {
            foreach trt in `trt4_`i'' {
                local trtdate0 = `i'
                local trtdateA = `i' - 1
                local trtdateB = `i' - 2
                local trtdateC = `i' - 3
                local trtdateD = `i' - 4

                use if `j'~=. using cps_depvar, clear
                synth `j' `j'(`trtdateA') `j'(`trtdateB') `j'(`trtdateC') `j'(`trtdateD'), trperiod(`trtdate0') trunit(`trt') counit(`ctl4_`i'') xperiod(`trtdateD'/`trtdateA') mspeperiod(`trtdateD'/`trtdateA') customV(1 1 1 1) keep(sc_`trtdate0'_`trt') replace
                use sc_`trtdate0'_`trt', clear
                rename _Co_Number state
                rename _W_Weight weight
                drop _Y_treated _Y_synthetic _time
                gen trtdate = `trtdate0'
                gen trtunit = `trt'
                save sc_`trtdate0'_`trt', replace
                use weights_`j', clear
                append using sc_`trtdate0'_`trt'
                save weights_`j', replace
                erase sc_`trtdate0'_`trt'.dta
            }
        }
    }

    use weights_`j', clear
    drop in 1
    drop null
    drop if weight==.
    save weights_`j', replace
}

/* only use if want to be strict about no previous MW increase
foreach j in d1lnepteen {
    foreach i in `trtdates5' {
        foreach trt in `trt5_`i'' {
            local trtdate0 = `i'
            local trtdateA = `i' - 1
            local trtdateB = `i' - 2
            local trtdateC = `i' - 3
            local trtdateD = `i' - 4

            use if d1lnepteen~=. using cps_depvar, clear
            synth `j' `j'(`trtdateA') `j'(`trtdateB') `j'(`trtdateC') `j'(`trtdateD'), trperiod(`trtdate0') trunit(`trt') counit(`ctl5_`i'') xperiod(`trtdateD'/`trtdateA') mspeperiod(`trtdateD'/`trtdateA') customV(1 1 1 1) keep(sc_`trtdate0'_`trt') replace
            use sc_`trtdate0'_`trt', clear
            rename _Co_Number state
            rename _W_Weight weight
            drop _Y_treated _Y_synthetic _time
            gen trtdate = `trtdate0'
            gen trtunit = `trt'
            save sc_`trtdate0'_`trt', replace
            use weights_`j', clear
            append using sc_`trtdate0'_`trt'
            save weights_`j', replace
            erase sc_`trtdate0'_`trt'.dta
        }
    }

    use weights_`j', clear
    drop in 1
    drop null
    drop if weight==.
    save weights_`j', replace
}

foreach j in d4lnepteen {
    foreach i in `trtdates8' {
        foreach trt in `trt8_`i'' {
            local trtdate0 = `i'
            local trtdateA = `i' - 1
            local trtdateB = `i' - 2
            local trtdateC = `i' - 3
            local trtdateD = `i' - 4

            use if d4lnepteen~=. using cps_depvar, clear
            synth `j' `j'(`trtdateA') `j'(`trtdateB') `j'(`trtdateC') `j'(`trtdateD'), trperiod(`trtdate0') trunit(`trt') counit(`ctl8_`i'') xperiod(`trtdateD'/`trtdateA') mspeperiod(`trtdateD'/`trtdateA') customV(1 1 1 1) keep(sc_`trtdate0'_`trt') replace
            use sc_`trtdate0'_`trt', clear
            rename _Co_Number state
            rename _W_Weight weight
            drop _Y_treated _Y_synthetic _time
            gen trtdate = `trtdate0'
            gen trtunit = `trt'
            save sc_`trtdate0'_`trt', replace
            use weights_`j', clear
            append using sc_`trtdate0'_`trt'
            save weights_`j', replace
            erase sc_`trtdate0'_`trt'.dta
        }
    }

    use weights_`j', clear
    drop in 1
    drop null
    drop if weight==.
    save weights_`j', replace
}
*/

foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    use weights_`j', clear
    merge m:1 state_fips using censusdivkey, nogen keepusing(censusdiv) keep(master match)
    rename censusdiv ctldiv
    rename state_fips state
    rename trtunit state_fips
    merge m:1 state_fips using censusdivkey, nogen keepusing(censusdiv) keep(master match)
    rename censusdiv trtdiv
    rename state_fips trtunit
    gen samediv = ctldiv==trtdiv
    gen donors = 1
    collapse (sum) weight donors, by(trtdate trtunit samediv)
    gen avewt = weight/donors
    reshape wide weight donors avewt, i(trtdate trtunit) j(samediv)
    drop if donors1==.
    gen donors = donors0 + donors1
    gen higher = avewt1 > avewt0
    sort trtunit trtdate
    gen depvar = "`j'"
    save weights_`j'_one, replace
    preserve
        rename trtunit state_fips
        merge m:1 state_fips using censusdivkey, nogen keepusing(censusdiv) keep(master match)
        gen one = 1
        collapse (mean) weight? avewt? donors donors1 (sum) higher one, by(censusdiv)
        gen depvar = "`j'"
        save weights_`j'_div, replace
    restore
    gen one = 1
    collapse (mean) weight? avewt? donors donors1 (sum) higher one
    gen depvar = "`j'"
    save weights_`j'_all, replace
}

clear
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    append using weights_`j'_one
    erase weights_`j'_one.dta
}
save weights_one, replace
reshape wide weight? avewt? donors donors? higher, i(trtdate trtunit) j(depvar) string
order trtdate trtunit weight1resid1 weight1resid2 weight1lnepteen weight1d1lnepteen weight1d4lnepteen avewt1resid1 avewt1resid2 avewt1lnepteen avewt1d1lnepteen avewt1d4lnepteen avewt0resid1 avewt0resid2 avewt0lnepteen avewt0d1lnepteen avewt0d4lnepteen higherresid1 higherresid2 higherlnepteen higherd1lnepteen higherd4lnepteen donorsresid1 donorsresid2 donorslnepteen donorsd1lnepteen donorsd4lnepteen donors1resid1 donors1resid2 donors1lnepteen donors1d1lnepteen donors1d4lnepteen weight0resid1 weight0resid2 weight0lnepteen weight0d1lnepteen weight0d4lnepteen donors0resid1 donors0resid2 donors0lnepteen donors0d1lnepteen donors0d4lnepteen
sort trtdate trtunit
outsheet using $folder\output\weights_one.csv, comma replace

clear
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    append using weights_`j'_div
    erase weights_`j'_div.dta
}
save weights_div, replace
reshape wide weight? avewt? donors donors? higher one, i(censusdiv) j(depvar) string
order censusdiv weight1resid1 weight1resid2 weight1lnepteen weight1d1lnepteen weight1d4lnepteen avewt1resid1 avewt1resid2 avewt1lnepteen avewt1d1lnepteen avewt1d4lnepteen avewt0resid1 avewt0resid2 avewt0lnepteen avewt0d1lnepteen avewt0d4lnepteen higherresid1 higherresid2 higherlnepteen higherd1lnepteen higherd4lnepteen donorsresid1 donorsresid2 donorslnepteen donorsd1lnepteen donorsd4lnepteen donors1resid1 donors1resid2 donors1lnepteen donors1d1lnepteen donors1d4lnepteen oneresid1 oneresid2 onelnepteen oned1lnepteen oned4lnepteen weight0resid1 weight0resid2 weight0lnepteen weight0d1lnepteen weight0d4lnepteen
sort censusdiv
outsheet using $folder\output\Table3a_weights_div.csv, comma replace

clear
foreach j in resid1 resid2 lnepteen d1lnepteen d4lnepteen {
    append using weights_`j'_all
    erase weights_`j'_all.dta
}
save weights_all, replace
gen x = 0
reshape wide weight? avewt? donors donors? higher one, i(x) j(depvar) string
order x weight1resid1 weight1resid2 weight1lnepteen weight1d1lnepteen weight1d4lnepteen avewt1resid1 avewt1resid2 avewt1lnepteen avewt1d1lnepteen avewt1d4lnepteen avewt0resid1 avewt0resid2 avewt0lnepteen avewt0d1lnepteen avewt0d4lnepteen higherresid1 higherresid2 higherlnepteen higherd1lnepteen higherd4lnepteen donorsresid1 donorsresid2 donorslnepteen donorsd1lnepteen donorsd4lnepteen donors1resid1 donors1resid2 donors1lnepteen donors1d1lnepteen donors1d4lnepteen oneresid1 oneresid2 onelnepteen oned1lnepteen oned4lnepteen weight0resid1 weight0resid2 weight0lnepteen weight0d1lnepteen weight0d4lnepteen
outsheet using $folder\output\Table3b_weights_all.csv, comma replace
