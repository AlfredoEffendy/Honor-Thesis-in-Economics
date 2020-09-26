
global folder "C:\JMIS\UCI\RA-work\MW\ILRR\data"

cd $folder


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
    by trtdate trtcounty: gen rank_`j' = _n
    keep trtdate trtcounty ctlcounty contig mspe_`j' rank_`j'
    gen depvar = "`j'"
    save trtctl_`j'_50, replace
restore
}

* change to trtctl5_cty2 if want to expand selection window without MW changes
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
    by trtdate trtcounty: gen rank_`j' = _n
    gen depvar = "`j'"
    save trtctl_`j'_50, replace
}

* change to trtctl8_cty2 if want to expand selection window without MW changes
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
    by trtdate trtcounty: gen rank_`j' = _n
    gen depvar = "`j'"
    save trtctl_`j'_50, replace
}


foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    use trtctl_`j'_50, clear
    * take out 9013 and 17001 only for resid
    drop if mspe_`j'==.
    keep if rank_`j'<=50 | contig==1
    keep trtdate trtcounty ctlcounty depvar
    save trtctl_`j'_50, replace
}

clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    append using trtctl_`j'_50
    erase trtctl_`j'_50.dta
}
save trtctl_50, replace

use trtctl_50, clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    qui levelsof trtdate if depvar=="`j'", local(trtdates_`j')
    foreach i in `trtdates_`j'' {
        qui levelsof trtcounty if depvar=="`j'" & trtdate==`i', local(trt_`j'_`i')
        foreach k in `trt_`j'_`i'' {
            qui levelsof ctlcounty if depvar=="`j'" & trtdate==`i' & trtcounty==`k', local(ctl_`j'_`i'_`k') separate(,)
        }
    }

    file open  myfile using top50county_`j'.R, write replace
    file write myfile `"library(foreign)"' _n
    file write myfile `"library(LowRankQP)"' _n
    file write myfile `"library(optimx)"' _n
    file write myfile `"library(rgenoud)"' _n
    file write myfile `"library(Synth)"' _n _n
    file write myfile `"# set working directory"' _n
    file write myfile `"setwd("c:/JMIS/UCI/RA-work/MW/synth2/cty")"' _n
    file write myfile `"V=c(1,1,1,1)"' _n _n _n _n
    file close myfile

    foreach i in `trtdates_`j'' {
        local ctlp4 = `i' - 4
        local ctlp3 = `i' - 3
        local ctlp2 = `i' - 2
        local ctlp1 = `i' - 1
        foreach k in `trt_`j'_`i'' {

            file open  myfile using top50county_`j'.R, write append
            file write myfile `"# ## top50 donor pool for county `k', date `i'"' _n
            file write myfile `"tests.stata <- read.dta("qcew_depvar.dta")"' _n
            file write myfile `"dataprep.out <- dataprep("' _n
            file write myfile _tab `"foo=tests.stata,"' _n
            file write myfile _tab `"unit.variable="county","' _n
            file write myfile _tab `"time.variable="date","' _n
            file write myfile _tab `"dependent="`j'","' _n
            file write myfile _tab `"treatment.identifier=`k',"' _n
            file write myfile _tab `"controls.identifier=c(`ctl_`j'_`i'_`k''),"' _n
            file write myfile _tab `"time.predictors.prior=c(`ctlp4':`ctlp1'),"' _n
            file write myfile _tab `"time.optimize.ssr=c(`ctlp4':`ctlp1'),"' _n
            file write myfile _tab `"special.predictors=list("' _n
            file write myfile _tab _tab `"list("`j'",`ctlp4',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp3',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp2',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp1',"mean"))"' _n
            file write myfile `")"' _n
            file write myfile `"# synth.out <- synth(data.prep.obj=dataprep.out,verbose=TRUE)"' _n
            file write myfile `"# write.csv(synth.out\$solution.w, "county_`j'_`i'_`k'_top50.csv")"' _n
            file write myfile `"synth.out <- synth(data.prep.obj=dataprep.out,custom.v=V,verbose=TRUE)"' _n
            file write myfile `"write.csv(synth.out\$solution.w, "county_`j'_`i'_`k'_top50_vi.csv")"' _n _n _n _n
            file close myfile
*"
        }
    }
}

* *******************************************************
* run produced files in R!!!
* *******************************************************
* afterwards, continue below
* *******************************************************

use trtctl_50, clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    qui levelsof trtdate if depvar=="`j'", local(trtdates_`j')
    foreach i in `trtdates_`j'' {
        qui levelsof trtcounty if depvar=="`j'" & trtdate==`i', local(trt_`j'_`i')
    }
}

* create folder ...
mkdir $folder\county
* ... and then move new CSV files inside it

cd $folder\county

* convert Excel weights file to Stata
clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    foreach i in `trtdates_`j'' {
        foreach k in `trt_`j'_`i'' {
            insheet using county_`j'_`i'_`k'_top50_vi.csv, comma clear
            rename v1 county
            rename wweight weight
            gen depvar = "`j'"
            gen trtdate = `i'
            gen trtunit = `k'
            compress
            save wtcounty_`j'_`i'_`k'_top50_vi, replace
        }
    }
}

foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    clear
    foreach i in `trtdates_`j'' {
        foreach k in `trt_`j'_`i'' {
            append using wtcounty_`j'_`i'_`k'_top50_vi
            erase wtcounty_`j'_`i'_`k'_top50_vi.dta
            erase county_`j'_`i'_`k'_top50_vi.csv
        }
    }
    save wtcounty_top50_`j', replace
}

clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    append using wtcounty_top50_`j'
    erase wtcounty_top50_`j'.dta
}
save wtcounty_top50_depvar, replace

use wtcounty_top50_depvar, clear
rename county ctlcounty
rename trtunit trtcounty
merge m:1 trtcounty ctlcounty using borderpairs, keep(master match)
gen contig = _merge==3
drop _merge
* to fend off cases of no pairing
egen a = total(contig), by(depvar trtdate trtcounty)
drop if a==0
gen donors = 1
collapse (sum) weight donors, by(depvar trtdate trtcounty contig)
gen avewt = weight/donors
reshape wide weight donors avewt, i(depvar trtdate trtcounty) j(contig)
gen donors = donors0 + donors1
gen higher = avewt1 > avewt0 if avewt1<.
sort depvar trtdate trtcounty
save wtcounty_top50_qcew_one, replace
gen one = 1
collapse (mean) weight? avewt? donors donors1 (sum) higher one, by(depvar)
save wtcounty_top50_qcew_all, replace

use wtcounty_top50_qcew_one, replace
reshape wide weight? avewt? donors donors? higher, i(trtdate trtcounty) j(depvar) string
forval i = 1/6 {
    gen a`i' = .
}
order trtdate trtcounty weight1resid1 weight1resid2 weight1lneprest weight1d1lneprest weight1d4lneprest a1 avewt1resid1 avewt1resid2 avewt1lneprest avewt1d1lneprest avewt1d4lneprest a2 avewt0resid1 avewt0resid2 avewt0lneprest avewt0d1lneprest avewt0d4lneprest a3 higherresid1 higherresid2 higherlneprest higherd1lneprest higherd4lneprest a4 donorsresid1 donorsresid2 donorslneprest donorsd1lneprest donorsd4lneprest a5 donors1resid1 donors1resid2 donors1lneprest donors1d1lneprest donors1d4lneprest a6 weight0resid1 weight0resid2 weight0lneprest weight0d1lneprest weight0d4lneprest donors0resid1 donors0resid2 donors0lneprest donors0d1lneprest donors0d4lneprest
sort trtdate trtcounty
outsheet using $folder\output\Table5a_wtcounty_top50_qcew_one.csv, comma replace

use wtcounty_top50_qcew_all, clear
gen x = 0
reshape wide weight? avewt? donors donors? higher one, i(x) j(depvar) string
forval i = 1/7 {
    gen a`i' = .
}
order x weight1resid1 weight1resid2 weight1lneprest weight1d1lneprest weight1d4lneprest a1 avewt1resid1 avewt1resid2 avewt1lneprest avewt1d1lneprest avewt1d4lneprest a2 avewt0resid1 avewt0resid2 avewt0lneprest avewt0d1lneprest avewt0d4lneprest a3 higherresid1 higherresid2 higherlneprest higherd1lneprest higherd4lneprest a4 donorsresid1 donorsresid2 donorslneprest donorsd1lneprest donorsd4lneprest a5 donors1resid1 donors1resid2 donors1lneprest donors1d1lneprest donors1d4lneprest a6 oneresid1 oneresid2 onelneprest oned1lneprest oned4lneprest a7 weight0resid1 weight0resid2 weight0lneprest weight0d1lneprest weight0d4lneprest
outsheet using $folder\output\Table5a_wtcounty_top50_qcew_all.csv, comma replace







cd $folder

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


use trtctl_all, clear
gen ctlstate = floor(ctlcounty/1000)
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    qui levelsof trtdate if depvar=="`j'", local(trtdates_`j')
    foreach i in `trtdates_`j'' {
        qui levelsof trtcounty if depvar=="`j'" & trtdate==`i', local(trt_`j'_`i')
        foreach k in `trt_`j'_`i'' {
            qui levelsof ctlstate if depvar=="`j'" & trtdate==`i' & trtcounty==`k', local(ctl_`j'_`i'_`k')
            foreach m in `ctl_`j'_`i'_`k'' {
                qui levelsof ctlcounty if depvar=="`j'" & trtdate==`i' & trtcounty==`k' & ctlstate==`m', local(ctl_`j'_`i'_`k'_`m') separate(,)
            }
        }
    }

    file open  myfile using fullcounty_`j'.R, write replace
    file write myfile `"library(foreign)"' _n
    file write myfile `"library(LowRankQP)"' _n
    file write myfile `"library(optimx)"' _n
    file write myfile `"library(rgenoud)"' _n
    file write myfile `"library(Synth)"' _n _n
    file write myfile `"# set working directory"' _n
    file write myfile `"setwd("c:/JMIS/UCI/RA-work/MW/synth2/cty")"' _n
    file write myfile `"V=c(1,1,1,1)"' _n _n _n _n
    file close myfile

    foreach i in `trtdates_`j'' {
        local ctlp4 = `i' - 4
        local ctlp3 = `i' - 3
        local ctlp2 = `i' - 2
        local ctlp1 = `i' - 1
        foreach k in `trt_`j'_`i'' {
            local M ""
            foreach m in `ctl_`j'_`i'_`k'' {
                local M "`M',_`ctl_`j'_`i'_`k'_`m''"
            }

            file open  myfile using fullcounty_`j'.R, write append
            file write myfile `"# ## full donor pool for county `k', date `i'"' _n
            file write myfile `"tests.stata <- read.dta("qcew_depvar.dta")"' _n
            file write myfile `"dataprep.out <- dataprep("' _n
            file write myfile _tab `"foo=tests.stata,"' _n
            file write myfile _tab `"unit.variable="county","' _n
            file write myfile _tab `"time.variable="date","' _n
            file write myfile _tab `"dependent="`j'","' _n
            file write myfile _tab `"treatment.identifier=`k',"' _n
            file write myfile _tab `"controls.identifier=c(`M'"' _n
            file write myfile _tab `"),"' _n
            file write myfile _tab `"time.predictors.prior=c(`ctlp4':`ctlp1'),"' _n
            file write myfile _tab `"time.optimize.ssr=c(`ctlp4':`ctlp1'),"' _n
            file write myfile _tab `"special.predictors=list("' _n
            file write myfile _tab _tab `"list("`j'",`ctlp4',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp3',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp2',"mean"),"' _n
            file write myfile _tab _tab `"list("`j'",`ctlp1',"mean"))"' _n
            file write myfile `")"' _n
            file write myfile `"# synth.out <- synth(data.prep.obj=dataprep.out,verbose=TRUE)"' _n
            file write myfile `"# write.csv(synth.out\$solution.w, "county_`j'_`i'_`k'_full.csv")"' _n
            file write myfile `"synth.out <- synth(data.prep.obj=dataprep.out,custom.v=V,verbose=TRUE)"' _n
            file write myfile `"write.csv(synth.out\$solution.w, "county_`j'_`i'_`k'_full_vi.csv")"' _n _n _n _n
            file close myfile
*"
        }
    }
}

* ******************************************
* using extended search mode,
* find and replace: (,_ -> (\r\n\t
*                    ,_ -> ,\r\n\t
* ******************************************

* run produced files in R!!!

* create folder ...
mkdir $folder\state
* ... and then move new CSV files inside it




use $folder\state\trtctl_all, clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    qui levelsof trtdate if depvar=="`j'", local(trtdates_`j')
    foreach i in `trtdates_`j'' {
        qui levelsof trtcounty if depvar=="`j'" & trtdate==`i', local(trt_`j'_`i')
    }
}

cd $folder\state

* convert Excel weights file to Stata
clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    foreach i in `trtdates_`j'' {
        foreach k in `trt_`j'_`i'' {
            insheet using county_`j'_`i'_`k'_full_vi.csv, comma clear
            rename v1 county
            rename wweight weight
            gen depvar = "`j'"
            gen trtdate = `i'
            gen trtunit = `k'
            compress
            save wtcounty_`j'_`i'_`k'_full_vi, replace
        }
    }
}

foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    clear
    foreach i in `trtdates_`j'' {
        foreach k in `trt_`j'_`i'' {
            append using wtcounty_`j'_`i'_`k'_full_vi
            erase wtcounty_`j'_`i'_`k'_full_vi.dta
            erase county_`j'_`i'_`k'_full_vi.csv
        }
    }
    save wtcounty_`j', replace
}

clear
foreach j in resid1 resid2 lneprest d1lneprest d4lneprest {
    append using wtcounty_`j'
    erase wtcounty_`j'.dta
}
save wtcounty_depvar, replace

use wtcounty_depvar, clear
rename county ctlcounty
rename trtunit trtcounty
merge m:1 trtcounty ctlcounty using borderpairs, keep(master match)
gen contig = _merge==3
drop _merge
* to fend off cases of no pairing
egen a = total(contig), by(depvar trtdate trtcounty)
drop if a==0
gen donors = 1
collapse (sum) weight donors, by(depvar trtdate trtcounty contig)
gen avewt = weight/donors
reshape wide weight donors avewt, i(depvar trtdate trtcounty) j(contig)
gen donors = donors0 + donors1
gen higher = avewt1 > avewt0 if avewt1<.
sort depvar trtdate trtcounty
save wtcounty_qcew_one, replace
gen one = 1
collapse (mean) weight? avewt? donors donors1 (sum) higher one, by(depvar)
save wtcounty_qcew_all, replace

use wtcounty_qcew_one, clear
reshape wide weight? avewt? donors donors? higher, i(trtdate trtcounty) j(depvar) string
forval i = 1/6 {
    gen a`i' = .
}
order trtdate trtcounty weight1resid1 weight1resid2 weight1lneprest weight1d1lneprest weight1d4lneprest a1 avewt1resid1 avewt1resid2 avewt1lneprest avewt1d1lneprest avewt1d4lneprest a2 avewt0resid1 avewt0resid2 avewt0lneprest avewt0d1lneprest avewt0d4lneprest a3 higherresid1 higherresid2 higherlneprest higherd1lneprest higherd4lneprest a4 donorsresid1 donorsresid2 donorslneprest donorsd1lneprest donorsd4lneprest a5 donors1resid1 donors1resid2 donors1lneprest donors1d1lneprest donors1d4lneprest a6 weight0resid1 weight0resid2 weight0lneprest weight0d1lneprest weight0d4lneprest donors0resid1 donors0resid2 donors0lneprest donors0d1lneprest donors0d4lneprest
sort trtdate trtcounty
outsheet using $folder\output\Table5b_wtcounty_qcew_one.csv, comma replace
* create 

use wtcounty_qcew_all, clear
gen x = 0
reshape wide weight? avewt? donors donors? higher one, i(x) j(depvar) string
forval i = 1/7 {
    gen a`i' = .
}
order x weight1resid1 weight1resid2 weight1lneprest weight1d1lneprest weight1d4lneprest a1 avewt1resid1 avewt1resid2 avewt1lneprest avewt1d1lneprest avewt1d4lneprest a2 avewt0resid1 avewt0resid2 avewt0lneprest avewt0d1lneprest avewt0d4lneprest a3 higherresid1 higherresid2 higherlneprest higherd1lneprest higherd4lneprest a4 donorsresid1 donorsresid2 donorslneprest donorsd1lneprest donorsd4lneprest a5 donors1resid1 donors1resid2 donors1lneprest donors1d1lneprest donors1d4lneprest a6 oneresid1 oneresid2 onelneprest oned1lneprest oned4lneprest a7 weight0resid1 weight0resid2 weight0lneprest weight0d1lneprest weight0d4lneprest
outsheet using $folder\output\Table5b_wtcounty_qcew_all.csv, comma replace



