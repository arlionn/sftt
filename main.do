cscript sftt
version 17  // Stata version must be >= 14
pwd
dir

// ----- Load commands -----
adopath + "./src"


// ----- Section 5.1 -----
// The benchmark 2TSF model
// model estimation
sjlog using ./output/output_simu_ori1, replace
clear
set seed 999
quietly set obs 1600
generate x1 = invnormal(uniform())
generate x2 = invnormal(uniform())
generate ue = invexponential(0.6, uniform())
generate we = invexponential(1.4, uniform())
generate v = invnormal(uniform())
generate y = x1 + 2 * x2 - ue + we + v
sftt y x1 x2, nocons
sjlog close, replace
// efficiency analysis
sjlog using ./output/output_simu_ori2, replace
sftt sigs
sftt eff
sum _u_hat_exp _w_hat_exp _wu_diff_exp
sum _wu_diff_exp, detail
sjlog close, replace


// ----- Section 5.2 -----
// 2TSF model with scaling property
// First result - no initial values
sjlog using ./output/output_simu_scal1, replace
clear
set seed 999
quietly set obs 10000
matrix C = (1, 0.1, 0.1 \ 0.1, 1, 0.1 \ 0.1, 0.1, 1)
drawnorm x zu zw, corr(C)
generate ui = invexponential(1, uniform())
generate wi = invexponential(1, uniform())
generate vi = invnormal(uniform())
generate y = x - exp(0.6 * zu) * ui + exp(0.8 * zw) * wi + vi
sftt y x, scal sigmau(zu) sigmaw(zw) robust nocons
sjlog close, replace
// Second result - with initial values
sjlog using ./output/output_simu_scal2, replace
sftt y x, scal sigmau(zu) sigmaw(zw) robust nocons ///
         initial(delta_x 1 du_zu 0.6 mu_u 1 dw_zw 0.8 mu_w 1)
sjlog close, replace


// ----- Section 6.1 -----
// Replicate the results in Kumbhakar and Parmeter (2009)
sjlog using ./output/output_exmp1_1, replace
set seed 20220612
use https://sftt.oss-cn-hangzhou.aliyuncs.com/kp09.dta, clear
sftt lwage iq educ educ2 exper exper2 tenure tenure2 age married south ///
         urban black sibs brthord meduc feduc
sftt sigs
sjlog close, replace
// Second result - efficiency
sjlog using ./output/output_exmp1_2, replace
sftt eff
tabstat _w_hat _u_hat _wu_diff, by(black) stat(mean p25 p50 p75) format(%6.3f) c(s)
tabstat _w_hat_exp _u_hat_exp _wu_diff_exp, by(black) stat(mean p25 p50 p75) ///
        format(%6.3f) c(s)
sjlog close, replace


// ----- Section 6.2 -----
// Replicate the results in Kumbhakar and Parmeter (2010)
sjlog using ./output/output_exmp2_1, replace
use https://sftt.oss-cn-hangzhou.aliyuncs.com/kp10.dta, clear
sftt lprn lsf unitsftc bathstot roomsn sfan sfdn          ///
          agelt5 age510 age1015 agegte30                  ///
          cencityn urbsubn urbann riuraln inadeq degreen  ///
          s87 s88 s89 s90 s91 s92 s93                     ///
          verylg large siz1to3 small,                     ///
     sigmaw(outbuy firstbuy incbuy busbuy agebuy          ///
            blkbuy marbuy sfbuy edubuy kidbuy)            ///
     sigmau(incsell bussell agesell blksell marsell       ///
            sfsell edusell kidsell)                       ///
     hnormal seed(6)
sjlog close, replace


// ----- Section 6.3 -----
// Replicate the results in Lu et al. (2011)
// First result - estimation
sjlog using ./output/output_exmp3_1, replace
set seed 20220612
use https://sftt.oss-cn-hangzhou.aliyuncs.com/lu11.dta, clear
sftt lnprice lnage symp urban education job endurance insur i.province i.year
sjlog close, replace

// Second result - efficiency
sjlog using ./output/output_exmp3_2, replace
sftt eff, exp
sum _u_hat_exp _w_hat_exp _wu_diff_exp
sum _wu_diff_exp, detail
sjlog close, replace

// Third result - histogram
sjlog using ./output/output_exmp3_3, replace
histogram _u_hat_exp, percent title(Percent, place(10) size(*0.7))                 ///
        ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by patients (%)")  ///
        xscale(titlegap(3) outergap(-2))
histogram _w_hat_exp, percent title(Percent, place(10) size(*0.7))                 ///
        ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by doctors (%)")   ///
        xscale(titlegap(3) outergap(-2))
histogram _wu_diff_exp, percent title(Percent, place(10) size(*0.7))               ///
        ylabel(,angle(0)) ytitle("") xtitle("Net Surplus (%)")                    ///
        xscale(titlegap(3) outergap(-2))
sjlog close, replace

// Export the figures into files
histogram _u_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
       ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by patients (%)")  ///
       xscale(titlegap(3) outergap(-2)) scheme(sj)
graph export output/patients.eps, replace
histogram _w_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
       ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by doctors (%)")   ///
       xscale(titlegap(3) outergap(-2)) scheme(sj)
graph export output/doctors.eps, replace
histogram _wu_diff_exp, percent title(Percent, place(10) size(*0.7))             ///
       ylabel(,angle(0)) ytitle("") xtitle("Net Surplus (%)")                    ///
       xscale(titlegap(3) outergap(-2)) scheme(sj)
graph export output/netsurplus.eps, replace

/*
  DISTRIBUTION COMPARISON
  The following code is prepared to generate the results of distribution comparison
  Some format details need manually adjustments to generate the exact results presented in manuscript.
  The original results used to prepare Table 2 & 3 are not sjlogged,
  while the bundled results are esttabed to tex files.
*/
use https://sftt.oss-cn-hangzhou.aliyuncs.com/lu11.dta, clear
quietly {
    // OLS
    reg lnprice lnage symp urban education job endurance insur i.province i.year, r
    est store res0
    // 2TSF - exponential specification
    sftt lnprice lnage symp urban education job endurance insur i.province i.year, seed(20220613)
    est store res1
    noisily sftt sigs
    sftt eff
    foreach var in u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
        rename _`var' _`var'_e
    }
    // 2TSF - half-normal specification
    sftt lnprice lnage symp urban education job endurance insur i.province i.year, hnormal seed(20220613)
    est store res2
    noisily sftt sigs
    sftt eff
}
esttab res0 res1 res2 using ./output/output_empirical_cmp_raw.tex, drop(i_* *.year *.province) ///
        title(Estimation results with different distributions) replace   // Export results into a .tex file

sum _u_hat_e _u_hat
sum _wu_net_effect _wu_net_effect_e
sum _wu_diff_exp_e _wu_diff_exp
sort _u_hat_e
gen rnk_u_e = _n
sort _u_hat
gen rnk_u_n = _n
sort _w_hat_e
gen rnk_w_e = _n
sort _w_hat
gen rnk_w_n = _n
scatter rnk_w_n rnk_w_e, xtitle("Exponential") ytitle("Half-Normal") scheme(sj)
graph export output/wi_rank.eps, replace  // Export figures into .eps files
scatter rnk_u_n rnk_u_e, xtitle("Exponential") ytitle("Half-Normal") scheme(sj)
graph export output/ui_rank.eps, replace
