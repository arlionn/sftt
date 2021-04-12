
// ----- Load commands -----
cd /Users/rang/sftt
adopath + "./code"


// ----- Section 4.1 -----
// Original 2TSF with simulated data
// First result
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
sftt y x1 x2, check search nocons
sjlog close, replace
// Second result
sjlog using ./output/output_simu_ori2, replace
sftt_sigs
sftt_eff
sum u_hat_exp w_hat_exp wu_diff_exp
sum wu_diff_exp, detail
sjlog close, replace


// ----- Section 4.2 -----
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
// Third result - postestimation
sjlog using ./output/output_simu_scal3, replace
sftt_sigs
sjlog close, replace

// ----- Section 4.3 -----
// I run [scaling_mc.do] and collect the data manually.


// ----- Section 5.1 -----
sjlog using ./output/output_exmp1_1, replace
use kp09, clear
sftt lwage iq educ educ2 exper exper2 tenure tenure2 age married south ///
         urban black sibs brthord meduc feduc, check search seed(999)
sftt_sigs
sjlog close, replace
// Second result - efficiency
sjlog using ./output/output_exmp1_2, replace
sftt_eff
tabstat w_hat u_hat wu_diff, by(black) stat(mean p25 p50 p75) format(%6.3f) c(s)
tabstat w_hat_exp u_hat_exp wu_diff_exp, by(black) stat(mean p25 p50 p75) ///
        format(%6.3f) c(s)
sjlog close, replace


// ----- Section 5.2 -----
// First result - estimation
sjlog using ./output/output_exmp2_1, replace
use lu11, clear
sftt lnprice lnage symp urban education job endurance insur i.province i.year, search check
sjlog close, replace
// Second result - efficiency
sjlog using ./output/output_exmp2_2, replace
sftt_eff
sum u_hat_exp w_hat_exp wu_diff_exp
sum wu_diff_exp, detail
sjlog close, replace
// Third result - histogram
sjlog using ./output/output_exmp2_3, replace
histogram u_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by patients (%)") ///
	   xscale(titlegap(3) outergap(-2))
histogram w_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by doctors (%)")  ///
	   xscale(titlegap(3) outergap(-2))
histogram wu_diff_exp, percent title(Percent, place(10) size(*0.7))             ///
	   ylabel(,angle(0)) ytitle("") xtitle("Net Surplus (%)")                   ///
	   xscale(titlegap(3) outergap(-2))
sjlog close, replace
// Generate the figures used in the paper
histogram u_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by patients (%)") ///
	   xscale(titlegap(3) outergap(-2)) scheme(sj) 
graph export output/patients.eps, replace
histogram w_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by doctors (%)")  ///
	   xscale(titlegap(3) outergap(-2)) scheme(sj) 
graph export output/doctors.eps, replace	   
histogram wu_diff_exp, percent title(Percent, place(10) size(*0.7))             ///
	   ylabel(,angle(0)) ytitle("") xtitle("Net Surplus (%)")                   ///
	   xscale(titlegap(3) outergap(-2)) scheme(sj) 
graph export output/netsurplus.eps, replace	   
drop u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect
quietly {
	// OLS
	reg lnprice lnage symp urban education job endurance insur i.province i.year, r
	est store res0
	// 2TSF - exponential specification
	sftt lnprice lnage symp urban education job endurance insur i.province i.year, search check
	est store res1
	noisily sftt_sigs
	sftt_eff
	foreach var in u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
		rename `var' `var'_e
	}  
	// 2TSF - half-normal specification
	sftt lnprice lnage symp urban education job endurance insur i.province i.year, search check hnormal
	est store res2
	noisily sftt_sigs
	sftt_eff
}
esttab res0 res1 res2, drop(i_* *.year *.province)

sum u_hat_e u_hat
sum wu_net_effect wu_net_effect_e
sum wu_diff_exp_e wu_diff_exp
sort u_hat_e
gen rnk_u_e = _n
sort u_hat
gen rnk_u_n = _n
sort w_hat_e
gen rnk_w_e = _n
sort w_hat
gen rnk_w_n = _n

scatter rnk_w_n rnk_w_e, xtitle("Exponential") ytitle("Half-Normal") scheme(sj) 
graph export output/wi_rank.eps, replace	   
scatter rnk_u_n rnk_u_e, xtitle("Exponential") ytitle("Half-Normal") scheme(sj) 
graph export output/ui_rank.eps, replace	   





