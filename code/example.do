// ----- Section 4.1 -----
// Original 2TSF with simulated data
// First result
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
// Second result
sftt_sigs
sftt_eff
sum u_hat_exp w_hat_exp wu_diff_exp
sum wu_diff_exp, detail


// ----- Section 4.2 -----
// 2TSF model with scaling property
// without initial values
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
// with initial values
sftt y x, scal sigmau(zu) sigmaw(zw) robust nocons ///
         initial(delta_x 1 du_zu 0.6 mu_u 1 dw_zw 0.8 mu_w 1)
// postestimation
sftt_sigs


// ----- Section 5.1 -----
// reproduce result in KP09
use kp09, clear
sftt lwage iq educ educ2 exper exper2 tenure tenure2 age married south ///
         urban black sibs brthord meduc feduc, seed(999)
sftt lwage iq educ educ2 exper exper2 tenure tenure2 age married south ///
         urban black sibs brthord meduc feduc, seed(666) hnorm
sftt_sigs
sftt_eff
tabstat w_hat u_hat wu_diff, by(black) stat(mean p25 p50 p75) format(%6.3f) c(s)
tabstat w_hat_exp u_hat_exp wu_diff_exp, by(black) stat(mean p25 p50 p75) ///
        format(%6.3f) c(s)



// ----- Section 5.2 -----
// estimation
use lu11, clear
sftt lnprice lnage symp urban education job endurance insur i.province i.year, 
// inefficiency estimation
sftt_eff
sum u_hat_exp w_hat_exp wu_diff_exp
sum wu_diff_exp, detail
// histogram
histogram u_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by patients (%)") ///
	   xscale(titlegap(3) outergap(-2))
histogram w_hat_exp, percent title(Percent, place(10) size(*0.7))               ///
	   ylabel(,angle(0)) ytitle("") xtitle("Surplus extracted by doctors (%)")  ///
	   xscale(titlegap(3) outergap(-2))
histogram wu_diff_exp, percent title(Percent, place(10) size(*0.7))             ///
	   ylabel(,angle(0)) ytitle("") xtitle("Net Surplus (%)")                   ///
	   xscale(titlegap(3) outergap(-2))



