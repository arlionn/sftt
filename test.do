cd C:\Users\lc_zd\dev\sftt
adopath + "./code"

// ---------------------------
clear
set seed 888
quietly set obs 1600
generate x1 = invnormal(uniform())
generate x2 = invnormal(uniform())
generate ue = invexponential(0.6, uniform())
generate we = invexponential(1.4, uniform())
generate v = invnormal(uniform())
generate y = x1 + 2 * x2 - ue + we + v
sftt y x1 x2, nocons
sftt y x1 x2, nocons hnormal

// ---------------------
clear
set seed 888
quietly set obs 1600
generate x1 = invnormal(uniform())
generate x2 = invnormal(uniform())
generate ue = invexponential(0.6, uniform())
generate we = invexponential(1.4, uniform())
generate ss = invnormal(uniform()) + 1
replace ss = -ss if ss < 0
generate ues = invexponential(ss, uniform())
generate v = invnormal(uniform())
generate ys = x1 + 2 * x2 - ues + we + v
sftt ys x1 x2, nocons sigmau(ss)
sftt ys x1 x2, nocons sigmau(ss) hnormal

// --------------------
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
sftt y x, scal sigmau(zu) sigmaw(zw) robust nocons ///
         initial(delta_x 1 du_zu 0.6 mu_u 1 dw_zw 0.8 mu_w 1)
sftt_sigs

