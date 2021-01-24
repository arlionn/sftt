clear
local simus 10000
local obs 1600
set seed 20200517


// ----- Load commands -----
adopath + "./code"


// ----- settings -----
set obs `simus'
generate beta = .
generate delta_u = .
generate delta_w = .
generate mu_u = .
generate mu_w = .
local non_converge = 0

local delta_u = 0.8
local delta_w = 1.2

foreach obs in 6400  {
	foreach delta_u in 0.6 {
		foreach delta_w in 1.4 {
			forvalues i = 1 / `simus' {
				preserve
					clear
					quietly set obs `obs'
					gen id = mod(_n, 20) + 1
					bysort id: gen t = _n
					matrix C = (1, 0.1, 0.1 \ .1, 1, 0.1 \ 0.1, 0.1, 1)
					drawnorm x1 x2 x3, corr(C)

					generate ui = invexponential(1, uniform())
					generate wi = invexponential(1, uniform())
					generate vi = invnormal(uniform())
					generate y = x1 - exp(`delta_u' * x2) * ui + exp(`delta_w' * x3) * wi + vi

					quietly xtset id t
					capture SFA2tier2 y x1, scal zu(x2) zw(x3) robust nocons iter(100) ///
					        initial(dw_x3 `delta_w' mu_w 1 beta_x1 1 du_x2 `delta_u' mu_u 1)
					matrix b = e(b)
				restore

				if e(converged) {
					quietly {
						replace beta = b[1, 1] if _n == `i'
						replace delta_u = b[1, 2] if _n == `i'
						replace mu_u = b[1, 3] if _n == `i'
						replace delta_w = b[1, 4] if _n == `i'
						replace mu_w = b[1, 5] if _n == `i'
					}
				}
				else {
					quietly display "Not converged: " `++non_converge'
				}
				
				if mod(`i', 100) == 0 {
				    di "`i' tests finished."
				}
				
			}
			display as error "scaling_`delta_u'_`delta_w'_`obs': "  `non_converge' " tests not converged."
			save ./result/scaling_`delta_u'_`delta_w'_`obs'.dta, replace
		}
	}
}





// --- load MC results ---
local obs 1600
local delta_u 0.8
local delta_w 1.2
use	./result/scaling_`delta_u'_`delta_w'_`obs'.dta, clear
// drop unconverged samples
foreach var in beta delta_u delta_w mu_u mu_w {
	drop if `var' == .
}
// drop extrem values
winsor2 beta delta_u delta_w mu_u mu_w, replace cut(0.5 99.5) trim
foreach var in beta delta_u delta_w mu_u mu_w {
	drop if `var' == .
}

generate beta_e = (beta - 1)
generate delta_u_e = (delta_u - `delta_u')
generate delta_w_e = (delta_w - `delta_w')
generate mu_u_e = (mu_u - 1)
generate mu_w_e = (mu_w - 1)

generate beta_b = (beta - 1)^2
generate delta_u_b = (delta_u - `delta_u')^2
generate delta_w_b = (delta_w - `delta_w')^2
generate mu_u_b = (mu_u - 1)^2
generate mu_w_b = (mu_w - 1)^2

sum *_e *_b





// --- export MC results ---
local delta_w 1.4
local output

foreach delta_u in 0.6 0.8 1.1 {
local output
	foreach obs in 400 1600 6400 {
		quietly {
			use	./result/scaling_`delta_u'_`delta_w'_`obs'.dta, clear

			// drop unconverged samples
			foreach var in beta delta_u delta_w mu_u mu_w {
				drop if `var' == .
			}

			// drop extrem values
			winsor2 beta delta_u delta_w mu_u mu_w, replace cut(0.5 99.5) trim
			foreach var in beta delta_u delta_w mu_u mu_w {
				drop if `var' == .
			}

			generate beta_e = (beta - 1)
			generate delta_u_e = (delta_u - `delta_u')
			generate delta_w_e = (delta_w - `delta_w')
			generate mu_u_e = (mu_u - 1)
			generate mu_w_e = (mu_w - 1)

			generate beta_b = (beta - 1)^2
			generate delta_u_b = (delta_u - `delta_u')^2
			generate delta_w_b = (delta_w - `delta_w')^2
			generate mu_u_b = (mu_u - 1)^2
			generate mu_w_b = (mu_w - 1)^2

			foreach var in beta mu_u delta_u mu_w delta_w {
				sum `var'_e
				local tmp = round(r(mean), 0.0001)
				local output `output' `tmp'
			}

			foreach var in beta mu_u delta_u mu_w delta_w {
				sum `var'_b
				local tmp = round(r(mean), 0.0001)
				local output `output' `tmp'
			}
		}
	}
di "`output'"
}













