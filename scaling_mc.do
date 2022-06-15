// ----- Load commands -----
adopath + "./code"


// ----- settings -----
clear
local simus 10000
set seed 20220613
set obs `simus'
generate beta = .
generate delta_u = .
generate delta_w = .
generate mu_u = .
generate mu_w = .
local non_converge = 0

foreach obs in 100 400 1600 6400 {
    foreach delta_u in 0.6 0.8 1.1 {
        foreach delta_w in 0.5 1.2 1.4 {
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
                    capture sftt y x1, scal sigmau(x2) sigmaw(x3) robust nocons iter(100) ///
                            initial(dw_x3 `delta_w' mu_w 1 delta_x1 1 du_x2 `delta_u' mu_u 1)
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
            local test_name "scaling_`delta_u'_`delta_w'_`obs'_`simus'"
            display as error "`test_name': " `non_converge' " test(s) not converged."
            save ./mc_results/`test_name'.dta, replace
        }
    }
}


// --- export MC results ---
local simus 10000
foreach delta_w in 0.5 1.2 1.4 {
    foreach delta_u in 0.6 0.8 1.1 {
        foreach obs in 100 400 1600 6400 {
            quietly {
                noisily di as result "delta_u: `delta_u', delta_w: `delta_w', obs: `obs'"
                use	./mc_results/scaling_`delta_u'_`delta_w'_`obs'_`simus'.dta, clear

                // drop unconverged samples
                foreach var in beta delta_u delta_w mu_u mu_w {
                    drop if `var' == .
                }

                // drop extrem values
                winsor2 beta delta_u delta_w mu_u mu_w, replace cut(0.5 99.5) trim
                foreach var in beta delta_u delta_w mu_u mu_w {
                    drop if `var' == .
                }

                // bias
                generate beta_bias = beta - 1
                generate delta_u_bias = delta_u - `delta_u'
                generate delta_w_bias = delta_w - `delta_w'
                generate mu_u_bias = mu_u - 1
                generate mu_w_bias = mu_w - 1
                // mean square error
                generate beta_mse = (beta - 1)^2
                generate delta_u_mse = (delta_u - `delta_u')^2
                generate delta_w_mse = (delta_w - `delta_w')^2
                generate mu_u_mse = (mu_u - 1)^2
                generate mu_w_mse = (mu_w - 1)^2

                local output
                foreach var in beta mu_u delta_u mu_w delta_w {
                    sum `var'_bias
                    local tmp = round(r(mean), 0.0001)
                    local output "`output' & `tmp'"
                }
                noi di as text "BIAS: `output'"
                local output
                foreach var in beta mu_u delta_u mu_w delta_w {
                    sum `var'_mse
                    local tmp = round(r(mean), 0.0001)
                    local output "`output' & (`tmp')"
                }
                noi di as text " MSE: `output'"
                noi di ""
            }
        }
    }
}
