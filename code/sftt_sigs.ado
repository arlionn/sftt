program define sftt_sigs, rclass
    version 13

    preserve
        // recover factor variables
        foreach var in `e(factor_variables)' {
            local tmp_var_name = subinstr("`var'", ".", "_", .)
            fvrevar `var', stub(`tmp_var_name'_)
        }

        // sigs
        if e(sftt_scaling) == "1" {
            sftt_sigs_scal
        }
        else {
            sftt_sigs_ori
        }
    restore
end

program define sftt_sigs_ori, rclass
    version 13
    marksample touse
    markout `touse' `e(depvar)' `e(rhs)' `e(sigmaw)' `e(sigmau)'

    if e(title) == "Two-tier SF Model (2TSF) : HOMO exponential" | ///
       e(title) == "Two-tier SF Model (2TSF) : HOMO half-normal" {
        local _sigv = exp([sigma_v]_cons)
        local _sigu = exp([sigma_u]_cons)
        local _sigw = exp([sigma_w]_cons)
    }
    else {
        tempvar sigu sigw
        local _sigv = exp([sigma_v]_cons)
        quietly predict `sigu' if `touse', eq(sigma_u)
        quietly replace `sigu' = exp(`sigu')
        quietly sum `sigu'
        local _sigu = r(mean)
        quietly predict `sigw' if `touse', eq(sigma_w)
        quietly replace `sigw' = exp(`sigw')
        quietly sum `sigw'
        local _sigw = r(mean)
    }

    local uw_diff   = `_sigu' - `_sigw'
    local sigs_sum  = `_sigv'^2 + `_sigu'^2 + `_sigw'^2
    local sigs_uw_r = (`_sigu'^2 + `_sigw'^2)/`sigs_sum'
    local sigs_u_r  = `_sigu'^2/(`_sigu'^2 + `_sigw'^2)
    local sigs_w_r  = 1 - `sigs_u_r'

    display
    display in g "               Variance Estimation          "
    display in g in smcl "{hline 47}"
    display in g "sigma_v    : " _col(20) in y %6.4f `_sigv'
    display in g "sigma_u    : " _col(20) in y %6.4f `_sigu'
    display in g "sigma_w    : " _col(20) in y %6.4f `_sigw'
    display in g "sigma_v_sq : " _col(20) in y %6.4f `_sigv'^2
    display in g "sigma_u_sq : " _col(20) in y %6.4f `_sigu'^2
    display in g "sigma_w_sq : " _col(20) in y %6.4f `_sigw'^2
    display in g in smcl "{hline 47}"
    display in g "               Variance Analysis          "
    display in g in smcl "{hline 47}"
    display in g "Total sigma_sqs     :  " _col(20) in y %6.4f `sigs_sum'
    display in g "(sigu2+sigw2)/Total :  " _col(20) in y %6.4f `sigs_uw_r'
    display in g "sigu2/(sigu2+sigw2) :  " _col(20) in y %6.4f `sigs_u_r'
    display in g "sigw2/(sigu2+sigw2) :  " _col(20) in y %6.4f `sigs_w_r'
    if `uw_diff' < 0 {
        display in g "sig_u - sig_w       : " _col(20) in y %6.4f `uw_diff'
    }
    else {
        display in g "sig_u - sig_w       :  " _col(20) in y %6.4f `uw_diff'
    }
    display in g in smcl "{hline 47}"

    return scalar sigma_u = `_sigu'
    return scalar sigma_w = `_sigw'
    return scalar sigma_v = `_sigv'
    return scalar sigma_u_sq = `_sigu'^2
    return scalar sigma_w_sq = `_sigw'^2
    return scalar sigma_v_sq = `_sigv'^2
    return scalar totoal_sigma_sq = `sigs_sum'

    return scalar u2_U2plusW2 = `sigs_u_r'
    return scalar w2_U2plusW2 = `sigs_w_r'
    return scalar U2plusW2_total_sq = `sigs_uw_r'
end


program define sftt_sigs_scal, rclass
    version 13
    marksample touse
    markout `touse' `e(depvar)' `e(rhs)' `e(zu)' `e(zw)'

    tempname coe
    matrix `coe' = e(b)

    // calculate zu * delta_u
    tempvar du_zu
    generate double `du_zu' = 0
    foreach var in `e(zu)' {
        local idx = 0
        foreach para in `e(params)' {
            local idx = `idx' + 1
            if "du_`var'" == "`para'" {
                continue, break
            }
        }
        quietly replace `du_zu' = `du_zu' + `var' * `coe'[1, `idx']
    }
    // calculate zw * delta_w
    tempvar dw_zw
    generate double `dw_zw' = 0
    foreach var in `e(zw)' {
        local idx = 0
        foreach para in `e(params)' {
            local idx = `idx' + 1
            if "dw_`var'" == "`para'" {
                continue, break
            }
        }
        quietly replace `dw_zw' = `dw_zw' + `var' * `coe'[1, `idx']
    }
    // find out mu_u & mu_w
    local idx = 0
    foreach para in `e(params)' {
        local idx = `idx' + 1
        if "`para'" == "mu_u" {
            local mu_u = `coe'[1, `idx']
        }
        if "`para'" == "mu_w" {
            local mu_w = `coe'[1, `idx']
        }
    }

    // calculate sigams
    tempvar _sigu _sigw
    generate double `_sigu' = exp(`du_zu') * `mu_u'
    quietly sum `_sigu' if `touse'
    local _sigu = r(mean)

    generate double `_sigw' = exp(`dw_zw') * `mu_w'
    quietly sum `_sigw' if `touse'
    local _sigw = r(mean)

    local _sigv = (`e(rss)' / (`e(N)' - 1))^0.5

    local uw_diff   = `_sigu' - `_sigw'
    local sigs_sum  = `_sigv'^2 + `_sigu'^2 + `_sigw'^2
    local sigs_uw_r = (`_sigu'^2 + `_sigw'^2)/`sigs_sum'
    local sigs_u_r  = `_sigu'^2/(`_sigu'^2 + `_sigw'^2)
    local sigs_w_r  = 1 - `sigs_u_r'

    display
    display in g "               Variance Estimation          "
    display in g in smcl "{hline 47}"
    display in g "sigma_v    : " _col(20) in y %6.4f `_sigv'
    display in g "sigma_u    : " _col(20) in y %6.4f `_sigu'
    display in g "sigma_w    : " _col(20) in y %6.4f `_sigw'
    display in g "sigma_v_sq : " _col(20) in y %6.4f `_sigv'^2
    display in g "sigma_u_sq : " _col(20) in y %6.4f `_sigu'^2
    display in g "sigma_w_sq : " _col(20) in y %6.4f `_sigw'^2
    display in g in smcl "{hline 47}"
    display in g "               Variance Analysis          "
    display in g in smcl "{hline 47}"
    display in g "Total sigma_sqs     :  " _col(20) in y %6.4f `sigs_sum'
    display in g "(sigu2+sigw2)/Total :  " _col(20) in y %6.4f `sigs_uw_r'
    display in g "sigu2/(sigu2+sigw2) :  " _col(20) in y %6.4f `sigs_u_r'
    display in g "sigw2/(sigu2+sigw2) :  " _col(20) in y %6.4f `sigs_w_r'
    if `uw_diff' < 0 {
        display in g "sig_u - sig_w       : " _col(20) in y %6.4f `uw_diff'
    }
    else {
        display in g "sig_u - sig_w       :  " _col(20) in y %6.4f `uw_diff'
    }
    display in g in smcl "{hline 47}"

    return scalar sigma_u = `_sigu'
    return scalar sigma_w = `_sigw'
    return scalar sigma_v = `_sigv'
    return scalar sigma_u_sq = `_sigu'^2
    return scalar sigma_w_sq = `_sigw'^2
    return scalar sigma_v_sq = `_sigv'^2
    return scalar totoal_sigma_sq = `sigs_sum'
    return scalar u2_U2plusW2 = `sigs_u_r'
    return scalar w2_U2plusW2 = `sigs_w_r'
    return scalar U2plusW2_total_sq = `sigs_uw_r'
end
