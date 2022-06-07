program nl_scaling_opt_fe
    syntax varlist if, at(name) [zu(varlist) zw(varlist) skipconstant id_var(varlist)]
    gettoken y xs: varlist

    // Retrieve parameters out of at matrix
    local para_idx = 1
    // xb
    local expr_xb
    if "`skipconstant'" != "skipconstant" {
        tempname beta0
        scalar `beta0' = `at'[1, `para_idx++']
        local expr_xb `beta0'
    }
    if "`xs'" != "" {
        local expr_xb `expr_xb' +
    }
    local first = 1
    foreach var in `xs' {
        tempname beta_`var'
        scalar `beta_`var'' = `at'[1, `para_idx++']
        if `first' {
            local expr_xb `expr_xb' `beta_`var'' * `var'
            local first = 0
        }
        else {
            local expr_xb `expr_xb' + `beta_`var'' * `var'
        }
    }
    // zu * delta_u
    local first = 1
    foreach var in `zu' {
        tempname du_`var'
        scalar `du_`var'' = `at'[1, `para_idx++']
        if `first' {
            local expr_u `du_`var'' * `var'
            local first = 0
        }
        else {
            local expr_u `expr_u' + `du_`var'' * `var'
        }
    }
    tempname mu_u
    scalar `mu_u' = `at'[1, `para_idx++']
    // zw * delta_w
    local first = 1
    foreach var in `zw' {
        tempname dw_`var'
        scalar `dw_`var'' = `at'[1, `para_idx++']
        if `first' {
            local expr_w `dw_`var'' * `var'
            local first 0
        }
        else {
            local expr_w `expr_w' + `dw_`var'' * `var'
        }
    }
    tempname mu_w
    scalar `mu_w' = `at'[1, `para_idx++']

    tempvar expu expw expu_b expw_b expu_bb expw_bb
    generate double `expu' = exp(`expr_u')
    generate double `expw' = exp(`expr_w')
    bysort `id_var': egen `expu_b' = mean(`expu')
    bysort `id_var': egen `expw_b' = mean(`expw')
    egen `expu_bb' = mean(`expu')
    egen `expw_bb' = mean(`expw')
    
    // Now fill in dependent variable
    replace `y' = `expr_xb' - (`expu' - `expu_b' + `expu_bb') * `mu_u' + ///
            (`expw' - `expw_b' + `expw_bb') * `mu_w'
end
