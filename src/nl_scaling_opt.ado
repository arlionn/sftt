program nl_scaling_opt
    version 13
    syntax varlist if, at(name) [zu(varlist) zw(varlist) skipconstant]
    gettoken y xs: varlist

    // Retrieve parameters out of at matrix
    local para_idx = 1
    // xb
    local expr_xb
    if "`skipconstant'" != "skipconstant" {
        tempname beta0
        scalar `beta0' = `at'[1, `para_idx++']
        local expr_xb `beta0' +
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

    // Now fill in dependent variable
    replace `y' = `expr_xb' - exp(`expr_u') * `mu_u' + exp(`expr_w') * `mu_w'
end
