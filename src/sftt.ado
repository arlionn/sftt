*! sftt 1.0.0 06Dec2023

/*
SFTT MAIN FUNCTION
Invoke the appropriate subfunction according to the subcommand.

Related command:
    sftt

Syntax:
    - sftt version
        sftt, version
        sftt, ver
    - sftt regression
        sftt y x1 x2, nocons
        sftt y x1 x2 x3, hnormal
        sftt y x, scal sigmau(zu) sigmaw(zw) robust
    - sftt variance decomposition
        sftt sigs
    - sftt efficiency analysis
        sftt eff           // all measures of inefficiency, default names
        sftt eff, level    // only level specification of inefficiencies, default names
        sftt eff, exp      // only logarithmic specification of inefficiencies
        sftt eff, relative // only relative measures of inefficiencies
        sftt, level u_hat(my_u_hat)  // assign variable name
*/
program sftt
    version 13
    gettoken before_comma after_comma: 0, parse(,)
    gettoken subcmd after_subcmd: before_comma

    if strtrim("`after_comma'") == "version" | strtrim("`after_comma'") == "ver"  {
        display _n(1) "{bf:{ul:Version}}"
        display _n(1) "{txt}{sf}    sftt version 1.0.0"
        display "    December 06, 2023"
        exit 0
    }

    if "`0'" == "sigs" {
        _sftt_sigs
    }
    else {
        if "`subcmd'" == "eff" {
            _sftt_eff `after_subcmd' `after_comma'
        }
        else {
            _sftt_regression `0'
        }
    }
end


/*
SFTT REGRESSION
Estimate sftt models, invoked without subcommand.

Related commands:
    _sftt_regression
    _sftt_regression_scaling
    _sftt_regression_original
*/
program _sftt_regression, eclass
    version 13
    syntax varlist(min=1 fv) [if] [in] [, sigmau(varlist) sigmaw(varlist)    ///
                                  SCALing INITial(string) Robust vce(string) ///
                                  noCONStant ITERate(integer 1000)        ///
                                  HNORMal findseed seed(string)]
    preserve
        // Clear ereturn macro
        ereturn clear
        // classify variables
        local normal_varlist
        local factor_varlist
        local ori_factor_variables
        foreach var in `varlist' {
            local tmp_var_name = subinstr("`var'", ".", "_", .)
            fvrevar `var', stub(`tmp_var_name'_)
            if "`var'" == "`r(varlist)'" {
                local normal_varlist `normal_varlist' `var'
            }
            else {
                local factor_varlist `factor_varlist' `r(varlist)'
                local ori_factor_variables `ori_factor_variables' `var'
            }
        }
        // get command options
        gettoken first options: 0, parse(",")
        // recover variable list
        local varlist `normal_varlist' `factor_varlist'
        // Find seed
        if "`findseed'" == "findseed" {
            if "`seed'" != "" {
                display as error "Random seed provided, skipping seed finding. Remove seed() to use this option."
                local seed_using `seed'
            }
            else {
                display "Finding seeds ", _continue
                forvalues i = 1 / 100 {
                    // start estimating
                    set seed `i'
                    capture {
                        if "`scaling'" == "scaling" {
                            _sftt_regression_scaling `varlist' `if' `options' findingseedmode
                            ereturn local sftt_scaling = 1
                        }
                        else {
                            _sftt_regression_original `varlist' `if' `options' findingseedmode
                            ereturn local sftt_scaling = 0
                        }
                    }
                    if _rc == 1 {
                        display ""
                        display "Exiting seed finding."
                        continue, break
                    }
                    if "`e(converged)'" == "1" {
                        display "."
                        display "Using seed `i'."
                        local seed_using `i'
                        set seed `seed_using'
                        continue, break
                    }
                    display ".", _continue
                }
            }
        }
        if "`seed'" != "" {
            set seed `seed'
        }
        // start estimating
        if "`scaling'" == "scaling" {
            _sftt_regression_scaling `varlist' `if' `options'
            ereturn local sftt_scaling = 1
        }
        else {
            _sftt_regression_original `varlist' `if' `options'
            ereturn local sftt_scaling = 0
            if "`sigmau'" == "" & "`sigmaw'" == "" {
                ereturn local sftt_model "homogeneous"
            }
            else {
                ereturn local sftt_model "heterogeneous"
            }
            if "`hnormal'" == "hnormal" {
                ereturn local sftt_specification "half-normal"
            }
            else {
                ereturn local sftt_specification "exponential"
            }
        }
        ereturn local factor_variables `ori_factor_variables'
    restore
end


// sftt with scaling property
program _sftt_regression_scaling, eclass
    version 13
    syntax varlist(min=2) [if] [in] [, SCALing noCONStant INITial(string) ///
            sigmau(varlist) sigmaw(varlist) Robust ITERate(integer 1000)     ///
            findseed seed(string) findingseedmode vce(string)]
    local zu `sigmau'
    local zw `sigmaw'
    gettoken y xs: varlist
    marksample touse, strok
    markout `touse' `y' `xs' `zu' `zw'

    // if there is any dependent
    if "`xs'" == "" & "`constant'" == "noconstant" {
        display as error "No independent variable found."
        exit 102
    }
    // When invoked in finding seed mode, we only iterate 100 times
    if "`findingseedmode'"  == "findingseedmode" {
        local iterate 100
    }
    // if y is constant
    quietly _rmcoll `y'
    if r(k_omitted) {
        display as error "Dependent variable cannot be constant."
        exit 198
    }
    // remove collinearity in linear model
    if "`constant'" == "noconstant" {
        capture noisily _rmdcoll `y' `xs' if `touse', nocons
    }
    else {
        capture noisily _rmdcoll `y' `xs' if `touse'
    }
    if _rc == 459 {
        exit _rc
    }
    local xs
    local linear_omit
    foreach var in `r(varlist)' {
        if substr("`var'", 1, 2) == "o." {
            local linear_omit = "`linear_omit' " + substr("`var'", 3, .)
            continue
        }
        local xs `xs' `var'
    }
    // generate default initial value with reg / xtreg if not provided
    if "`initial'" == "" {
        if "`constant'" == "noconstant" {
            quietly reg `y' `xs', nocons
        }
        else {
            quietly reg `y' `xs'
        }
        matrix b = e(b)
        local coe_idx = 1
        foreach var in `xs' {
            local tmp_coe = b[1, `coe_idx++']
            local initial `initial' delta_`var' `tmp_coe'
        }
        if "`constant'" != "noconstant" {
            local tmp_coe = b[1, `coe_idx++']
            local initial delta0 `tmp_coe' `initial'
        }
        foreach var in `zu' {
            local initial `initial' du_`var' 1
        }
        local initial `initial' mu_u 1
        foreach var in `zw' {
            local initial `initial' dw_`var' 1
        }
        local initial `initial' mu_w 1

    }
    display "initial value: `initial'"
    // generate parameters list
    local para_list
    if "`constant'" != "noconstant" {
        local para_list delta0
    }
    foreach var in `xs' {
        local para_list `para_list' delta_`var'
    }
    foreach var in `zu' {
        local para_list `para_list' du_`var'
    }
    local para_list `para_list' mu_u
    foreach var in `zw' {
        local para_list `para_list' dw_`var'
    }
    local para_list `para_list' mu_w
    // nls estimation
    local skipconstant
    if "`constant'" == "noconstant" {
        local skipconstant skipconstant
    }
    capture noisily nl _scaling_opt @ `y' `xs' if `touse', zu(`zu') zw(`zw') `skipconstant'  ///
            `robust' vce(`vce') iterate(`iterate')                                ///
            parameters(`para_list') initial(`initial')                            ///
            title("{bf:Two-tier stochastic frontier model with scaling property}")
    if _rc {
        display as error "Convergence not achieved, you can adjust the initial values, or remove {bf:scaling} to estimate with distributional assumptions. See {help sftt}."
        exit _rc
    }
    ereturn local zu `zu'
    ereturn local zw `zw'
end


// original sftt without scaling property
program define _sftt_regression_original, eclass
    version 8.2
    syntax varlist(min=1 ts fv) [if] [in] [,noCONStant      ///
            sigmau(string) sigmaw(string) Check SEarch   ///
            Plot Robust vce(string) ITERate(integer 1000)   ///
            HNORMal findseed seed(string) findingseedmode]

    if "`constant'" != "" {
        local nocns=", `constant'"
    }
    if "`robust'" != ""{
        local robust "robust"
    }
    // When invoked in finding seed mode, we only iterate 100 times
    if "`findingseedmode'"  == "findingseedmode" {
        local iterate 100
    }
    marksample touse
    markout `touse' `varlist' `sigmau' `sigmaw'
    gettoken lhs varlist: varlist
    if "`varlist'" == "" & "`constant'" != "" {
        error 102
    }

    tsunab lhs : `lhs'
    /* check `lhs' not constant */
    quietly _rmcoll `lhs'
    if "`r(varlist)'" == "" {
        display as err "dependent variable cannot be constant"
        exit 198
    }

    markout `touse' `ivar' `tvar'  /* iis does not allow string */
    quietly count if `touse' == 1
    if r(N) == 0 {
        error 2000
    }

    // remove collinearity
    capture noi _rmdcoll `lhs' `varlist' if `touse'  `nocns'
    if _rc {
        exit _rc
    }
    local varlist `r(varlist)'
    local names `varlist'

    // time-series operator
    local eq1 : subinstr local lhs "." "_", all
    tsrevar `varlist'
    local varlist `r(varlist)'
    local lhsname `lhs'
    tsrevar `lhs'
    local lhs `r(varlist)'
    markout `touse' `varlist' `lhs'
    if "`hnormal'" != "hnormal" {
        ml model lf _sftt_kp09_ll (frontier_`eq1': `lhs'=`varlist' `nocns')    ///
                (ln_sig_v: )                                                   ///
                (ln_sig_u: `sigmau')                                           ///
                (ln_sig_w: `sigmaw')                                           ///
                if `touse', `robust' vce(`vce')                                ///
                title("{bf:Two-tier stochastic frontier model with exponential specification}")
    }
    else {
        ml model lf _sftt_pa15_ll (frontier_`eq1': `lhs'=`varlist' `nocns')    ///
                (ln_sig_v: )                                               ///
                (ln_sig_u: `sigmau')                                       ///
                (ln_sig_w: `sigmaw')                                       ///
                if `touse', `robust' vce(`vce')                            ///
                title("{bf:Two-tier stochastic frontier model with half-normal specification}")
    }
    quietly ml check
    quietly ml search
    capture noisily ml max, nooutput iterate(`iterate') difficult
    if _rc {
        display as error "Convergence not achieved, you can set another randon seed using the option {bf:seed(#)}, or, specify {bf:findseed} option. See {help sftt}."
        exit _rc
    }
    version 13
    ml display
    ereturn repost
    ereturn local cmd "sftt"
    ereturn local cmdbase "ml"
    ereturn local sigmau "`sigmau'"
    ereturn local sigmaw "`sigmaw'"
    ereturn local rhs `varlist'
    ereturn local good_seed "`seed_using'"
    estimates title: Model 2TSF
    estimates store Model2TSF
end


/*
SFTT SIGS
Calculate standard deviation and decompose the composite error term.

Related command:
    _sftt_sigs
    _sftt_sigs_original
    _sftt_sigs_scaling
*/
program define _sftt_sigs, rclass
    version 13
    if "`e(sftt_scaling)'" == "" {
        display as error "No available {bf:sftt} estimation result found."
        exit 301
    }
    preserve
        // recover factor variables
        foreach var in `e(factor_variables)' {
            local tmp_var_name = subinstr("`var'", ".", "_", .)
            fvrevar `var', stub(`tmp_var_name'_)
        }

        // sigs
        if e(sftt_scaling) == "1" {
            _sftt_sigs_scaling
        }
        else {
            _sftt_sigs_original
        }
    restore
end


program define _sftt_sigs_original, rclass
    version 13
    marksample touse
    markout `touse' `e(depvar)' `e(rhs)' `e(sigmaw)' `e(sigmau)'

    if e(sftt_model) == "homogeneous" {
        local _sigv = exp([ln_sig_v]_cons)
        local _sigu = exp([ln_sig_u]_cons)
        local _sigw = exp([ln_sig_w]_cons)
    }
    else {
        tempvar sigu sigw
        local _sigv = exp([ln_sig_v]_cons)
        quietly predict `sigu' if `touse', eq(ln_sig_u)
        quietly replace `sigu' = exp(`sigu')
        quietly sum `sigu'
        local _sigu = r(mean)
        quietly predict `sigw' if `touse', eq(ln_sig_w)
        quietly replace `sigw' = exp(`sigw')
        quietly sum `sigw'
        local _sigw = r(mean)
    }

    local wu_diff   = `_sigw' - `_sigu'
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
    if `wu_diff' < 0 {
        display in g "sig_w - sig_u       : " _col(20) in y %6.4f `wu_diff'
    }
    else {
        display in g "sig_w - sig_u       :  " _col(20) in y %6.4f `wu_diff'
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


program define _sftt_sigs_scaling, rclass
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

    local wu_diff   = `_sigw' - `_sigu'
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
    if `wu_diff' < 0 {
        display in g "sig_w - sig_u       : " _col(20) in y %6.4f `wu_diff'
    }
    else {
        display in g "sig_w - sig_u       :  " _col(20) in y %6.4f `wu_diff'
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


/*
SFTT EFF
Efficiency analysis.
This command is not available for 2TSF model with scaling property.

Related commands:
    _sftt_eff
*/
program define _sftt_eff
    version 13
    syntax [if] [, LEVel exp ABSolute RELative replace                   ///
              u_hat(string) w_hat(string) wu_diff(string)  ///
              u_hat_exp(string) w_hat_exp(string)          ///
              wu_diff_exp(string) wu_net_effect(string)]
    if "`e(sftt_scaling)'" == "" {
        display as error "No available {bf:sftt} estimation result found."
        exit 301
    }
    if e(sftt_scaling) == "1" {
        display as error "{bf:sftt eff} is not available for 2TSF model with scaling property."
        exit 499
    }

    // syntax check
    if "`level'" == "level" & "`exp'" == "exp" {
        display as error "Option {bf:level} and {bf:exp} should not be used simultaneously."
        exit 198
    }
    if "`absolute'" == "absolute" & "`relative'" == "relative" {
        display as error "Option {bf:absolute} and {bf:relative} should not be used simultaneously."
        exit 198
    }
    // assign varname
    foreach var in u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
        if "``var''" == "" {
            local `var' = "_`var'"
        }
    }
    if "`level'" == "level" {
        foreach var in u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
            local `var' = ""
        }
    }
    if "`exp'" == "exp" {
        foreach var in u_hat w_hat wu_diff {
            local `var' = ""
        }
    }

    // check if variable already defined
    foreach var in u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
        if "``var''" != "" {
            quietly capture confirm variable ``var'', exact
            if _rc == 0 {
                if "`replace'" == "replace" {
                    drop ``var''
                }
                else {
                    display as error "Variable {bf:``var''} already defined, please assign another name or use {bf:replace} to overwrite."
                    exit 110
                }
            }
        }
    }

    // recover factor variables
    foreach var in `e(factor_variables)' {
        local tmp_var_name = subinstr("`var'", ".", "_", .)
        fvrevar `var', stub(`tmp_var_name'_)
    }

    // mark sample
    marksample touse
    markout `touse' `e(depvar)' `e(rhs)' `e(sigmaw)' `e(sigmau)'
    tempvar yhat ehat
    predict double `yhat' if `touse'
    generate double `ehat' = `e(depvar)' - `yhat' if `touse'

    local sftt_eff_model_located 0
    if e(sftt_model) == "homogeneous" & e(sftt_specification) == "exponential" {
        local sig_u = exp([ln_sig_u]_cons)
        local sig_w = exp([ln_sig_w]_cons)
        local sig_v = exp([ln_sig_v]_cons)
        local lamda  = 1 / `sig_u' + 1 / `sig_w'
        local lamda1 = 1 / `lamda'
        local lamda2 = `lamda' / (1 + `lamda')
        local sftt_eff_model_located 1
    }
    if e(sftt_model) == "heterogeneous" & e(sftt_specification) == "exponential" {
        tempvar sig_u sig_w lamda lamda1 lamda2
        local sig_v = exp([ln_sig_v]_cons)
        quietly predict `sig_u' if `touse', eq(ln_sig_u)
        quietly replace `sig_u' = exp(`sig_u')
        quietly predict `sig_w' if `touse', eq(ln_sig_w)
        quietly replace `sig_w' = exp(`sig_w')
        generate double `lamda'  = 1 / `sig_u' + 1 / `sig_w'
        generate double `lamda1' = 1 / `lamda'
        generate double `lamda2' = `lamda' / (1+`lamda')
        local sftt_eff_model_located 1
    }
    if e(sftt_model) == "homogeneous" & e(sftt_specification) == "half-normal" {
        local sig_u = exp([ln_sig_u]_cons)
        local sig_w = exp([ln_sig_w]_cons)
        local sig_v = exp([ln_sig_v]_cons)
        local sftt_eff_model_located 1
    }
    if e(sftt_model) == "heterogeneous" & e(sftt_specification) == "half-normal" {
        tempvar sig_u sig_w lamda lamda1 lamda2
        local sig_v = exp([ln_sig_v]_cons)
        quietly predict `sig_u' if `touse', eq(ln_sig_u)
        quietly replace `sig_u' = exp(`sig_u')
        quietly predict `sig_w' if `touse', eq(ln_sig_w)
        quietly replace `sig_w' = exp(`sig_w')
        local sftt_eff_model_located 1
    }
    if `sftt_eff_model_located' == 0 {
        display as error "No available {bf:sftt} estimation result found."
        exit 301
    }

    if e(sftt_specification) == "exponential" {
        tempvar aa bb betahat etahat Eta1 Eta2
        quietly generate double `aa'      = `ehat' / `sig_u' + `sig_v'^2 / (2 * `sig_u'^2)
        quietly generate double `bb'      = `ehat' / `sig_v' - `sig_v' / `sig_w'
        quietly generate double `etahat'  = (0.5 * `sig_v'^2) / (`sig_w'^2) - `ehat' / `sig_w'
        quietly generate double `betahat' = -`ehat' / `sig_v' - `sig_v' / `sig_u'
        quietly generate double `Eta1'    = normal(`bb') + exp(`aa' - `etahat') * normal(`betahat')
        quietly generate double `Eta2'    = exp(`etahat' - `aa') * `Eta1'

        tempvar pp2 pp3
        quietly generate double `pp2' = normalden(-`betahat') + `betahat' * normal(`betahat')
        quietly generate double `pp3' = normalden(-`bb') + `bb' * normal(`bb')
        if "`u_hat'" != "" {
            quietly generate double `u_hat' = `lamda1' + (`sig_v' * `pp2') / `Eta2' `if'     /* E(u_i | epselon_i) */
        }
        if "`w_hat'" != "" {
            quietly generate double `w_hat' = `lamda1' + (`sig_v' * `pp3') / `Eta1' `if'     /* E(w_i | epselon_i) */
        }
        if "`wu_diff'" != "" {
            quietly generate double `wu_diff' = `w_hat' - `u_hat' `if'
        }

        tempvar qq1 qq2 qq3 qq4
        quietly generate double `qq1' = `sig_v'^2/2 - `sig_v'*`betahat'
        quietly generate double `qq2' = normal(`bb') + exp(`aa' - `etahat') * exp(`qq1') * normal(`betahat' - `sig_v')
        quietly generate double `qq3' = (`sig_v'^2) / 2 - `sig_v' * `bb'
        quietly generate double `qq4' = normal(`betahat') + exp(`etahat' - `aa') * exp(`qq3') * normal(`bb' - `sig_v')
        if "`u_hat_exp'" != "" {
            di "quietly generate double `u_hat_exp' = 1 - `lamda2' * `qq2' / `Eta1' `if'"
            quietly generate double `u_hat_exp' = 1 - `lamda2' * `qq2' / `Eta1' `if'     /* E(1-e^{-u} | epselon) */
        }
        if "`w_hat_exp'" != "" {
            quietly generate double `w_hat_exp' = 1 - `lamda2' * `qq4' / `Eta2' `if'     /* E(1-e^{-w} | epselon) */
        }
        if "`wu_diff_exp'" != "" {
            quietly generate double `wu_diff_exp' = `w_hat_exp' - `u_hat_exp' `if'       /* E(e^{-w}-e^{-u}) | epselon */
        }

        // net effect
        tempvar ne1 ne2 ne_p1 ne_p2
        quietly generate double `ne1' = (1 + `sig_u') * (`aa' + `sig_v'^2 / (2 * `sig_u'))
        quietly generate double `ne2' = (1 - `sig_w') * (`etahat' - `sig_v'^2 / (2 * `sig_w'))
        quietly generate double `ne_p1' = exp(`ne1') * normal(`betahat' - `sig_v') + exp(`ne2') * normal(`bb' + `sig_v')
        quietly generate double `ne_p2' = exp(`aa') * normal(`betahat') + exp(`etahat') * normal(`bb')
        if "`wu_net_effect'" != "" {
            quietly generate double `wu_net_effect' = `ne_p1' / `ne_p2' - 1 `if'
        }
    }
    else if e(sftt_specification) == "half-normal" {
        tempvar theta1 theta2 tmp_sqrt s
        quietly generate double `theta1'   = `sig_w' / `sig_v'
        quietly generate double `theta2'   = `sig_u' / `sig_v'
        quietly generate double `tmp_sqrt' = sqrt(1 + `theta1'^2 + `theta2'^2)
        quietly generate double `s'        = `sig_v' * `tmp_sqrt'

        tempvar omega1 omega2 lambda1 lambda2
        quietly generate double `omega1'   = `s' * sqrt(1 + `theta2'^2) / `theta1'
        quietly generate double `omega2'   = `s' * sqrt(1 + `theta1'^2) / `theta2'
        quietly generate double `lambda1'  = `theta2' / `theta1' * `tmp_sqrt'
        quietly generate double `lambda2'  = `theta1' / `theta2' * `tmp_sqrt'
        tempvar g1 g2 rho1 rho2 G1 G2
        quietly generate double `g1'   = 2 / `omega1' * normalden(`ehat' / `omega1') * normal(-`lambda1' * `ehat' / `omega1')
        quietly generate double `g2'   = 2 / `omega2' * normalden(`ehat' / `omega2') * normal(`lambda2' * `ehat' / `omega2')
        quietly generate double `rho1' = `lambda1' / sqrt(1 + `lambda1'^2)
        quietly generate double `rho2' = - `lambda2' / sqrt(1 + `lambda2'^2)
        quietly generate double `G1'   = 2 * binormal(`ehat' / `omega1', 0, `rho1')
        quietly generate double `G2'   = 2 * binormal(`ehat' / `omega2', 0, `rho2')

        tempvar psi1 psi2 psi
        quietly generate double `psi1' = `g1' / (`G1' - `G2')
        quietly generate double `psi2' = `g2' / (`G1' - `G2')
        quietly generate double `psi'  = `psi1' - `psi2'
        if "`u_hat'" != "" {
            quietly generate double `u_hat'   = `s'^2 * `psi2' - `sig_u'^2 / `s'^2 * (`ehat' - `s'^2 * `psi') `if'
        }
        if "`w_hat'" != "" {
            quietly generate double `w_hat'   = `s'^2 * `psi1' + `sig_w'^2 / `s'^2 * (`ehat' - `s'^2 * `psi') `if'
        }
        if "`wu_diff'" != "" {
            quietly generate double `wu_diff' = `w_hat' - `u_hat' `if'
        }

        tempvar s1 s2 omega_w omega_u
        quietly generate double `s1'      = sqrt(`sig_w'^2 + `sig_v'^2)
        quietly generate double `s2'      = sqrt(`sig_u'^2 + `sig_v'^2)
        quietly generate double `omega_u' = `sig_u' * `s1' / `s'
        quietly generate double `omega_w' = `sig_w' * `s2' / `s'

        tempvar G_diff rho
        quietly generate double `G_diff' = `G1' - `G2'
        quietly generate double `rho'    = - `sig_w' * `sig_u' / (`s1' * `s2')

        tempvar tmp_w1 tmp_w2 tmp_w3 tmp_u1 tmp_u2 tmp_u3
        quietly generate double `tmp_w1' = `omega_w'^2 / 2 - `omega_w' / `omega1' * `ehat'
        quietly generate double `tmp_w2' = - (`ehat' + `sig_w'^2) / `omega2'
        quietly generate double `tmp_w3' = `omega_w' - `ehat' / `omega1'
        quietly generate double `tmp_u1' = `omega_u'^2 / 2 + `omega_u' / `omega2' * `ehat'
        quietly generate double `tmp_u2' = (`ehat' - `sig_u'^2) / `omega1'
        quietly generate double `tmp_u3' = `omega_u' + `ehat' / `omega2'
        if "`u_hat_exp'" != "" {
            quietly generate double `u_hat_exp' = 1 - 2 / `G_diff' * exp(`tmp_u1') * (normal(`tmp_u2') - ///
                                                  binormal(`tmp_u2', `tmp_u3', `rho')) `if'
        }
        if "`w_hat_exp'" != "" {
            quietly generate double `w_hat_exp' = 1 - 2 / `G_diff' * exp(`tmp_w1') * (normal(`tmp_w2') - ///
                                                  binormal(`tmp_w2', `tmp_w3', `rho')) `if'
        }
        if "`wu_diff_exp'" != "" {
            quietly generate double `wu_diff_exp' = `w_hat_exp' - `u_hat_exp' `if'  /*E(e^{-w}-e^{-u}) | epselon*/
        }

        tempvar ne1 ne2 ne3 ne4 ne5 ne_rho1 ne_rho2
        quietly generate double `ne1'     = ((`sig_w'^2 + `sig_u'^2) / `s'^2) * (`sig_v'^2 / 2 + `ehat')
        quietly generate double `ne2'     = (`ehat' + `sig_v'^2) / `omega1'
        quietly generate double `ne3'     = (`ehat' + `sig_v'^2) / `omega2'
        quietly generate double `ne4'     = `ehat' / `omega1'
        quietly generate double `ne5'     = `ehat' / `omega2'
        quietly generate double `ne_rho1' = `lambda1' / sqrt(1 + `lambda1'^2)
        quietly generate double `ne_rho2' = - `lambda2' / sqrt(1 + `lambda2'^2)
        if "`wu_net_effect'" != "" {
            quietly generate `wu_net_effect' = exp(`ne1') * ///
                                              (binormal(`ne2', 0, `ne_rho1') - binormal(`ne3', 0, `ne_rho2')) /   ///
                                              (binormal(`ne4', 0, `ne_rho1') - binormal(`ne5', 0, `ne_rho2')) - 1 ///
                                              `if'
        }
    }

    if "`relative'" == "relative" {
        foreach var in u_hat w_hat u_hat_exp w_hat_exp {
            if "``var''" != "" {
                drop ``var''
                local `var' = ""
            }
        }
    }
    if "`absolute'" == "absolute" {
        foreach var in wu_diff wu_diff_exp wu_net_effect {
            if "``var''" != "" {
                drop ``var''
                local `var' = ""
            }
        }
    }
    display as text "The following variables have been generated:"
    foreach var in u_hat w_hat wu_diff u_hat_exp w_hat_exp wu_diff_exp wu_net_effect {
        if "``var''" != "" {
            display as result "``var''"
        }
    }
end
