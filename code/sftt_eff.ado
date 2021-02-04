program define sftt_eff
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
	gen double `ehat' = `e(depvar)' - `yhat' if `touse'

	
	if e(title) == "Two-tier SF Model (2TSF) : HOMO exponential" {
		local sig_u = exp([sigma_u]_cons)
		local sig_w = exp([sigma_w]_cons)
		local sig_v = exp([sigma_v]_cons)
		local lamda  = 1 / `sig_u' + 1 / `sig_w'
		local lamda1 = 1 / `lamda'
		local lamda2 = `lamda' / (1 + `lamda')
	}
	else if e(title) == "Two-tier SF Model (2TSF) : HET exponential" {
		tempvar sig_u sig_w lamda lamda1 lamda2
		local sig_v = exp([sigma_v]_cons)
		qui predict `sig_u' if `touse', eq(sigma_u) 
		qui replace `sig_u' = exp(`sig_u')
		qui predict `sig_w' if `touse', eq(sigma_w)
		qui replace `sig_w' = exp(`sig_w')
		gen double `lamda'  = 1 / `sig_u' + 1 / `sig_w'
		gen double `lamda1' = 1 / `lamda'
		gen double `lamda2' = `lamda' / (1+`lamda')
	}
	else if e(title) == "Two-tier SF Model (2TSF) : HOMO half-normal" {
	    local sig_u = exp([sigma_u]_cons)
		local sig_w = exp([sigma_w]_cons)
		local sig_v = exp([sigma_v]_cons)
	}
	else if e(title) == "Two-tier SF Model (2TSF) : HET half-normal" {
		tempvar sig_u sig_w lamda lamda1 lamda2
		local sig_v = exp([sigma_v]_cons)
		qui predict `sig_u' if `touse', eq(sigma_u) 
		qui replace `sig_u' = exp(`sig_u')
		qui predict `sig_w' if `touse', eq(sigma_w)
		qui replace `sig_w' = exp(`sig_w')
	}
    else {
	    if e(sftt_scaling) == "1" {
			display as error "sftt_eff doesn't support 2TSF model with scaling property."
		}
		else {
			display as error "No qualified estimation result found."
		}
		exit 0
	}
	
	
	if e(title) == "Two-tier SF Model (2TSF) : HOMO exponential" | ///
	   e(title) == "Two-tier SF Model (2TSF) : HET" {
		tempvar aa bb betahat etahat Eta1 Eta2
		quietly generate double `aa'      =  `ehat'/`sig_u' + `sig_v'^2/(2*`sig_u'^2)
		quietly generate double `bb'      =  `ehat'/`sig_v' - `sig_v'/`sig_w'
		quietly generate double `etahat'  =  (0.5*`sig_v'^2)/(`sig_w'^2) - `ehat'/`sig_w'
		quietly generate double `betahat' = - `ehat'/`sig_v' - `sig_v'/`sig_u' 
		quietly generate double `Eta1'    = normal(`bb') + exp(`aa' - `etahat')*normal(`betahat')
		quietly generate double `Eta2'    = exp(`etahat' - `aa')*`Eta1'

		tempvar pp2 pp3
		quietly generate double `pp2' = normalden(-`betahat') + `betahat'*normal(`betahat')
		quietly generate double `pp3' = normalden(-`bb') + `bb'*normal(`bb')

		quietly generate u_hat = `lamda1' + (`sig_v'*`pp2')/`Eta2'            /*E(u_i | epselon_i)*/ 
		quietly generate w_hat = `lamda1' + (`sig_v'*`pp3')/`Eta1'            /*E(w_i | epselon_i)*/
		quietly generate wu_diff = w_hat - u_hat
		
		tempvar qq1 qq2 qq3 qq4 
		quietly generate double `qq1' = `sig_v'^2/2 - `sig_v'*`betahat'
		quietly generate double `qq2' = normal(`bb') + exp(`aa' - `etahat')*exp(`qq1')*normal(`betahat' - `sig_v')
		quietly generate double `qq3' = `sig_v'^2/2 - `sig_v'*`bb'
		quietly generate double `qq4' = normal(`betahat') + exp(`etahat' - `aa')*exp(`qq3')*normal(`bb' - `sig_v')

		quietly generate double u_hat_exp = 1 - `lamda2' * `qq2' / `Eta1'  /*E(1-e^{-u} | epselon)*/
		quietly generate double w_hat_exp = 1 - `lamda2' * `qq4' / `Eta2'  /*E(1-e^{-w} | epselon)*/
		quietly generate double wu_diff_exp = w_hat_exp - u_hat_exp        /*E(e^{-w}-e^{-u}) | epselon*/ 
		// quietly generate double uw_diff_exp2= 1 - (1-u_hat_exp)/(1-w_hat_exp)  /*E[e^{w-u}-1], see KP05b, pp.16*/
		
		// net effect
		tempvar ne1 ne2 ne_p1 ne_p2
		quietly generate double `ne1' = (1 + `sig_u') * (`aa' + `sig_v'^2 / (2 * `sig_u'))
		quietly generate double `ne2' = (1 - `sig_w') * (`etahat' - `sig_v'^2 / (2 * `sig_w'))
		quietly generate double `ne_p1' = exp(`ne1') * normal(`betahat' - `sig_v') + exp(`ne2') * normal(`bb' + `sig_v')
		quietly generate double `ne_p2' = exp(`aa') * normal(`betahat') + exp(`etahat') * normal(`bb')
		quietly generate wu_net_effect = `ne_p1' / `ne_p2' - 1
    }
	else if e(title) == "Two-tier SF Model (2TSF) : HOMO half-normal" {
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
		
		quietly generate u_hat   = `s'^2 * `psi2' - `sig_u'^2 / `s'^2 * (`ehat' - `s'^2 * `psi')
		quietly generate w_hat   = `s'^2 * `psi1' + `sig_w'^2 / `s'^2 * (`ehat' - `s'^2 * `psi')
		quietly generate wu_diff = w_hat - u_hat
		
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
		
		quietly generate u_hat_exp = 1 - 2 / `G_diff' * exp(`tmp_u1') * (normal(`tmp_u2') - binormal(`tmp_u2', `tmp_u3', `rho'))
		quietly generate w_hat_exp = 1 - 2 / `G_diff' * exp(`tmp_w1') * (normal(`tmp_w2') - binormal(`tmp_w2', `tmp_w3', `rho'))		
		quietly generate double wu_diff_exp = w_hat_exp - u_hat_exp        /*E(e^{-w}-e^{-u}) | epselon*/ 
		
		tempvar ne1 ne2 ne3 ne4 ne5 ne_rho1 ne_rho2
		quietly generate double `ne1'     = ((`sig_w'^2 + `sig_u'^2) / `s'^2) * (`sig_v'^2 / 2 + `ehat')
		quietly generate double `ne2'     = (`ehat' + `sig_v'^2) / `omega1'
		quietly generate double `ne3'     = (`ehat' + `sig_v'^2) / `omega2'
		quietly generate double `ne4'     = `ehat' / `omega1'
		quietly generate double `ne5'     = `ehat' / `omega2'
		quietly generate double `ne_rho1' = `lambda1' / sqrt(1 + `lambda1'^2)
		quietly generate double `ne_rho2' = - `lambda2' / sqrt(1 + `lambda2'^2)

		quietly generate wu_net_effect = exp(`ne1') * ///
		                                 (binormal(`ne2', 0, `ne_rho1') - binormal(`ne3', 0, `ne_rho2')) / ///
										 (binormal(`ne4', 0, `ne_rho1') - binormal(`ne5', 0, `ne_rho2')) - 1
	}
end
