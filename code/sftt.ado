*! Author: Lian Yu-jun
*! E-mail: arlionn@163.com

program sftt, eclass		
	syntax varlist(min=1 fv) [if] [in] [, sigmau(varlist) sigmaw(varlist) /// 
								  Check SEarch SCALing INITial(string) ///
								  Robust vce(string) /// 
								  noCONStant fe ITERate(integer 1000) ///
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
			
		// process fe
		if "`fe'" == "fe" {
			_xt, trequired				
			quietly xtset
			local id_var `r(panelvar)'
			local t_var `r(timevar)'
			sort `id_var' `t_var'
			foreach var in `normal_varlist' {
				bysort `id_var': egen `var'_b = mean(`var') 
				egen `var'_bb = mean(`var')
				replace `var' = `var' - `var'_b + `var'_bb
				drop `var'_b `var'_bb
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
				forvalues i = 1 / 10 {
					// start estimating
					set seed `i'
					capture {
						if "`scaling'" == "scaling" {
							sftt_scaling `varlist' `options' findingseedmode
							ereturn local sftt_scaling = 1
						}
						else {
							sftt_ori `varlist' `options' findingseedmode
							ereturn local sftt_scaling = 0
						}
					}
					if _rc == 1 {
						di ""
						di "Exiting seed finding."
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
			sftt_scaling `varlist' `options'
			ereturn local sftt_scaling = 1
		}
		else {
			sftt_ori `varlist' `options'
			ereturn local sftt_scaling = 0
		}
		ereturn local factor_variables `ori_factor_variables'
	restore
end


// sftt with scaling property
// capture program drop sftt_scaling
program sftt_scaling, eclass
	syntax varlist(min=2) [if] [in] [, SCALing noCONStant fe INITial(string) ///
			sigmau(varlist) sigmaw(varlist) Robust ITERate(integer 1000) ///
			findseed findingseedmode vce(string)]
	local zu `sigmau'
	local zw `sigmaw'
			
	gettoken y xs: varlist  		
	
	marksample touse, strok
	markout `touse' `y' `xs' `zu' `zw'
	
	// if there is any dependent
	if "`xs'" == "" & "`constant'" == "noconstant" {
		display as error "No independent variable found."
		error 102
	}
	
	// When invoked in finding seed mode, we only iterate 100 times
	if "`findingseedmode'"  == "findingseedmode" {
		local iterate 200
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
	if "`fe'" == "fe" {
		local skipconstant 
		if "`constant'" == "noconstant" {
			local skipconstant skipconstant
		}
		quietly xtset
		nl scaling_opt_fe @ `y' `xs', zu(`zu') zw(`zw') `skipconstant' ///
				id_var(`r(panelvar)') parameters(`para_list')  ///
				iterate(`iterate') initial(`initial') ///
				`robust' vce(`vce') ///
				title(Two-tiered SF Model (2TSF) : Scaling)        		
	}
	else {
		local skipconstant 
		if "`constant'" == "noconstant" {
			local skipconstant skipconstant
		}
		nl scaling_opt @ `y' `xs', zu(`zu') zw(`zw') `skipconstant' ///
				parameters(`para_list') `robust' vce(`vce') iterate(`iterate') ///
				initial(`initial') title(Two-tiered SF Model (2TSF) : Scaling)
	}	
	
	ereturn local zu `zu'
	ereturn local zw `zw'
end

// original sftt
// capture program drop sftt_ori
program define sftt_ori, eclass
version 8.2 
		if replay() {
		
				if "`e(cmd)'" != "SFAsw" {
			error 301
		}
				Replay `0'
				exit
		}

	   syntax varlist(min=1 ts fv) [if] [in] [,noCONStant    ///
			   sigmau(string) sigmaw(string) FE Check SEarch ///
			   Plot Robust vce(string) ITERate(integer 1000) ///
			   HNORMal findseed seed(string) findingseedmode]
		
		if "`constant'" != "" {
			local nocns=", `constant'"
		}  
		if "`robust'" != ""{
			local robust "robust"
		}
		
		// When invoked in finding seed mode, we only iterate 100 times
		if "`findingseedmode'"  == "findingseedmode" {
			local iterate 200
		}
		  
		marksample touse
		markout `touse' `varlist' `sigmau' `sigmaw'
		gettoken lhs varlist: varlist  
		
		if "`varlist'" == "" & "`constant'" != "" {
			 error 102
		}
		   
		tsunab lhs : `lhs'
		/* check `lhs' not constant */
		qui _rmcoll `lhs'
		if "`r(varlist)'" == "" {
			 di as err "dependent variable cannot be constant"
			 exit 198
		}
		
				
		/*
		if "`firmeff'"!="" | "`timeeff'"!=""{
		   qui tsset
		   local id   "`r(panelvar)'"
		   local t    "`r(timevar)'"   
		}
		*/
			
		
		markout `touse' `ivar' `tvar'  /* iis does not allow string */

		qui count if `touse' == 1
		if r(N) == 0 {
			  error 2000
		}   
		


	/*
	** Firm fixed effects and Time fixed effects
			
	  if "`firmeff'"!="" | "`timeeff'"!=""{
		   qui xtdes
		   local NN = r(N)
		   local TT = r(max)
	  }

	  if "`firmeff'" != ""{
		cap drop _dum_fe*
		qui tab `id', gen(_dum_fe)
		drop _dum_fe1
		local varlist `varlist' _dum_fe2-_dum_fe`NN'
		local mlmax "ml max"
		local dropf "drop _dum_fe*"
	  }
	  if "`timeeff'" != ""{
		cap drop _dum_t*
		qui tab `t', gen(_dum_t)
		drop _dum_t1 _dum_t2
		local varlist `varlist' _dum_t3-_dum_t`TT'
		local mlmax "ml max"
		local dropt "drop _dum_t*"
	  }   
	*/
	

	** remove collinearity 
	cap noi _rmdcoll `lhs' `varlist' if `touse'  `nocns'
	if _rc {
		exit _rc
	}
	local varlist `r(varlist)'
	local names `varlist'

	** time-series operator 
	local eq1 : subinstr local lhs "." "_", all
	tsrevar `varlist'
	local varlist `r(varlist)'
	local lhsname `lhs'
	tsrevar `lhs'
	local lhs `r(varlist)'
	markout `touse' `varlist' `lhs'  
	if `"`sigmau'"' == "" & `"`sigmaw'"' == ""{
		if "`hnormal'" != "hnormal" {
			ml model lf sftt_KP09_ll (`eq1': `lhs'=`varlist' `nocns') ///
					(sigma_v: )                                           ///
					(sigma_u: `sigmau')                                   ///
					(sigma_w: `sigmaw' )                                  ///
					if `touse', `robust' vce(`vce') technique(bfgs)       ///
					title(Two-tiered SF Model (2TSF) : HOMO exponential)  
		}
		else {
			ml model lf sftt_Pa15_ll (`eq1': `lhs'=`varlist' `nocns') ///
					 (sigma_v: )                                          ///
					 (sigma_u: `sigmau')                                  ///
					 (sigma_w: `sigmaw' )                                 ///
					 if `touse', `robust' vce(`vce') technique(bfgs)      ///
					 title(Two-tiered SF Model (2TSF) : HOMO half-normal) 
		}
	}
	else {
		if "`hnormal'" != "hnormal" {
			ml model lf sftt_KP09_ll (`eq1': `lhs'=`varlist' `nocns') ///
					 (sigma_v: )                                          ///
					 (sigma_u: `sigmau')                                  ///
					 (sigma_w: `sigmaw' )                                 ///
					 if `touse', `robust' vce(`vce') technique(bfgs)      ///
					 title(Two-tiered SF Model (2TSF) : HET exponential)
		}
		else {
			ml model lf sftt_Pa15_ll (`eq1': `lhs'=`varlist' `nocns') ///
					 (sigma_v: )                                          ///
					 (sigma_u: `sigmau')                                  ///
					 (sigma_w: `sigmaw' )                                 ///
					 if `touse', `robust' vce(`vce') technique(bfgs)      ///
					 title(Two-tiered SF Model (2TSF) : HET half-normal)
		}
	}		
	
// 	if "`check'" != ""{
		quietly ml check
// 	}
// 	if "`search'" != ""{
		quietly ml search   
// 	}
// 	if "`plot'" != ""{             
// 		ml plot _cons 
// 	}       
	ml max, iterate(`iterate') difficult
	//	`dropf'
	//	`dropt'
	ereturn repost
	ereturn local sigmau "`sigmau'"
	ereturn local sigmaw "`sigmaw'"
	ereturn local rhs `varlist'
	ereturn local good_seed "`seed_using'"
end
