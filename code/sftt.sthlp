{smcl}
{* *! version 1.1.2  25nov2011}{...}
{cmd:help sftt}{right:also see:  {help sftt_sigs}, {help sftt_eff}}
{hline}

{title:Title}

{p2colset 5 18 23 2}{...}
{p2col :{hi:sftt} {hline 2}}Two-tier stochastic frontier model{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}
{cmd:sftt}
{depvar}
[{indepvars}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 33 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Estimation}
{synopt :{opt nocons:tant}}suppress constant term{p_end}
{synopt :{opt hnorm:al}}uses the half-normal specification rather than the benchmark exponential specification{p_end}
{synopt :{opt scal:ing}}estimates the 2TSF model with scaling property by NLS{p_end}

{syntab :Ancillary variables}
{synopt :{cmdab:sigmau(}{it:{help varlist:varlist_u}}{cmd:)}}specify explanatory
variables for the lower (nonpositive) inefficiency variance function{p_end}
{synopt :{cmdab:sigmaw(}{it:{help varlist:varlist_w}}{cmd:)}}specify explanatory
variables for the upper (nonnegative) inefficiency variance function{p_end}

{syntab :Other options}
{synopt :{opt seed(#)}}sets a random seed before estimating to ensure that the results are reproducible{p_end}
{synopt :{opt findseed}}loops through 100 estimations, during which the random seed was set from 1 to 100. In each attampt, we iterate at most 200 times{p_end}
{synopt :{opt iter:ate(#)}}specifies the maximum number of iterations, default is {opt iterate(1000)}. In most cases, the optimization should be converged in less than 1000 iterations.{p_end}
{synopt :{opt fe}}uses the fix-effect estimator. In {opt fe} mode, noconstant should not be added to absorb the constant term{p_end}

{syntab :SE}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt oim}, {opt opg}, {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt boot:strap}, or {opt jack:knife}{p_end}
{synopt :{opt r:obust}}synonym for {cmd:vce(robust)}{p_end}

{p2colreset}{...}
{p 4 6 2}
{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}
Note that when using {opt hnormal}, the estimation might not converge because
of flat derivatives or missing values. 
Rerun the command might help, you can also use {opt findseed} to find a usable random seed. 
Or just removing this option to use the exponential specification{p_end}
{p 4 6 2}


{title:Description}

{pstd}
{opt sftt} fits two-tier stochastic frontier (2TSF) models; the default
is a 2TSF model with inefficiency terms assumed to be distributed exponential. 
It provides estimators for the parameters of
a linear model with a disturbance that is assumed to be a mixture of three components, 
which have a strictly nonnegative, nonpositive and symmetric distribution, respectively.
{opt sftt} can fit models in which the nonnegative / nonpositive distributed component
(a measurement of inefficiency) is assumed to be from a half-normal or exponential, 
or even without distributional assumption. In the latter case, maximization is performed
through nonlinear least-squares estimation.

{title:Options}

{dlgtab:Estimation}

{phang}
{opt noconstant}; see
{helpb estimation options##noconstant:[R] estimation options}.

{phang}
{opt hnormal} uses the half-normal specification rather than the benchmark exponential specification. Note that when using this option, the estimation might not converge because
of flat derivatives or missing values. Rerun the command might help, you can also use {opt findseed} to find a usable random seed. Or just removing this option to use the exponential specification.

{phang}
{opt scaling} estimates the 2TSF model with scaling property by NLS, 
distributional assumption for the inefficiency terms are not needed. 

{dlgtab:Ancillary equations}

{phang}
{cmd:sigmau(}{help varlist:varlist_u} [,{opt noconstant}]{cmd:)}
specifies that the lower (nonpositive) inefficiency component is heteroskedastic,
with the variance expressed as a function of the covariates defined in 
{it:varlist_u}.

{phang}
{cmd:sigmaw(}{help varlist:varlist_u} [,{opt noconstant}]{cmd:)}
specifies that the lower (nonnegative) inefficiency component is heteroskedastic,
with the variance expressed as a function of the covariates defined in 
{it:varlist_u}.

{dlgtab:Other options}

{phang}
{opt seed(#)} sets a random seed before estimating to ensure that the results are reproducible.

{phang}
{opt findseed} loops through 100 estimations, during which the random seed was set from 1 to 100. In each attampt, we iterate at most 200 times.

{phang}
{opt iterate(#)} specifies the maximum number of iterations, default is {opt iterate(1000)}. In most cases, the optimization should be converged in less than 1000 iterations.

{dlgtab:SE}

{phang}
{opt vce(vcetype)} specifies the type of standard error reported,
which includes types that are derived from asymptotic theory and
that use bootstrap or jackknife methods; see 
{helpb vce_option:[R] {it:vce_option}}.

{phang}
{opt robust} synonym for {cmd:vce(robust)}

{title:Examples}

    {hline}
    Setup
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 996}{p_end}
{phang2}{cmd:. set obs 1600}{p_end}
{phang2}{cmd:. generate x1 = invnormal(uniform())}{p_end}
{phang2}{cmd:. generate x2 = invnormal(uniform())}{p_end}
{phang2}{cmd:. generate ue = invexponential(0.6, uniform())}{p_end}
{phang2}{cmd:. generate we = invexponential(1.4, uniform())}{p_end}
{phang2}{cmd:. generate v = invnormal(uniform())}{p_end}
{phang2}{cmd:. generate generate y = x1 + 2 * x2 - ue + we + v}{p_end}

{pstd}Two-tier stochastic frontier model with exponential distribution for
inefficiency terms{p_end}
{phang2}{cmd:. sftt y x1 x2, nocons}{p_end}

{pstd}Two-tier stochastic frontier model with half-normal distribution for
inefficiency terms{p_end}
{phang2}{cmd:. sftt y x1 x2, nocons hnormal}{p_end}

{pstd}Use {cmd: findseed} to find a usable seed{p_end}
{phang2}{cmd:. sftt y x1 x2, nocons hnormal findseed}{p_end}

{pstd}Identify the variances of inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. sftt_sigs}{p_end}

{pstd}Decompose the residual into inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. sftt_eff}{p_end}


    {hline}
    Setup
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 996}{p_end}
{phang2}{cmd:. set obs 10000}{p_end}
{phang2}{cmd:. matrix C = (1, 0.1, 0.1 \ 0.1, 1, 0.1 \ 0.1, 0.1, 1)}{p_end}
{phang2}{cmd:. drawnorm x zu zw, corr(C)}{p_end}
{phang2}{cmd:. generate ui = invexponential(1, uniform())}{p_end}
{phang2}{cmd:. generate wi = invexponential(1, uniform())}{p_end}
{phang2}{cmd:. generate vi = invnormal(uniform())}{p_end}
{phang2}{cmd:. generate y = x - exp(0.6 * zu) * ui + exp(0.8 * zw) * wi + vi}{p_end}

{pstd}Two-tier stochastic frontier model with scaling assumption for
inefficiency terms{p_end}
{phang2}{cmd:. sftt y x, scal sigmau(zu) sigmaw(zw) robust}{p_end}

{pstd}Identify the variances of inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. sftt_sigs}{p_end}
    {hline}


{title:Saved results}

{pstd}
{cmd:sftt} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test; usually stored{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(mss)}}model sum of squares{p_end}
{synopt:{cmd:e(mms)}}model mean square{p_end}
{synopt:{cmd:e(msr)}}residual mean square{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(dev)}}residual deviance{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(df_t)}}total degrees of freedom{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(cj)}}position of constant in {cmd:e(b)} or {cmd:0} if no constant{p_end}
{synopt:{cmd:e(gm_2)}}square of geometric mean of (y-k) if {cmd:lnlsq},
{cmd:1} otherwise{p_end}
{synopt:{cmd:e(delta)}}relative change used to compute derivatives{p_end}
{synopt:{cmd:e(log_t)}}{cmd:1} if {cmd:lnlsq} specified, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(sftt_scaling)}}{cmd:1} if using {cmd:scaling}, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(rhs)}}names of independent variables; always stored{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ml} or {cmd:nl}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}; usually stored{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared test; usually stored{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(opt)}}type of optimization; always stored{p_end}
{synopt:{cmd:e(title)}}title in estimation output; usually stored by commands using {cmd:ml} or {cmd:nl}{p_end}
{synopt:{cmd:e(depvar)}}names of dependent variables; always stored{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method; always stored by commands using {cmd:ml}{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(zu)}}explanatory variables for the lower inefficiency variance function{p_end}
{synopt:{cmd:e(zw)}}explanatory variables for the upper inefficiency variance function{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(marginsok)}}predictions allowed by {cmd:margins}{p_end}
{synopt:{cmd:e(marginsnotok)}}predictions disallowed by {cmd:margins}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(funcprog)}}function evaluator program{p_end}
{synopt:{cmd:e(type)}}{cmd:1} = interactively entered expression{p_end}
{synopt:}{cmd:2} = substitutable expression program{p_end}
{synopt:}{cmd:3} = function evaluator program{p_end}
{synopt:{cmd:e(cmdline)}}options as typed{p_end}
{synopt:{cmd:e(params)}}names of parameters{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(V_modelbased)}}model-based variance{p_end}

{title:Authors}

{pstd}Yujun Lian{p_end}
{pstd}Lingnan College, Sun Yat-sen University{p_end}
{pstd}Guangzhou, China{p_end}
{pstd}arlionn@163.com{p_end}

{pstd}Chang Liu{p_end}
{pstd}Lingnan College, Sun Yat-sen University{p_end}
{pstd}Guangzhou, China{p_end}
{pstd}liuch288@mail2.sysu.edu.cn{p_end}

{pstd}Christopher F. Parmeter{p_end}
{pstd}Department of Economics, University of Miami{p_end}
{pstd}Miami, FL, USA{p_end}
{pstd}cparmeter@bus.miami.edu{p_end}


{title:Also see}

{psee}
{space 2}Help: {help sftt_sigs}, {help sftt_eff}, {help owenst}
{p_end}
