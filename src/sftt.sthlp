{smcl}
{* *! version 0.0.1  13Jun2022}{...}
{vieweralsosee "[R] frontier" "help frontier"}{...}
{vieweralsosee "sfkk" "help sfkk"}{...}
{viewerjumpto "Syntax" "sftt##syntax"}{...}
{viewerjumpto "Description" "sftt##description"}{...}
{viewerjumpto "Estimate Options" "sftt##estimate_options_full"}{...}
{viewerjumpto "Efficiency Options" "sftt##estimate_efficiency_full"}{...}
{viewerjumpto "Stored Results" "sftt##stored_results"}{...}
{viewerjumpto "Program Authors" "sftt##program_authors"}{...}
{viewerjumpto "Reference" "sftt##reference"}{...}
{viewerjumpto "Acknowledgments" "sftt##acknowledgments"}{...}
{title:Title}

{p2colset 9 17 23 2}{...}
{p2col :{hi:sftt} {hline 2}}Two-tier stochastic frontier model{p_end}
{p2colreset}{...}


{title:Contents}

{pstd}{help sftt##syntax:Syntax}{p_end}
{pstd}{help sftt##description:Description}{p_end}
{pstd}{help sftt##estimate_options_full:Estimate Options}{p_end}
{pstd}{help sftt##estimate_efficiency_full:Efficiency Options}{p_end}
{pstd}{help sftt##examples:Examples}{p_end}
{pstd}{help sftt##stored_results:Stored Results}{p_end}
{pstd}{help sftt##program_authors:Program Authors}{p_end}
{pstd}{help sftt##reference:Reference}{p_end}
{pstd}{help sftt##acknowledgments:Acknowledgments}{p_end}


{marker syntax}{...}
{title:Syntax}

    Estimation Syntax

{p 8 17 2}
{cmd:sftt}
{depvar}
[{indepvars}]
{ifin}
[{cmd:,} {it:{help sftt##estimate_options:estimate_options}}]

    Variance decomposition subcommand

{p 8 17 2}
{cmd:sftt sigs}

    Efficiency analysis subcommand

{p 8 17 2}
{cmd:sftt eff}
[{cmd:,} {it: {help sftt##efficiency_options:efficiency_options}}]


{synoptset 28 tabbed}{...}
{marker estimate_options}{...}
{synopthdr:estimate_options}
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
{synopt :{opt iter:ate(#)}}specifies the maximum number of iterations, default is {opt iterate(1000)}. In most cases, the optimization should be converged in less than 1000 iterations{p_end}
{synopt :{opt seed(#)}}sets a random seed before estimating to ensure that the results are reproducible{p_end}
{synopt :{opt findseed}}loops through 100 estimations, during which the random seed was set from 1 to 100. In each attampt, we iterate at most 100 times{p_end}
{synopt :{cmdab:init:ial(}{it:{help nl##initial_values:initial_values}}{cmd:)}}specifies the initial values to begin the NLS estimation. This option is only used when estimating with {opt scaling}{p_end}

{syntab :SE}
{synopt :{opth vce(vcetype)}}specifies the type of standard error reported; {it:vcetype} may be {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt boot:strap}, or {opt jack:knife}{p_end}
{synopt :{opt r:obust}}synonym for {cmd:vce(robust)}{p_end}
{synoptline}


{synoptset 28 tabbed}{...}
{marker efficiency_options}{...}
{synopthdr:efficiency_options}
{synoptline}
{syntab: Select variable type}
{synopt: {opt lev:el}} only generate inefficiency terms in level specification{p_end}
{synopt: {opt exp}} only generate inefficiency terms in logarithmic specification, {opt level} and {opt exp} should not be used simultaneously{p_end}
{synopt: {opt abs:olute}} only generate absolute measures of inefficiency{p_end}
{synopt: {opt rel:ative}} only generate relative measures of inefficiency, {opt absolute} and {opt relative} should not be used simultaneously{p_end}

{syntab: Assign variable name}
{synopt: {opth u_hat(newvar)}} variable name of the estimated lower one-sided inefficiency term in level specification, default is {it: _u_hat}{p_end}
{synopt: {opth w_hat(newvar)}} variable name of the estimated upper one-sided inefficiency term in level specification, default is {it: _w_hat}{p_end}
{synopt: {opth wu_diff(newvar)}} variable name of the estimated net surplus in level specification, default is {it: _wu_diff}{p_end}
{synopt: {opth u_hat_exp(newvar)}} variable name of the estimated lower one-sided inefficiency term in logarithmic specification, default is {it: _u_hat_exp}{p_end}
{synopt: {opth w_hat_exp(newvar)}} variable name of the estimated upper one-sided inefficiency term in logarithmic specification, default is {it: _w_hat_exp}{p_end}
{synopt: {opth wu_diff_exp(newvar)}} variable name of the estimated net surplus in logarithmic specification, default is {it: _wu_diff_exp}{p_end}
{synopt: {opth wu_net_effect(newvar)}} variable name of the estimated net effect, default is {it: _wu_net_effect}{p_end}

{syntab: Other option}
{synopt: {opt replace}} permits {opt sftt} to overwrite existing variables{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{opt sftt} fits two-tier stochastic frontier (2TSF) models with multiple model settings. {opt sftt} 
provides estimators for the parameters of a linear model with a disturbance that is assumed to be a 
mixture of three components: two measures of inefficiency which are strictly nonnegative and 
nonpositive respectively, and a two-sided error term from a symmetric distribution. {opt sftt} can
fit 2TSF models with distributional assumption. When using distributional assumption mode, 
this command is applicable to estimate models in exponential/exponential/normal specification 
following {help sftt##KP09:{bind:Kumbhakar and Parmeter (2009)}} and models in half-normal/half-normal/normal 
specification following {help sftt##PP15:{bind:Papadopoulos (2015)}}. {opt sftt} also fits models
with scaling assumption following {help sftt##PA18:{bind:Parmeter (2018)}}.

{pstd}
{opt sftt sigs} identifies the distribution of each component in the composite error term.

{pstd}
{opt sftt eff} decomposes the residual and generate measures of inefficiency.


{marker estimate_options_full}{...}
{title:Estimate Options}

{dlgtab:Model setting}

{phang}
{opt noconstant} suppresses the constant term (intercept) in the frontier.

{phang}
{opt hnormal} uses the {it:half-normal/half-normal/normal} specification rather than the benchmark {it:exponential/exponential/normal} specification. Note that when using this option, the estimation might not converge because of flat derivatives or missing values, set another random seed using the option {opt seed(#)} might help. Users can also specify the {opt findseed} option to find a usable random seed.

{phang}
{opt scaling} estimates the 2TSF model with scaling property. Note that the results might be very sensitive to the initial values if the models are complex.

{dlgtab:Ancillary equations}

{phang}
{cmd:sigmau(}{help varlist:varlist_u} [,{opt noconstant}]{cmd:)}
specifies that the lower (nonpositive) inefficiency component is heteroskedastic,
with the variance expressed as a function of the covariates defined in 
{it:varlist_u}.

{phang}
{cmd:sigmaw(}{help varlist:varlist_w} [,{opt noconstant}]{cmd:)}
specifies that the lower (nonnegative) inefficiency component is heteroskedastic,
with the variance expressed as a function of the covariates defined in 
{it:varlist_w}.

{dlgtab:Other options}

{phang}
{opt seed(#)} sets a random seed before estimating to ensure that the results are reproducible.

{phang}
{opt findseed} loops through 100 estimations, during which the random seed was set from 1 to 100. In each attampt, we iterate at most 200 times.

{phang}
{opt iterate(#)} specifies the maximum number of iterations, default is {opt iterate(1000)}. In most cases, the optimization should be converged in less than 1000 iterations.

{phang}
{cmd:initial(}{it:{help nl##initial_values:initial_values}}{cmd:)} specifies the initial values to begin the NLS estimation. This option is only used when estimating with {opt scaling}.

{dlgtab:SE}

{phang}
{opt vce(vcetype)} specifies the type of standard error reported, which includes types thatare robust to some kinds of misspecification ({opt robust}), 
that allow for intragroup correlation ({opt cluster} {it:clustvar}), that use bootstrap or jackknife methods ({opt bootstrap}, {opt jackknife}) for estimation with {opt scaling}; 
see {helpb vce_option:[R] {it:vce_option}}.

{phang}
{opt robust} is the synonym for {cmd:vce(robust)}.


{marker estimate_efficiency_full}{...}
{title:Efficiency Options}

{dlgtab:Select variable type}

{phang}
{opt lev:el} only generates inefficiency terms in level specification.

{phang}
{opt exp} only generates inefficiency terms in logarithmic specification, {opt level} and {opt exp} should not be used simultaneously.

{phang}
{opt abs:olute} only generates absolute measures of inefficiency.

{phang}
{opt rel:ative} only generates relative measures of inefficiency, {opt absolute} and {opt relative} should not be used simultaneously.

{dlgtab:Assign variable name}

{phang}
{opth u_hat(newvar)} sets the variable name of the estimated lower one-sided inefficiency term in level specification, default is {it: _u_hat}.{p_end}

{phang}
{opth w_hat(newvar)} sets the variable name of the estimated upper one-sided inefficiency term in level specification, default is {it: _w_hat}.{p_end}

{phang}
{opth wu_diff(newvar)} sets the variable name of the estimated net surplus in level specification, default is {it: _wu_diff}.{p_end}

{phang}
{opth u_hat_exp(newvar)} sets the variable name of the estimated lower one-sided inefficiency term in logarithmic specification, default is {it: _u_hat_exp}.{p_end}

{phang}
{opth w_hat_exp(newvar)} sets the variable name of the estimated upper one-sided inefficiency term in logarithmic specification, default is {it: _w_hat_exp}.{p_end}

{phang}
{opth wu_diff_exp(newvar)}} sets the variable name of the estimated net surplus in logarithmic specification, default is {it: _wu_diff_exp}.{p_end}

{phang}
{opth wu_net_effect(newvar)}} sets the variable name of the estimated net effect, default is {it: _wu_net_effect}.{p_end}

{dlgtab:Other}

{phang}
{opt replace} permits {opt sftt} to overwrite existing variables.{p_end}


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}{bf:Two-tier stochastic frontier model with distributional assumptions}{p_end}
    Setup
{phang2}{cmd:. {stata "use https://sftt.oss-cn-hangzhou.aliyuncs.com/sftt_dist.dta, clear"}}{p_end}
{phang2}{cmd:. {stata "sftt, version"}}{p_end}

{pstd}Estimate model in exponential specification{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2"}}{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2, nocons"}}{p_end}

{pstd}Estimate model in half-normal specification{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2, hnormal"}}{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2, nocons hnormal"}}{p_end}

{pstd}Set random {cmd: seed} to make sure the estimation results remain constant{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2, nocons hnormal seed(20220612)"}}{p_end}

{pstd}Use {cmd: findseed} to find a usable seed{p_end}
{phang2}{cmd:. {stata "sftt y x1 x2, nocons hnormal findseed"}}{p_end}

{pstd}Identify the distributions of inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. {stata "sftt sigs"}}{p_end}

{pstd}Generate the measures of inefficiency{p_end}
{pmore}To generate all measures of inefficiency with default variable names:{p_end}
{phang2}{cmd:. {stata "sftt eff"}}{p_end}
{pmore}To generate measures of inefficiency in level specification with default variable names:{p_end}
{phang2}{cmd:. {stata "sftt eff, level replace"}}{p_end}
{pmore}To generate relative measures of inefficiency in exponential specification with customized variable names:{p_end}
{phang2}{cmd:. {stata "sftt eff, exp relative wu_diff_exp(nsurplus) wu_net_effect(neffect)"}}{p_end}

    {hline}
{pstd}{bf:Two-tier stochastic frontier model with scaling assumptions}{p_end}
    Setup
{phang2}{cmd:. {stata "use https://sftt.oss-cn-hangzhou.aliyuncs.com/sftt_scal.dta, clear"}}{p_end}

{pstd}Estimate model with scaling assumptions{p_end}
{phang2}{cmd:. {stata "sftt y x, scal sigmau(zu) sigmaw(zw) robust"}}{p_end}

{pstd}Identify the distributions of inefficiency terms{p_end}
{phang2}{cmd:. {stata "sftt sigs"}}{p_end}

    {hline}


{marker stored_results}{...}
{title:Stored results}

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


{marker program_authors}{...}
{title:Program Authors}

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


{marker reference}{...}
{title:Reference}

{marker KP09}{...}
{pstd}Kumbhakar, S.C., Parmeter, C.F. 2009.{p_end}
{pstd}{browse "https://doi.org/10.1007/s11123-008-0117-3":The effects of match uncertainty and bargaining on labor market outcomes: evidence from firm and worker specific estimates}.{p_end}
{pstd}Journal of Productivity Analysis 31(1): 1-14.{p_end}

{marker PP15}{...}
{pstd}Papadopoulos, A. 2015.{p_end}
{pstd}{browse "https://doi.org/10.1007/s11123-014-0389-8":The half-normal specification for the two-tier stochastic frontier model}.{p_end}
{pstd}Journal of Productivity Analysis 43(2): 225-230.{p_end}

{marker PA18}{...}
{pstd}Parmeter, C. F. 2018.{p_end}
{pstd}{browse "https://doi.org/10.1007/s11123-017-0520-8":Estimation of the two-tiered stochastic frontier model with the scaling property}.{p_end}
{pstd}Journal of Productivity Analysis 49(1): 37-47.{p_end}


{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}We thank Alecos Papadopoulos for his amazing support.


{title:Also see}

{psee}
{space 2}Help: {help sftt_sigs}, {help sftt_eff}, {help owenst}
{p_end}
