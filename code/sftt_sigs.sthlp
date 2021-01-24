{smcl}
{* *! version 1.1.2  25nov2011}{...}
{cmd:help sftt_sigs}{right:also see:  {help sftt}, {help sftt_eff}}
{hline}

{title:Title}

{p2colset 5 22 23 2}{...}
{p2col :{hi:sftt_sigs} {hline 2}}Variance identification in two-tiered stochastic frontier model{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}
{cmd:sftt_sigs}


{title:Description}

{pstd}
{opt sftt_sigs} identifies the variances of each component in the composite error term.


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

{pstd}Two-tiered stochastic frontier model with exponential distribution for
inefficiency terms{p_end}
{phang2}{cmd:. sftt y x1 x2, check search nocons}{p_end}

{pstd}Identify the variances of inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. sftt_sigs}{p_end}

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

{pstd}Two-tiered stochastic frontier model with scaling assumption for
inefficiency terms{p_end}
{phang2}{cmd:. sftt y x, scal sigmau(zu) sigmaw(zw) robust}{p_end}

{pstd}Identify the variances of inefficiency terms and stochastic noise{p_end}
{phang2}{cmd:. sftt_sigs}{p_end}
    {hline}

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
{space 2}Help: {help sftt}, {help sftt_eff}, {help owenst}
{p_end}
