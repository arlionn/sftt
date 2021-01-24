{smcl}
{* *! version 1.1.2  25nov2011}{...}
{cmd:help owenst}{right:also see:  {help sftt}, {help binormal}}
{hline}

{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{hi:owenst} {hline 2}}Calculate Owen's T function{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}
{cmd:owenst h a t}

{title:Description}

{pstd}
{opt owenst} calculates Owen's T function T({cmd: h}, {cmd: a}), and save the result in {cmd: t}.

{title:Examples}

    {hline}
    Setup
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 5}{p_end}
{phang2}{cmd:. generate h = _n}{p_end}
{phang2}{cmd:. generate a = 2}{p_end}
{phang2}{cmd:. generate t = .}{p_end}

{pstd}Calculate Owen's T function values{p_end}
{phang2}{cmd:. owenst h a t}{p_end}
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
{space 2}Help: {help binormal}
{p_end}
