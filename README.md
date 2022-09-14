# SFTT

The repository of Stata command `sftt`. For details, see:

> Yujun Lian, Chang Liu, Christopher F. Parmeter. Two-tier stochastic frontier analysis using Stata. **Stata Joural**, 2022, forthcoming.


## Description
`sftt` fits two-tier stochastic frontier (**2TSF**) models with multiple model settings. 

The 2TSF model consists of a linear model with a disturbance that is assumed to be a mixture of three components: 
two measures of inefficiency which are strictly nonnegative and nonpositive respectively,
and a two-sided error term from a symmetric distribution.

`sftt` can fit 2TSF models with distributional assumption.
When using distributional assumption mode, 
this command is applicable to estimate 
- models in **exponential/exponential/normal** specification
following [Kumbhakar and Parmeter (2009)](https://doi.org/10.1007/s11123-008-0117-3) 
- models in **half-normal/half-normal/normal** specification following
[Papadopoulos (2015)](https://doi.org/10.1007/s11123-014-0389-8).

This command also fits 
- models with **scaling property** following
[Parmeter (2018)](https://doi.org/10.1007/s11123-017-0520-8).

## The sftt commands 
- `sftt` estimates two-tier SF models listed above.
- `sftt sigs` identifies the distribution of each component in the composite error term.
- `sftt eff` decomposes the residual and generate measures of inefficiency.


## Acknowledgments
We thank Alecos Papadopoulos for his amazing support.


## Program Authors
- [Yujun Lian](mailto:arlionn@163.com) (repo owner).
Lingnan College, Sun Yat-sen University.
Guangzhou, China. [-Link-](https://www.lianxh.cn) 
- [Chang Liu](mailto:liuch288@mail2.sysu.edu.cn).
Lingnan College, Sun Yat-sen University
Guangzhou, China.
- [Christopher F. Parmeter](mailto:cparmeter@bus.miami.edu).
Department of Economics, University of Miami
Miami, FL, USA.
