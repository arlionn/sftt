# SFTT

The repository of Stata command `sftt`. For details, see:

> Lian, Y., Liu, C., & Parmeter, C. F. (2023). Two-tier stochastic frontier analysis using Stata. The Stata Journal, 23(1), 197–229. [link](https://journals.sagepub.com/doi/abs/10.1177/1536867X231162033), [-PDF-](
https://file.lianxh.cn/Refs/LianPub/Lian-2023-SJ-sftt-Two-tier-SFA.pdf)


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

## Install
You can always type `search sftt` in Stata's command window to get access to package. 
Or, you can use the following commands to download it directly.
```stata
net install st0705.pkg, replace
net get     st0705.pkg, replace  // to get main.do file
```
Then you can read the help document to get more detailed information:
```stata
help sftt
```

## Example
 First add directory to end of ado-path.
 ```
 adopath + "./src"
 ```

Load the data used in [Kumbhakar and Parmeter (2009)](https://doi.org/10.1007/s11123-008-0117-3) and replicate their results.

```
use "https://sftt.oss-cn-hangzhou.aliyuncs.com/kp09.dta", clear
sftt lwage iq educ educ2 exper exper2 tenure tenure2  ///
     age married south urban black sibs brthord meduc feduc
```

Finally, you can use the post-estimation commands `sftt sigs` and `sftt eff` to assist your efficiency analysis.

```
sftt sigs
sftt eff, replace
```

You can use `help sftt` to see more detailed instructions.

## Files
- `./mc_results` stores the results generated by the Monte Carlo simulation using `scaling_mc.do`.
- `./output` stores the results shown in Lian et al. (2022).
- `./src` stores the source code of `sftt`.
- `main.do` provides example of the applications of `sftt`, the results will be saved in `./output`.
- `scaling_mc.do` implements the Monte-Carlo simulation following [Parmeter (2018)](https://doi.org/10.1007/s11123-017-0520-8), the results will be saved in `./mc_results`.


## References

- Kumbhakar, S. C., and C. F. Parmeter.  2009.  The effects of match uncertainty and bargaining on labor market outcomes: Evidence from firm and worker specific estimates.  Journal of Productivity Analysis 31: 1–14.  https://doi.org/10.1007/s11123-008-0117-3.
- Papadopoulos, A. A.  2015.  The half-normal specification for the two-tier stochastic frontier model.  Journal of Productivity Analysis 43: 225–230.  https://doi.org/10.1007/s11123-014-0389-8.
- Parmeter, C. F.  2018.  Estimation of the two-tiered stochastic frontier model with the scaling property.  Journal of Productivity Analysis 49: 37–47.  https://doi.org/10.1007/s11123-017-0520-8.


## Acknowledgments
We thank Dr. Jenkins and the anonymous reviewer for their valuable and insightful comments.

We also thank Alecos Papadopoulos for his helpful support.


## Authors
- [Yujun Lian](mailto:arlionn@163.com) (repo owner).
Lingnan College, Sun Yat-sen University.
Guangzhou, China. 
[profile](https://lingnan.sysu.edu.cn/en/faculty/LianYujun),
[blog](https://www.lianxh.cn) 
- [Chang Liu](mailto:liuch288@mail2.sysu.edu.cn).
Lingnan College, Sun Yat-sen University
Guangzhou, China.
- [Christopher F. Parmeter](mailto:cparmeter@bus.miami.edu).
Department of Economics, University of Miami
Miami, FL, USA.
[profile](https://people.miami.edu/profile/c.parmeter@miami.edu)
