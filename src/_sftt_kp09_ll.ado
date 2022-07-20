/* KP09 exp-exp-normal log likelihood */
program define _sftt_kp09_ll
    version 13
    args lnf xb sigma_v sigma_u sigma_w

    tempvar e a b eta beta tmpsum
    quietly generate double    `e' = $ML_y1 - `xb'
    quietly generate double    `a' = `e' / exp(`sigma_u') + exp(`sigma_v')^2 / (2 * exp(`sigma_u')^2)
    quietly generate double    `b' = `e' / exp(`sigma_v') - exp(`sigma_v') / exp(`sigma_w')
    quietly generate double  `eta' = exp(`sigma_v')^2 / (2 * exp(`sigma_w')^2) - `e' / exp(`sigma_w')
    quietly generate double `beta' = - `e' / exp(`sigma_v') - exp(`sigma_v') / exp(`sigma_u')
    quietly replace `lnf' = -ln(exp(`sigma_u') + exp(`sigma_w')) + ///
                             ln(exp(`a') * normal(`beta') + exp(`eta') * normal(`b'))
end
