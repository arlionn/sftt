/* Pa15 hn-hn-normal log likelihood */
program define _sftt_pa15_ll
    version 13
    args lnf xb sigma_v sigma_u sigma_w

    tempvar e theta1 theta2 s omega1 omega2 lambda1 lambda2 tmp_sqrt
    quietly {
        generate double   `e' = $ML_y1 - `xb'
        generate double `theta1'   = exp(`sigma_w') / exp(`sigma_v')
        generate double `theta2'   = exp(`sigma_u') / exp(`sigma_v')
        generate double `tmp_sqrt' = sqrt(1 + `theta1'^2 + `theta2'^2)
        generate double `s'        = exp(`sigma_v') * `tmp_sqrt'
        generate double `omega1'   = (`s' * sqrt(1 + `theta2'^2)) / `theta1'
        generate double `omega2'   = (`s' * sqrt(1 + `theta1'^2)) / `theta2'
        generate double `lambda1'  = (`theta2' / `theta1') * `tmp_sqrt'
        generate double `lambda2'  = (`theta1' / `theta2') * `tmp_sqrt'
    }

    // generate Gs with bivariate normal
    tempvar G1 G2 rho1 rho2
    quietly {
        generate double `rho1' = `lambda1' / sqrt(1 + `lambda1'^2)
        generate double `rho2' = - `lambda2' / sqrt(1 + `lambda2'^2)
        generate double `G1'   = 2 * binormal(`e' / `omega1', 0, `rho1')
        generate double `G2'   = 2 * binormal(`e' / `omega2', 0, `rho2')
    }
    quietly replace `lnf' = - ln(`s') - (1 / (2 * `s'^2)) * `e'^2 + ln(`G1' - `G2')
end
