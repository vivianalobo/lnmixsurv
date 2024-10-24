# print method works

    Code
      mod
    Output
      survival_ln_mixture
       formula: survival::Surv(y, delta) ~ x
       observations: 10000
       predictors: 2
       mixture groups: 2
      ------------------
                     estimate   std.error  cred.low cred.high
      (Intercept)_1 4.0446935 0.006659593 3.9118220 4.0574478
      x1_1          0.8062953 0.010040927 0.7861331 0.8271481
      (Intercept)_2 3.4235683 0.020073543 3.2987428 3.4611767
      x1_2          0.4865902 0.021463598 0.4158180 0.5315404
      
      Auxiliary parameter(s):
              estimate  std.error cred.low  cred.high
      phi_1 26.2490951 1.47347833 6.037256 29.1153792
      phi_2  3.1936586 0.11375225 2.901734  3.7560002
      eta_1  0.5089456 0.01261907 0.489533  0.7300081

