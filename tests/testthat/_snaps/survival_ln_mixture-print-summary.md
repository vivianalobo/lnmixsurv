# print method works

    Code
      mod
    Output
      survival_ln_mixture
       formula:      survival::Surv(y, delta) ~ x
       observations: 10000
       predictors:   2
      ------
                     estimate   std.error  conf.low conf.high
      (Intercept)_a 4.0448076 0.006356691 4.0364529 4.0528428
      x1_a          0.8095204 0.010001643 0.7970927 0.8221297
      (Intercept)_b 3.4264007 0.022319628 3.4015256 3.4553048
      x1_b          0.4888544 0.019117639 0.4615721 0.5144357
      
      Auxiliary parameter(s):
                estimate  std.error   conf.low  conf.high
      phi_a   26.5081615 1.37813559 24.9220177 28.5826875
      phi_b    3.1773781 0.12042904  3.0219889  3.3345268
      theta_a  0.5055948 0.01349581  0.4877678  0.5231282

