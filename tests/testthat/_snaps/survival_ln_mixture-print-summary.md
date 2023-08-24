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
      (Intercept)_1 4.0449468 0.006640402 4.0362772 4.0532077
      x1_1          0.8091331 0.009723673 0.7971570 0.8220558
      (Intercept)_2 3.4289818 0.017118112 3.4075577 3.4506387
      x1_2          0.4880462 0.021268204 0.4611461 0.5170355
      
      Auxiliary parameter(s):
             estimate   std.error  conf.low  conf.high
      phi_1 26.852790 1.250917196 25.255053 28.5227973
      phi_2  3.171878 0.090543019  3.055588  3.3002309
      eta_1  0.504822 0.009105166  0.492050  0.5176322

