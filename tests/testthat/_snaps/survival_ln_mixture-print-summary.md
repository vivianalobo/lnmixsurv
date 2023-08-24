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
      (Intercept)_1 4.0456457 0.006370159 4.0369899 4.0529927
      x1_1          0.8075922 0.009649999 0.7948419 0.8188464
      (Intercept)_2 3.4202430 0.020753769 3.3928802 3.4463243
      x1_2          0.4879574 0.022196693 0.4613753 0.5157443
      
      Auxiliary parameter(s):
              estimate  std.error   conf.low  conf.high
      phi_1 26.3151618 1.33872041 24.6999617 28.0830992
      phi_2  3.2160359 0.12065186  3.0759101  3.3784258
      eta_1  0.5115556 0.01311835  0.4948066  0.5286459

