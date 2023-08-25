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
      (Intercept)_1 4.0454719 0.005770805 4.0380584 4.0534575
      x1_1          0.8060473 0.009030251 0.7947023 0.8176975
      (Intercept)_2 3.4203686 0.021000882 3.3940826 3.4480849
      x1_2          0.4880618 0.021694714 0.4601262 0.5162373
      
      Auxiliary parameter(s):
              estimate  std.error   conf.low  conf.high
      phi_1 26.2318351 1.38610133 24.3733834 27.9517155
      phi_2  3.2016939 0.11635714  3.0615756  3.3523516
      eta_1  0.5105744 0.01318735  0.4933958  0.5263495

