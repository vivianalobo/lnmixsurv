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
      (Intercept)_1 4.0455659 0.006567584 4.0368490 4.0547803
      x1_1          0.8072690 0.009969275 0.7946591 0.8201348
      (Intercept)_2 3.4262658 0.020706785 3.3997479 3.4528735
      x1_2          0.4874025 0.023423408 0.4582146 0.5154283
      
      Auxiliary parameter(s):
              estimate  std.error   conf.low conf.high
      phi_1 26.6440591 1.41336921 24.9120464 28.452973
      phi_2  3.1884579 0.12109791  3.0417389  3.334675
      eta_1  0.5072736 0.01229135  0.4913598  0.523180

