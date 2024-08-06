# print method works

    Code
      mod
    Output
      survival_ln_mixture
       formula:      survival::Surv(y, delta) ~ x
       observations: 10000
       predictors:   2
      ------
                     estimate   std.error  cred.low cred.high
      (Intercept)_1 3.5340990 0.013925815 3.5158592 3.5519913
      x1_1          0.6849315 0.017822534 0.6611903 0.7086563
      (Intercept)_2 4.0389569 0.007380801 4.0297442 4.0482417
      x1_2          6.1613494 0.065996346 6.0998032 6.2806298
      
      Auxiliary parameter(s):
              estimate  std.error   cred.low  cred.high
      phi_1  2.5988739 0.05301561  2.5308990  2.6689762
      phi_2 31.7417973 1.90775508 29.4587221 34.2737305
      eta_1  0.5927116 0.00675987  0.5837615  0.6010928

