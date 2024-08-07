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
      (Intercept)_1 4.0452448 0.006677706 4.0370042 4.0532995
      x1_1          0.8083585 0.009417680 0.7960475 0.8206874
      (Intercept)_2 3.4279174 0.017201098 3.4044733 3.4514542
      x1_2          0.4865498 0.020273239 0.4593043 0.5144910
      
      Auxiliary parameter(s):
              estimate  std.error   cred.low  cred.high
      phi_1 26.6787610 1.34044730 25.0619765 28.5371591
      phi_2  3.1788775 0.11399542  3.0434570  3.3185873
      eta_1  0.5051823 0.01094656  0.4912011  0.5193892

