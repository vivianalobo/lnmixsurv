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
      (Intercept)_1 4.0445030 0.006522968 4.0353768 4.0527129
      x1_1          0.8087692 0.010299547 0.7970843 0.8217099
      (Intercept)_2 3.4282184 0.019235326 3.4060934 3.4539586
      x1_2          0.4877074 0.021489632 0.4575755 0.5158396
      
      Auxiliary parameter(s):
              estimate  std.error   cred.low  cred.high
      phi_1 26.7509482 1.38834169 25.1388200 28.7978371
      phi_2  3.1690732 0.11384819  3.0266862  3.3031788
      eta_1  0.5047119 0.01205289  0.4889579  0.5208052

