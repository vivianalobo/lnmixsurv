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
      (Intercept)_1 4.0446964 0.006796203 4.0359851 4.0533994
      x1_1          0.8095789 0.009936305 0.7969076 0.8222505
      (Intercept)_2 3.4277493 0.018647274 3.4033308 3.4498217
      x1_2          0.4868316 0.021226930 0.4593627 0.5118913
      
      Auxiliary parameter(s):
              estimate  std.error   conf.low  conf.high
      phi_1 26.5640273 1.21094188 25.1844674 28.1430092
      phi_2  3.1826322 0.11464682  3.0543650  3.3391176
      eta_1  0.5051077 0.01192346  0.4913414  0.5220004

