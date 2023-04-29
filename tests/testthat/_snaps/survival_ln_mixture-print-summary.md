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
      (Intercept)_a 4.0438332 0.006927102 4.0355844 4.0537543
      x1_a          0.8084604 0.009819686 0.7963534 0.8215608
      (Intercept)_b 3.4227421 0.021024321 3.3921034 3.4491069
      x1_b          0.4902553 0.019388447 0.4614154 0.5130385
      
      Auxiliary parameter(s):
                estimate  std.error   conf.low  conf.high
      phi_a   26.2474165 1.35605326 24.5316535 28.0389219
      phi_b    3.2016751 0.11714397  3.0565605  3.3568732
      theta_a  0.5090044 0.01367477  0.4936986  0.5278656

