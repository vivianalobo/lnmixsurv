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
      (Intercept)_a 4.0445757 0.006236653 4.0368830 4.0526401
      x1_a          0.8090734 0.009518335 0.7958502 0.8199592
      (Intercept)_b 3.4245488 0.019522311 3.4006760 3.4501877
      x1_b          0.4915648 0.021621560 0.4627559 0.5175077
      
      Auxiliary parameter(s):
                estimate  std.error   conf.low  conf.high
      phi_a   26.3862345 1.34253028 24.7298691 28.1202304
      phi_b    3.1816509 0.10961297  3.0511497  3.3148766
      theta_a  0.5088579 0.01305763  0.4908965  0.5243082

