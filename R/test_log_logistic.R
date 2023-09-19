X = flexsurv::rllogis(n = 10000, shape = 2, scale = 3)
C = flexsurv::rllogis(n = 10000, shape = 3, scale = 4)
W = X #pmin(X, C)
D = rep(1, 10000) #as.numeric(X <= C)
table(D)

fit = survival::survreg(formula = survival::Surv(time = W, event = D) ~ 1, 
                        dist = "loglogistic")
exp(-fit$coefficients/fit$scale)
1/fit$scale ## shape

temp = survival::psurvreg(q = X, 
                   mean = fit$linear.predictors, 
                   scale = fit$scale, 
                   distribution = "loglogistic")

alpha =  1 / fit$scale
lambda = exp(unique(fit$linear.predictors))
data.frame(X) |> 
  ggplot(aes(x = X)) +
  stat_function(fun = function(x) survival::psurvreg(q = x, 
                                                     mean = unique(fit$linear.predictors), 
                                                     scale = fit$scale, 
                                                     distribution = "loglogistic")) +
  stat_function(fun = function(x) flexsurv::pllogis(q = x, shape = 2, scale = 3), linetype = 2, color = "blue") +
  stat_function(fun = function(x) flexsurv::pllogis(q = x, shape =  1 / fit$scale, scale = exp(unique(fit$linear.predictors))), linetype = 2, color = "red") + 
  stat_function(fun = function(x) 1 - 1 / (1 + (x / lambda) ^ alpha), linetype = 2, color = "yellow")

summary((1 - 1 / (1 + (1:600 / lambda) ^ alpha)) - flexsurv::pllogis(q = 1:600, shape =  1 / fit$scale, scale = exp(unique(fit$linear.predictors))))
