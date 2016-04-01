# SLA-SEIR
Sequential Learning Algorithm for state space SEIR epidemiological model

This is an optimized R code  for sequential learning algorithm for state space SEIR epidemiological model used in Dukic, V., Lopes, H. F., & Polson, N. G. (2012). [Tracking epidemics with Google flu trends data and a state-space SEIR model](http://www.tandfonline.com/doi/abs/10.1080/01621459.2012.713876). Journal of the American Statistical Association, 107(500), 1410-1426. I commented and modified the code to facilitate my understanding to Dukic et al.'s paper. For the raw R code please contact the authors.

After reviewing Dukic et al.'s code, I find that there are some problems:

* In the implementation of `ar1plusnoise` function, Dukic et al. update $b_t$ by $b_t = (y_t - g_t)^2$ rather than by the iterative formula in their paper, $b_t = (y_t - g_t)^2/2$.