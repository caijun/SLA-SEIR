rtnorm <- function(N, mu, sd) {
  u <- runif(N)
  mu + sd * qnorm(u + (1 - u) * pnorm(-mu / sd))
}

PF <- function(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a) {
  y <- y / 100
  n <- length(yI)
  # a is the Liu-West shrinage factor
  # h is a controlling smoothing parameter
  h2 <- 1 - a ^ 2
  sig2 <- 1 / rgamma(M, a0, b0) # evolution variance, \sigma_g^2
  tau2 <- 1 / rgamma(M, c0, d0) # observation variance, \sigma_y^2
  alpha <- rtnorm(M, ma, sda)
  beta <- rtnorm(M, mb, sdb)
  gamma <- rtnorm(M, mg, sdg)
  S <- rep(1 - y[1], M)
  E <- rep(0, M)
  I <- rep(y[1], M)
  R <- rep(0, M)
  ss <- rep(b0, M) # parameter learning for \sigma_g^2
  vv <- rep(d0, M) # parameter learning for \sigma_y^2
  qs <- array(0, c(10, 3, n)) # quantiles of 10 parameters of interest
  pars <- cbind(alpha, beta, gamma)
  lpars <- log(pars)
  for (t in 1:n) {
    print(t)
    mp <- apply(lpars, 2, mean)
    vp <- var(lpars)
    lpars <- a * lpars + (1 - a) * matrix(mp, M, 3, byrow = TRUE)
    pars <- exp(lpars)
    # Resampling
    g <- -pars[, 3] + pars[, 1] * E / I
    w <- dnorm(yI[t], g, sqrt(sig2 + tau2), log = TRUE) # probabilities p are given as log(p)
    w1 <- exp(w - max(w))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    w <- w[k]
    S1 <- S[k]
    I1 <- I[k]
    E1 <- E[k]
    R1 <- R[k]
    sig2 <- sig2[k]
    ss1 <- ss[k]
    tau2 <- tau2[k]
    vv1 <- vv[k]
    lpars <- lpars[k, ]
    # Propagation
    lpars <- lpars + matrix(rnorm(M * 3), M, 3) %*% chol(h2 * vp)
    pars <- exp(lpars)
    g1 <- -pars[, 3] + pars[, 1] * E1 / I1
    var <- 1 / (1 / sig2 + 1 / tau2)
    mean <- var * (yI[t] / sig2 + g1 / tau2)
    x <- rnorm(M, mean, sqrt(var))
    I <- I1 * (1 + x)
    # Update (E, R, S)
    E <- pars[, 2] * I1 * S1 + (1 - pars[, 1]) * E1
    R <- R1 + pars[, 3] * I1
    S <- 1 - I - R - E
    # Reweighting
    w1 <- dnorm(yI[t], g1, sqrt(sig2 + tau2), log = TRUE) - w
    w1 <- exp(w1 - max(w1))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    S <- S[k]
    E <- E[k]
    I <- I[k]
    R <- R[k]
    sig2 <- sig2[k]
    ss1 <- ss1[k]
    tau2 <- tau2[k]
    vv1 <- vv1[k]
    lpars <- lpars[k, ]
    # Offline sampling fixed parameters
    ss <- ss1 + (yI[t] - x) ^ 2 / 2
    sig2 <- 1 / rgamma(M, a0 + t / 2, ss)
    vv <- vv1 + (x - g1) ^ 2 / 2
    tau2 <- 1 / rgamma(M, c0 + t / 2, vv)
    # Storage
    pars <- exp(lpars)
    # cbind is time consuming because of copying data.
    # pmat <- cbind(S, I, E, R, pars[, c(2, 1, 3)], sqrt(sig2), sqrt(tau2), x)
    # qs[, , t] <- t(apply(pmat, 2, quantile, probs = c(.05, .5, .95)))
    qs[1, , t] <- quantile(S, c(0.05, 0.5, 0.95))
    qs[2, , t] <- quantile(I, c(0.05, 0.5, 0.95))
    qs[3, , t] <- quantile(E, c(0.05, 0.5, 0.95))
    qs[4, , t] <- quantile(R, c(0.05, 0.5, 0.95))
    qs[5, , t] <- quantile(pars[, 2], c(0.05, 0.5, 0.95))
    qs[6, , t] <- quantile(pars[, 1], c(0.05, 0.5, 0.95))
    qs[7, , t] <- quantile(pars[, 3], c(0.05, 0.5, 0.95))
    qs[8, , t] <- quantile(sqrt(sig2), c(0.05, 0.5, 0.95))
    qs[9, , t] <- quantile(sqrt(tau2), c(0.05, 0.5, 0.95))
    qs[10, , t] <- quantile(x, c(0.05, 0.5, 0.95))
  }
  return(qs)
}

PFlike <- function(y, yI, alpha, beta, gamma, a0, b0, c0, d0, M) {
  y <- y / 100
  n <- length(yI)
  sig2 <- 1 / rgamma(M, a0, b0)
  tau2 <- 1 / rgamma(M, c0, d0)
  S <- rep(1 - y[1], M)
  E <- rep(0, M)
  I <- rep(y[1], M)
  R <- rep(0, M)
  ss <- rep(b0, M)
  vv <- rep(d0, M)
  loglike <- rep(0, n)
  for (t in 1:n) {
    print(t)
    # Resample (S, E, I, R)
    g <- -gamma + alpha * E / I
    w <- dnorm(yI[t], g, sqrt(sig2 + tau2), log = TRUE)
    w1 <- exp(w - max(w))
    loglike[t] <- log(mean(exp(w)))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    S1 <- S[k]
    E1 <- E[k]
    I1 <- I[k]
    R1 <- R[k]
    sig2 <- sig2[k]
    ss1 <- ss[k]
    tau2 <- tau2[k]
    vv1 <- vv[k]
    # Update I
    g1 <- g[k]
    var <- 1 / (1 / sig2 + 1 / tau2)
    mean <- var * (yI[t] / sig2 + g1 / tau2)
    x <- rnorm(M, mean, sqrt(var))
    I <- I1 * (1 + x)
    # Update (E, R, S)
    E <- beta * I1 * S1 + (1 - alpha) * E1
    R <- R1 + gamma * I1
    S <- 1 - I - R - E
    # Offline sampling fixed parameters
    ss <- ss1 + (yI[t] - x) ^ 2 / 2
    sig2 <- 1 / rgamma(M, a0 + t / 2, ss)
    vv <- vv1 + (x - g1) ^ 2 / 2
    tau2 <- 1 / rgamma(M, c0 + t / 2, vv)
  }
  return(loglike)
}
