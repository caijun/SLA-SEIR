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
    qs[1, , t] <- quantile(S, c(.05, .5, .95))
    qs[2, , t] <- quantile(I, c(.05, .5, .95))
    qs[3, , t] <- quantile(E, c(.05, .5, .95))
    qs[4, , t] <- quantile(R, c(.05, .5, .95))
    qs[5, , t] <- quantile(pars[, 2], c(.05, .5, .95))
    qs[6, , t] <- quantile(pars[, 1], c(.05, .5, .95))
    qs[7, , t] <- quantile(pars[, 3], c(.05, .5, .95))
    qs[8, , t] <- quantile(sqrt(sig2), c(.05, .5, .95))
    qs[9, , t] <- quantile(sqrt(tau2), c(.05, .5, .95))
    qs[10, , t] <- quantile(x, c(.05, .5, .95))
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

# one step ahead forecast
PF1 <- function(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a) {
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
  loglike <- rep(0, n)
  pred <- matrix(0, M, n) # predictions of the infected
  for (t in 1:n) {
    print(t)
    pred[, t] <- (1 - pars[, 3]) * I + pars[, 1] * E
    w <- dnorm(yI[t], -pars[, 3] + pars[, 1] * E / I, sqrt(sig2 + tau2), log = TRUE)
    loglike[t] <- log(mean(exp(w) / y[t]))
    ### same code from PF function ###
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
    qs[1, , t] <- quantile(S, c(.05, .5, .95))
    qs[2, , t] <- quantile(I, c(.05, .5, .95))
    qs[3, , t] <- quantile(E, c(.05, .5, .95))
    qs[4, , t] <- quantile(R, c(.05, .5, .95))
    qs[5, , t] <- quantile(pars[, 2], c(.05, .5, .95))
    qs[6, , t] <- quantile(pars[, 1], c(.05, .5, .95))
    qs[7, , t] <- quantile(pars[, 3], c(.05, .5, .95))
    qs[8, , t] <- quantile(sqrt(sig2), c(.05, .5, .95))
    qs[9, , t] <- quantile(sqrt(tau2), c(.05, .5, .95))
    qs[10, , t] <- quantile(x, c(.05, .5, .95))
  }
  return(list(qs = qs, loglike = loglike, pred = pred))
}

ar1plusnoise <- function(y, yI, q0, Q0, a0, b0, c0, d0, alpha, beta, V, W, g0) {
  y <- y / 100
  n <- length(yI)
  M <- length(g0)
  qs <- array(0, c(n, 6, 3))
  s <- matrix(0, M, 9) # sufficient statisitics
  s[, 1] <- 1.0 / Q0[1, 1] # Q_t^{-1}
  s[, 2] <- 0.0
  s[, 3] <- 1.0 / Q0[2, 2]
  s[, 4] <- q0[1] / Q0[1, 1] # Q_t^{-1}q_t = q_t / Q_t
  s[, 5] <- q0[2] / Q0[2, 2]
  s[, 6] <- c0 # c_t
  s[, 7] <- d0 # d_t
  s[, 8] <- a0 # a_t
  s[, 9] <- b0 # b_t
  qt1 <- rep(q0[1], M) # q_t, namely q_t = (mu, phi)
  qt2 <- rep(q0[2], M)
  I <- rep(y[1], M)
  loglike <- rep(0, n)
  pred <- matrix(0, n, M)
  for (t in 1:n) {
    print(t)
    g1 <- alpha + beta * g0 # g_t = mu + phi * g_{t-1}
    pred[t, ] <- I * (1 + g1)
    w <- dnorm(yI[t], g1, sqrt(V + W))
    loglike[t] <- log(mean(w / y[t]))
    k <- sample(1:M, size = M, replace = TRUE, prob = w)
    w <- w[k]
    g0o <- g0[k]
    g1 <- g1[k]
    alpha <- alpha[k]
    beta <- beta[k]
    W <- W[k]
    V <- V[k]
    so  <- s[k, ]
    qt1o <- qt1[k]
    qt2o <- qt2[k]
    var  <- 1 / (1 / V + 1 / W)
    mean <- var * (yI[t] / V + g1 / W)
    sd <- sqrt(var)
    g0 <- rnorm(M, mean, sd) # (g_t|s_t^k, epison) ~ N(m_t, C_t)
    s[, 1] <- so[, 1] + 1
    s[, 2] <- so[, 2] + g0o
    s[, 3] <- so[, 3] + g0o ^ 2
    s[, 4] <- so[, 4] + g0
    s[, 5] <- so[, 5] + g0 * g0o
    s[, 6] <- so[, 6] + 1 / 2 # c_t = c_{t-1} + 1/2
    m <- s[, 1] * s[, 3] - s[, 2] ^ 2 # determinant of Q_t^{-1}
    qt1 <- (s[, 3] * s[, 4] - s[, 2] * s[, 5]) / m
    qt2 <- (s[, 1] * s[, 5] - s[, 2] * s[, 4]) / m
    # y_t is replaced by g_t to avoid NA
    s[, 7] <- so[, 7] + (g0 - qt1 - qt2 * g0o) * g0 / 2 + 
      ((qt1o - qt1) * so[, 4] + (qt2o - qt2) * so[, 5]) / 2 # d_t = d_{t-1} + ...
    s[, 8] <- so[, 8] + 1 / 2 # a_t = a_{t-1} + 1/2
    s[, 9] <- so[, 9] + (yI[t] - g0) ^ 2 / 2 # b_t = b_{t-1} + (y_t - g_t)^2 / 2
    W <- 1 / rgamma(M, s[, 6], s[, 7]) # (W|g^t, X^t) ~ IG(c_t, d_t)
    std <- sqrt(W / m)
    norm <- cbind(rnorm(M, 0, std), rnorm(M, 0, std))
    alpha <- qt1 + sqrt(s[, 3]) * norm[, 1]
    beta <- qt2 - s[, 2] / sqrt(s[, 3]) * norm[, 1] + 
      sqrt(s[, 1] - s[, 2] ^ 2 / s[, 3]) * norm[, 2]
    V <- 1 / rgamma(M, s[, 8], s[, 9]) # (V|y^t, g^t) ~ IG(a_t, b_t)
    I <- I * (1 + g0)
    qs[t, 1, ] <- quantile(alpha, c(.05, .5, .95))
    qs[t, 2, ] <- quantile(beta, c(.05, .5, .95))
    qs[t, 3, ] <- quantile(W, c(.05, .5, .95))
    qs[t, 4, ] <- quantile(V, c(.05, .5, .95))
    qs[t, 5, ] <- quantile(g0, c(.05, .5, .95))
    qs[t, 6, ] <- quantile(I, c(.05, .5, .95))
  }
  return(list(qs = qs, loglike = loglike, pred = pred))
}
