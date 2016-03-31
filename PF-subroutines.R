rtnorm <- function(N, mu, sd) {
  u <- runif(N)
  mu + sd * qnorm(u + (1 - u) * pnorm(-mu / sd))
}

# ahead = 1, one step ahead forecast
PF <- function(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a, ahead = 0) {
  y <- y / 100
  n <- length(yI)
  # a is the Liu-West shrinage factor
  # h is a controlling smoothing parameter
  h2 <- 1 - a ^ 2
  sig_g2 <- 1 / rgamma(M, a0, b0) # evolution variance, \sigma_g^2
  sig_y2 <- 1 / rgamma(M, c0, d0) # observation variance, \sigma_y^2
  alpha <- rtnorm(M, ma, sda)
  beta <- rtnorm(M, mb, sdb)
  gamma <- rtnorm(M, mg, sdg)
  S_t <- rep(1 - y[1], M)
  E_t <- rep(0, M)
  I_t <- rep(y[1], M)
  R_t <- rep(0, M)
  b0_t <- rep(b0, M) # parameter learning for \sigma_g^2
  d0_t <- rep(d0, M) # parameter learning for \sigma_y^2
  pars <- cbind(alpha, beta, gamma)
  lpars <- log(pars)
  qs <- array(0, c(10, 3, n)) # quantiles of 10 parameters of interest
  loglike <- rep(0, n)
  pred <- matrix(0, M, n) # predictions of the infected
  for (t in 1:n) {
    pred[, t] <- (1 - pars[, 3]) * I_t + pars[, 1] * E_t # one step ahead forecast
    mu_g <- -pars[, 3] + pars[, 1] * E_t / I_t
    w <- dnorm(yI[t], mu_g, sqrt(sig_g2 + sig_y2), log = TRUE) # probabilities p are given as log(p)
    loglike[t] <- log(mean(exp(w) / y[t]))
    ### same code from PF function ###
    print(t)
    m_psi <- apply(lpars, 2, mean)
    V_psi <- var(lpars)
    lpars <- a * lpars + (1 - a) * matrix(m_psi, M, 3, byrow = TRUE) # a is eta
    pars <- exp(lpars)
    # Resampling
    mu_g <- -pars[, 3] + pars[, 1] * E_t / I_t
    w <- dnorm(yI[t], mu_g, sqrt(sig_g2 + sig_y2), log = TRUE) # probabilities p are given as log(p)
    w1 <- exp(w - max(w))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    w <- w[k]
    lpars <- lpars[k, ]
    sig_g2 <- sig_g2[k]
    sig_y2 <- sig_y2[k]
    S_t0 <- S_t[k]
    I_t0 <- I_t[k]
    E_t0 <- E_t[k]
    R_t0 <- R_t[k]
    b0_t0 <- b0_t[k]
    d0_t0 <- d0_t[k]
    # Propagation
    lpars <- lpars + matrix(rnorm(M * 3), M, 3) %*% chol(h2 * V_psi)
    pars <- exp(lpars)
    mu_g <- -pars[, 3] + pars[, 1] * E_t0 / I_t0
    B <- 1 / (1 / sig_g2 + 1 / sig_y2)
    b <- B * (yI[t] / sig_g2 + mu_g / sig_y2)
    g_t <- rnorm(M, b, sqrt(B))
    I_t <- I_t0 * (1 + g_t)
    # Update (E, R, S)
    E_t <- pars[, 2] * I_t0 * S_t0 + (1 - pars[, 1]) * E_t0
    R_t <- R_t0 + pars[, 3] * I_t0
    S_t <- 1 - I_t - R_t - E_t
    # Reweighting
    w1 <- dnorm(yI[t], mu_g, sqrt(sig_g2 + sig_y2), log = TRUE) - w
    w1 <- exp(w1 - max(w1))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    lpars <- lpars[k, ]
    sig_g2 <- sig_g2[k]
    sig_y2 <- sig_y2[k]
    S_t <- S_t[k]
    E_t <- E_t[k]
    I_t <- I_t[k]
    R_t <- R_t[k]
    b0_t0 <- b0_t0[k]
    d0_t0 <- d0_t0[k]
    # Offline sampling fixed parameters
    b0_t <- b0_t0 + (yI[t] - g_t) ^ 2 / 2
    sig_g2 <- 1 / rgamma(M, a0 + t / 2, b0_t)
    d0_t <- d0_t0 + (g_t - mu_g) ^ 2 / 2
    sig_y2 <- 1 / rgamma(M, c0 + t / 2, d0_t)
    # Storage
    pars <- exp(lpars)
    qs[1, , t] <- quantile(S_t, c(.05, .5, .95))
    qs[2, , t] <- quantile(I_t, c(.05, .5, .95))
    qs[3, , t] <- quantile(E_t, c(.05, .5, .95))
    qs[4, , t] <- quantile(R_t, c(.05, .5, .95))
    qs[5, , t] <- quantile(pars[, 2], c(.05, .5, .95))
    qs[6, , t] <- quantile(pars[, 1], c(.05, .5, .95))
    qs[7, , t] <- quantile(pars[, 3], c(.05, .5, .95))
    qs[8, , t] <- quantile(sqrt(sig_g2), c(.05, .5, .95))
    qs[9, , t] <- quantile(sqrt(sig_y2), c(.05, .5, .95))
    qs[10, , t] <- quantile(g_t, c(.05, .5, .95))
  }
  if (ahead == 0) return(qs)
  if (ahead == 1) return(list(qs = qs, loglike = loglike, pred = pred))
}

PFlike <- function(y, yI, alpha, beta, gamma, a0, b0, c0, d0, M) {
  y <- y / 100
  n <- length(yI)
  sig_g2 <- 1 / rgamma(M, a0, b0)
  sig_y2 <- 1 / rgamma(M, c0, d0)
  S_t <- rep(1 - y[1], M)
  E_t <- rep(0, M)
  I_t <- rep(y[1], M)
  R_t <- rep(0, M)
  b0_t <- rep(b0, M)
  d0_t <- rep(d0, M)
  loglike <- rep(0, n)
  for (t in 1:n) {
    print(t)
    # Resample (S, E, I, R)
    mu_g <- -gamma + alpha * E_t / I_t
    w <- dnorm(yI[t], mu_g, sqrt(sig_g2 + sig_y2), log = TRUE)
    # loglike[t] <- log(mean(exp(w)))
    # original formula, y_t is reduced when calculate log Bayes factor
    loglike[t] <- log(mean(exp(w) / y[t])) 
    w1 <- exp(w - max(w))
    k <- sample(1:M, size = M, replace = TRUE, prob = w1)
    sig_g2 <- sig_g2[k]
    sig_y2 <- sig_y2[k]
    S_t0 <- S_t[k]
    E_t0 <- E_t[k]
    I_t0 <- I_t[k]
    R_t0 <- R_t[k]
    b0_t0 <- b0_t[k]
    d0_t0 <- d0_t[k]
    # Update I
    mu_g <- mu_g[k]
    B <- 1 / (1 / sig_g2 + 1 / sig_y2)
    b <- B * (yI[t] / sig_g2 + mu_g / sig_y2)
    g_t <- rnorm(M, b, sqrt(B))
    I_t <- I_t0 * (1 + g_t)
    # Update (E, R, S)
    E_t <- beta * I_t0 * S_t0 + (1 - alpha) * E_t0
    R_t <- R_t0 + gamma * I_t0
    S_t <- 1 - I_t - R_t - E_t
    # Offline sampling fixed parameters
    b0_t <- b0_t0 + (yI[t] - g_t) ^ 2 / 2
    sig_g2 <- 1 / rgamma(M, a0 + t / 2, b0_t)
    d0_t <- d0_t0 + (g_t - mu_g) ^ 2 / 2
    sig_y2 <- 1 / rgamma(M, c0 + t / 2, d0_t)
  }
  return(loglike)
}

ar1plusnoise <- function(y, yI, q0, Q0, a0, b0, c0, d0, mu, phi, V, W, g0) {
  y <- y / 100
  n <- length(yI)
  M <- length(g0)
  # sufficient statisitics, s_t  = c(a_t, b_t, c_t, d_t, q_t, Q_t)
  a_t <- rep(a0, M) # a_t
  b_t <- rep(b0, M) # b_t
  c_t <- rep(c0, M) # c_t
  d_t <- rep(d0, M) # d_t
  q_t <- t(replicate(M, q0)) # q_t, namely q_t = (E(mu), E(phi))
  # Q_t^{-1}, inverse of Q_t is a symmetric matrix, 
  # only upper triangular matrix elements are stored.
  invQ0 <- solve(Q0)
  invQ_t <- matrix(0, M, 3)
  invQ_t[, 1] <- invQ0[1, 1]
  invQ_t[, 2] <- invQ0[1, 2]
  invQ_t[, 3] <- invQ0[2, 2]
  Qq_t <- t(replicate(M, q0 / diag(Q0))) # Q_t^{-1}q_t = inverse(Q_t) * q_t
  g_t <- g0
  I <- rep(y[1], M)
  qs <- array(0, c(6, 3, n)) # quantiles of 6 parameters of interest
  loglike <- rep(0, n)
  pred <- matrix(0, M, n)
  for (t in 1:n) {
    print(t)
    g_t1 <- mu + phi * g_t # g_t = mu + phi * g_{t-1}
    pred[, t] <- I * (1 + g_t1)
    w <- dnorm(yI[t], g_t1, sqrt(V + W))
    loglike[t] <- log(mean(w / y[t]))
    k <- sample(1:M, size = M, replace = TRUE, prob = w)
    w <- w[k]
    a_t0 <- a_t[k]
    b_t0 <- b_t[k]
    c_t0 <- c_t[k]
    d_t0 <- d_t[k]
    q_t0 <- q_t[k, ]
    invQ_t0 <- invQ_t[k, ]
    Qq_t0 <- Qq_t[k, ]
    g_t0 <- g_t[k]
    g_t1 <- g_t1[k]
    mu <- mu[k]
    phi <- phi[k]
    W <- W[k]
    V <- V[k]
    C_t  <- 1 / (1 / V + 1 / W)
    m_t <- C_t * (yI[t] / V + g_t1 / W)
    g_t <- rnorm(M, m_t, sqrt(C_t)) # (g_t|s_t^k, epison) ~ N(m_t, C_t)
    # x_t * x_t' = (1, g_{t-1}) * (1, g_{t-1})'
    invQ_t[, 1] <- invQ_t0[, 1] + 1
    invQ_t[, 2] <- invQ_t0[, 2] + g_t0
    invQ_t[, 3] <- invQ_t0[, 3] + g_t0 ^ 2
    x_t <- cbind(1, g_t0)
    Qq_t <- Qq_t0 + g_t * x_t
    m <- invQ_t[, 1] * invQ_t[, 3] - invQ_t[, 2] ^ 2 # determinant of Q_t^{-1}
    q_t[, 1] <- (invQ_t[, 3] * Qq_t[, 1] - invQ_t[, 2] * Qq_t[, 2]) / m
    q_t[, 2] <- (invQ_t[, 1] * Qq_t[, 2] - invQ_t[, 2] * Qq_t[, 1]) / m
    c_t <- c_t0 + 1 / 2 # c_t = c_{t-1} + 1/2
    # y_t is replaced by g_t to avoid NA
    d_t <- d_t0 + (g_t - rowSums(q_t * x_t)) * g_t / 2 + 
      rowMeans((q_t0 - q_t) * Qq_t0) # d_t = d_{t-1} + ...
    W <- 1 / rgamma(M, c_t, d_t) # (W|g^t, X^t) ~ IG(c_t, d_t)
    a_t <- a_t0 + 1 / 2 # a_t = a_{t-1} + 1/2
    # b_t <- b_t0 + (yI[t] - g_t) ^ 2 # equation from Dukic's code
    b_t <- b_t0 + (yI[t] - g_t) ^ 2 / 2 # b_t = b_{t-1} + (y_t - g_t)^2 / 2
    V <- 1 / rgamma(M, a_t, b_t) # (V|y^t, g^t) ~ IG(a_t, b_t)
    std <- sqrt(W / m)
    norm1 <- rnorm(M, 0, std)
    norm2 <- rnorm(M, 0, std)
    mu <- q_t[, 1] + sqrt(invQ_t[, 3]) * norm1
    phi <- q_t[, 2] - invQ_t[, 2] / sqrt(invQ_t[, 3]) * norm1 +
      sqrt(invQ_t[, 1] - invQ_t[, 2] ^ 2 / invQ_t[, 3]) * norm2
    I <- I * (1 + g_t)
    qs[1, , t] <- quantile(mu, c(.05, .5, .95))
    qs[2, , t] <- quantile(phi, c(.05, .5, .95))
    qs[3, , t] <- quantile(W, c(.05, .5, .95))
    qs[4, , t] <- quantile(V, c(.05, .5, .95))
    qs[5, , t] <- quantile(g_t, c(.05, .5, .95))
    qs[6, , t] <- quantile(I, c(.05, .5, .95))
  }
  return(list(qs = qs, loglike = loglike, pred = pred))
}
