setwd("~/Documents/github/SLA-SEIR")

source("PF-subroutines.R")

drawData <- function(season) {
  data <- read.table("Google-CDC-US-Feb2010.txt", header = TRUE)
  dateLabel <- format(as.Date(data$DateLabel), format = "%m/%d/%y")
  
  if (season == "2003-04") {
    from <- 1
    to <- 36
  }
  if (season == "2008-09") {
    from <- 245
    to <- 314
  }
  
  y <- data[from:to, 5]
  tickInd <- trunc(seq(from, to, length = 4))
  tickLab <- dateLabel[tickInd]
  tickInd <- tickInd - (from - 1)
  n <- length(y)
  yI <- diff(y) / y[1:(n - 1)] # the observed growth rate of the infectious population
  
  cat("SEIR model starts...\n")
  set.seed(1321654)
  a0 <- 1.1
  b0 <- (a0 - 1) * 0.05 # 0.005, evolution variance
  c0 <- 1.1
  d0 <- (c0 - 1) * 0.5 # 0.05, observation variance
  ma <- 2.0
  sda <- 0.5
  mb <- 1.1
  sdb <- 0.5
  mg <- 1.0
  sdg <- 0.5
  M <- 10000
  a <- 0.99
  
  pf <- PF(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a, ahead = 1)
  pred.pf <- pf$pred
  mpred.pf <- t(apply(pred.pf, 2, quantile, c(.05, .5, .95)))
  cat("SEIR model ends.\n")
  
  cat("AR(1) model starts...\n")
  set.seed(13121654)
  a0 <- 1.1
  b0 <- (a0 - 1) * 0.5 # 0.05, observation variance
  c0 <- 1.1
  d0 <- (c0 - 1) * 0.05 # 0.005, evolution variance
  q0 <- c(0, 1)
  Q0 <- diag(10, 2)
  m0 <- 0
  C0 <- 1
  V <- 1 / rgamma(M, a0, b0)
  W <- 1 / rgamma(M, c0, d0)
  mu <- rnorm(M, q0[1], sqrt(W * Q0[1, 1]))
  phi <- rnorm(M, q0[2], sqrt(W * Q0[2, 2]))
  g0 <- rnorm(M, m0, sqrt(C0))
  
  ssm <- ar1plusnoise(y, yI, q0, Q0, a0, b0, c0, d0, mu, phi, V, W, g0)
  pred.ssm <- ssm$pred
  mpred.ssm <- t(apply(pred.ssm, 2, quantile, c(.05, .5, .95)))
  cat("AR(1) model ends.\n")
  
  rmse <- var(y[2:n] - 100 * mpred.ssm[, 2]) / var(y[2:n] - 100 * mpred.pf[, 2])
  print(rmse)
  rdata <- paste0("forecast/draw", season, ".RData")
  save(y, tickInd, tickLab, mpred.pf, mpred.ssm, rmse, file = rdata)
}

drawData("2003-04")
drawData("2008-09")