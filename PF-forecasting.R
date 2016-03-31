rm(list = ls())
setwd("~/Documents/github/SLA-SEIR")

source("PF-subroutines.R")

# Data set
# --------
data <- read.table("Google-CDC-US-Feb2010.txt", header = TRUE)
y1 <- data[1:36, 5] # 2003/2004 season
y2 <- data[245:314, 5] # 2008/2009 season

# SEIR for growth rate
# --------------------
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
M <- 100000
a <- 0.99
n1 <- length(y1)
n2  <- length(y2)
y1I <- diff(y1) / y1[1:(n1 - 1)]
y2I <- diff(y2) / y2[1:(n2 - 1)]
n1 <- n1 - 1
n2 <- n2 - 1
pf1 <- PF1(y1, y1I, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a)
pf2 <- PF1(y2, y2I, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a)
loglike.pf1 <- pf1$loglike
loglike.pf2 <- pf2$loglike
pred.pf1 <- pf1$pred
pred.pf2 <- pf2$pred
qs.pf1 <- pf1$qs
qs.pf2 <- pf2$qs

# AR1 + noise for growth rate
# ---------------------------
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
ssm1 <- ar1plusnoise(y1, y1I, q0, Q0, a0, b0, c0, d0, mu, phi, V, W, g0)
ssm2 <- ar1plusnoise(y2, y2I, q0, Q0, a0, b0, c0, d0, mu, phi, V, W, g0)
loglike.ssm1 <- ssm1$loglike
loglike.ssm2 <- ssm2$loglike
pred.ssm1 <- ssm1$pred
pred.ssm2 <- ssm2$pred
ssm1 <- ssm1$qs
ssm2 <- ssm2$qs

# Comparison
# ----------
cs.ssm1 <- cumsum(loglike.ssm1)
cs.ssm2 <- cumsum(loglike.ssm2)
cs.pf1 <- cumsum(loglike.pf1)
cs.pf2 <- cumsum(loglike.pf2)

mpred.pf1 <- t(apply(pred.pf1, 2, quantile, c(.05, .5, .95)))
mpred.pf2 <- t(apply(pred.pf2, 2, quantile, c(.05, .5, .95)))
mpred.ssm1 <- t(apply(pred.ssm1, 1, quantile, c(.05, .5, .95)))
mpred.ssm2 <- t(apply(pred.ssm2, 1, quantile, c(.05, .5, .95)))
rmse1 <- var(y1[1:n1] - 100 * mpred.ssm1[, 2]) / var(y1[1:n1] - 100 * mpred.pf1[, 2])
rmse2 <- var(y2[1:n2] - 100 * mpred.ssm2[, 2]) / var(y2[1:n2] - 100 * mpred.pf2[, 2])

# only from the 2nd week, the observed growth rate of the infectious population
# can be calculated.
tickind1 <- trunc(seq(2, 36, length = 4)) - 1
ticklab1 <- c("10/5/03", "12/21/03", "3/7/04", "5/30/04")
tickind2 <- trunc(seq(246, 314, length = 4)) - 244 - 1
ticklab2 <- c("6/8/08", "11/9/08", "4/19/09", "9/27/09")

# Figure 10 of Dukic, Lopes and Polson (2012)
pdf(file = "Figure10.pdf", width = 10, height = 7)
par(cex.lab = 1.7, cex.axis = 1.3, cex.main = 1.9, mfrow = c(2, 2))

plot(y1[1:n1], ylim = c(0, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T)
title(xlab = "", ylab = "Infected %", main = "SEIR Model 2003/2004", cex.lab = 1.4)
axis(1, at = tickind1, lab = ticklab1)
lines(100 * mpred.pf1[, 1], col = grey(0.65), lwd = 2, lty = 2)
lines(100 * mpred.pf1[, 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
lines(100 * mpred.pf1[, 3], col = grey(0.65), lwd = 2, lty = 2)
points(y1[1:n1], pch = 16)

plot(y1[1:n1], ylim = c(0, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T)
title(xlab = "", ylab = "Infected %", main = "AR(1) plus noise model, 2003/2004",
      cex.lab = 1.4)
axis(1, at = tickind1, lab = ticklab1)
lines(100 * mpred.ssm1[, 1], col = grey(0.65), lwd = 2, lty = 2)
lines(100 * mpred.ssm1[, 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
lines(100 * mpred.ssm1[, 3], col = grey(0.65), lwd = 2, lty = 2)
points(y1[1:n1], pch = 16)

plot(y2[1:n2], ylim = c(0, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T)
title(xlab = "", ylab = "Infected %", main = "SEIR Model 2008/2009", cex.lab = 1.4)
axis(1, at = tickind2, lab = ticklab2)
lines(100 * mpred.pf2[, 1], col = grey(0.65), lwd = 2, lty = 2)
lines(100 * mpred.pf2[, 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
lines(100 * mpred.pf2[, 3], col = grey(0.65), lwd = 2, lty = 2)
points(y2[1:n2], pch = 16)

plot(y2[1:n2], ylim = c(0, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T)
title(xlab = "", ylab = "Infected %", main = "AR(1) plus noise model, 2008/2009",
      cex.lab = 1.4)
axis(1, at = tickind2, lab = ticklab2)
lines(100 * mpred.ssm2[, 1], col = grey(0.65), lwd = 2, lty = 2)
lines(100 * mpred.ssm2[, 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
lines(100 * mpred.ssm2[, 3], col = grey(0.65), lwd = 2, lty = 2)
points(y2[1:n2], pch = 16)

dev.off()
