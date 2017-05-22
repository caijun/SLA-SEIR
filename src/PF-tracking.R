rm(list = ls())
setwd("~/Documents/github/SLA-SEIR")

source("src/PF-subroutines.R")

# Data set
# --------
data <- read.table("Google-CDC-US-Feb2010.txt", header = TRUE)
y <- data[1:36, 5] # y is ILI percentage (%) of 2003/2004 season
# y <- data[245:314, 5] # y is ILI percentage (%) of 2008/2009 season

# Prior hyperparameters
# ---------------------
a0 <- 1.1
b0 <- (a0 - 1) * 0.05 # 0.005, evolution variance
c0 <- 1.1
d0 <- (c0 - 1) * 0.5 # 0.05, observation variance
ma <- 2.0
sda <- 0.5 # sd of latency parameter, \alpha
mb <- 1.5
sdb <- 0.5 # sd of transmission parameter, \beta
mg <- 1.0
sdg <- 0.5 # sd of recovery parameter, \gamma

# Sequential learning
# -------------------
M <- 100000 # number of particles used at each iteration
a <- 0.99 # Liu-West shrinkage factor
set.seed(1321654)
n <- length(y)
yI <- diff(y) / y[1:(n - 1)] # the observed growth rate of the infectious population
qs <- PF(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a)$qs
# loglike1 <- PFlike(y, yI, 2, 1.25, 1, a0, b0, c0, d0, M)
# loglike2 <- PFlike(y, yI, 2, 2.20, 1, a0, b0, c0, d0, M)
# The log Bayes factor calculated by PF is more reasonable than PFlike that is 
# used in  Dukic, Lopes and Polson (2012).
loglike1 <- PF(y, yI, a0, b0, c0, d0, 2, sda, 1.25, sdb, 1, sdg, M, a)$loglike
loglike2 <- PF(y, yI, a0, b0, c0, d0, 2, sda, 2.20, sdb, 1, sdg, M, a)$loglike
clbf <- cumsum(loglike1 - loglike2) # cummulative log Bayesian Factor

names <- c("Transmission", "Latency", "Recovery", "Obs St Dev", "Evo St Dev")
# only from the 2nd week, the observed growth rate of the infectious population
# can be calculated.
# tick index and label for 2003/2004 season
tickind <- trunc(seq(2, 36, length = 4)) - 1
ticklab <- c("10/5/03", "12/21/03", "3/7/04", "5/30/04")
# # tick index and label for 2008/2009 season
# tickind <- trunc(seq(246, 314, length = 4)) - 244 - 1
# ticklab <- c("6/8/08", "11/9/08", "4/19/09", "9/27/09")

# Figure 5 of Dukic, Lopes and Polson (2012)
# pdf(file = "fig/Figure5.pdf", width = 20, height = 10)
png(file = "fig/Figure5.png", width = 20, height = 10, units = "in", res = 300)
# # Figure 6 of Dukic, Lopes and Polson (2012)
# # pdf(file = "fig/Figure6.pdf", width = 20, height = 10)
# png(file = "fig/Figure6.png", width = 20, height = 10, units = "in", res = 300)
par(mar = c(5, 5, 4, 3), cex.lab = 3, cex.axis = 1.6, cex.main = 2.5,
    mfrow = c(2, 4))

plot(qs[1, 2, ], ylim = range(qs[1, , ]), yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T, type = "l", lwd = 2)
title(xlab = "", ylab = "Proportion", main = "S")
axis(1, at = tickind, lab = ticklab)
lines(qs[1, 1, ], lwd = 2, col = grey(0.5))
lines(qs[1, 3, ], lwd = 2, col = grey(0.5))

plot(qs[2, 2, ], ylim = range(qs[2, , ]), yaxs = "i", xaxs = "i", xaxt = "n",
     xlab = "", ylab = "", axes = T, type = "l", lwd = 2)
title(xlab = "", ylab = "Proportion", main = "I (and observations)")
axis(1, at = tickind, lab = ticklab)
lines(qs[2, 1, ], lwd = 2, col = grey(0.5))
lines(qs[2, 3, ], lwd = 2, col = grey(0.5))
lines(qs[2, 2, ], lwd = 2)
points(y[2:n] / 100, pch = 16)

for (i in 5:7) {
  miny <- min(qs[i, , ])
  maxy <- max(qs[i, , ])
  plot(qs[i, 2, ], ylim = c(miny, maxy), yaxs = "i", xaxs = "i", xaxt = "n",
       xlab = "", ylab = "", axes = T, type = "l", lwd = 2)
  title(xlab = "", ylab = "", main = names[i - 4])
  axis(1, at = tickind, lab = ticklab)
  lines(qs[i, 1, ], lwd = 2, col = grey(0.5))
  lines(qs[i, 3, ], lwd = 2, col = grey(0.5))
}

for (i in 8:9) {
  plot(qs[i, 2, ], ylim = range(qs[i, , ]), yaxs = "i", xaxs = "i", xaxt = "n",
       xlab = "", ylab = "", axes = T, type = "l", lwd = 2)
  title(xlab = "", ylab = "", main = names[i - 4])
  axis(1, at = tickind, lab = ticklab)
  lines(qs[i, 1, ], lwd = 2, col = grey(0.5))
  lines(qs[i, 3, ], lwd = 2, col = grey(0.5))
}
plot(clbf, yaxs = "i", xaxs = "i", xaxt = "n", xlab = "", ylab = "", axes = T,
     type = "l", lwd = 2)
title(xlab = "", ylab = "", main = "Log Bayes factor")
axis(1, at = tickind, lab = ticklab)
abline(h = 0, lwd = 2, lty = 2, col = grey(0.5))

dev.off()
