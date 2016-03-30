# setwd("~/Documents/github/SLA-SEIR")
# 
# source("PF-subroutines.R")
# 
# # Data set
# # --------
# data <- read.table("Google-CDC-US-Feb2010.txt", header = TRUE)
# y <- data[1:36, 5] # y is ILI percentage (%) of 2003-04 season
# 
# # Prior hyperparameters
# # ---------------------
# a0 <- 1.1
# b0 <- (a0 - 1) * 0.05 # 0.005, evolution variance
# c0 <- 1.1
# d0 <- (c0 - 1) * 0.5 # 0.05, observation variance
# ma <- 2.0
# sda <- 0.5 # sd of latency parameter, \alpha
# mb <- 1.5
# sdb <- 0.5 # sd of transmission parameter, \beta
# mg <- 1.0
# sdg <- 0.5 # sd of recovery parameter, \gamma
# 
# # Sequential learning
# # -------------------
# M <- 100000 # number of particles used at each iteration
# a <- 0.99 # Liu-West shrinkage factor
# set.seed(1321654)
# n <- length(y)
# yI <- diff(y) / y[1:(n - 1)] # the observed growth rate of the infectious population
# pf <- PF1(y, yI, a0, b0, c0, d0, ma, sda, mb, sdb, mg, sdg, M, a)
# pred.pf <- pf$pred
# mpred.pf <- t(apply(pred.pf, 2, quantile, c(.05, .5, .95)))
# save(y, mpred.pf, file = "forecast/draw.RData")

load("draw.RData")

tickind1 <- trunc(seq(1, 36, length = 4))
ticklab1 <- c("9/28/03", "12/14/03", "3/7/04", "5/30/04")

shinyServer(function(input, output) {
  output$plot <- renderPlot({
    i <- input$slider
    y <- c(y[1:i], rep(NA, 36 - i))
    mpred.pf <- rbind(mpred.pf, rep(NA, 3))
    plot(y, ylim = c(-1, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
         xlab = "", ylab = "", axes = T)
    points(y, pch = 16)
    title(xlab = "", ylab = "Infected %", main = paste0("SEIR Model ", input$sel.season), 
          cex.lab = 1.4)
    axis(1, at = tickind1, lab = ticklab1)
    points(i + 1, 100 * mpred.pf[i, 2], pch = 16, col = "red") # prediction for i + 1 week
    if (i > 1) {
      lines(2:i, 100 * mpred.pf[1:(i - 1), 1], col = grey(0.65), lwd = 2, lty = 2)
      lines(2:i, 100 * mpred.pf[1:(i - 1), 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
      lines(2:i, 100 * mpred.pf[1:(i - 1), 3], col = grey(0.65), lwd = 2, lty = 2)
      lines(i:(i + 1), 100 * mpred.pf[(i - 1):i, 1], col = "red", lwd = 2, lty = 2)
      lines(i:(i + 1), 100 * mpred.pf[(i - 1):i, 2], col = "red", lwd = 2)
      lines(i:(i + 1), 100 * mpred.pf[(i - 1):i, 3], col = "red", lwd = 2, lty = 2)
    }
    abline(h = 0, lty = 2)
    abline(v = i, col = "blue")
    legend("topright", c("Observation", "Previous Forecast", "Current Forecast"), 
           col = c("black", grey(0.65), "red"), pch = c(15, 16, 16))
  })
})