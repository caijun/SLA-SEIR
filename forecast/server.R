shinyServer(function(input, output) {
  output$selectInput <- renderUI({
    selectInput("sel.season", label = "Season", 
                choices = list("2003/2004", "2008/2009"), 
                selected = "2003-04")
  })
  
  output$sliderInput <- renderUI({
    if (input$sel.season == "2003/2004") maxVal <- 36 else maxVal <- 70
    
    sliderInput("slider", label = "Week number",
                min = 1, max = maxVal, value = 1,
                animate = TRUE, step = 1)
  })
  
  output$plot <- renderPlot({
    if (input$sel.season == "2003/2004") load("draw2003-04.RData") else load("draw2008-09.RData")
    
    i <- input$slider
    y <- c(y[1:i], rep(NA, length(y) - i))
    mpred.pf <- rbind(mpred.pf, rep(NA, 3))
    mpred.ssm <- rbind(mpred.ssm, rep(NA, 3))
    
    par(cex.lab = 1.7, cex.axis = 1.3, cex.main = 1.9, mfrow = c(2, 1))
    # SEIR plot
    plot(y, ylim = c(-1, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
         xlab = "", ylab = "", axes = T)
    points(y, pch = 16)
    title(xlab = "", ylab = "Infected %", main = paste0("SEIR Model ", input$sel.season), 
          cex.lab = 1.4)
    axis(1, at = tickInd, lab = tickLab)
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
    # AR(1) plot
    plot(y, ylim = c(-1, 20), pch = 15, yaxs = "i", xaxs = "i", xaxt = "n",
         xlab = "", ylab = "", axes = T)
    points(y, pch = 16)
    title(xlab = "", ylab = "Infected %", main = paste0("AR(1) plus noise model, ", input$sel.season), 
          cex.lab = 1.4)
    axis(1, at = tickInd, lab = tickLab)
    points(i + 1, 100 * mpred.ssm[i, 2], pch = 16, col = "red") # prediction for i + 1 week
    if (i > 1) {
      lines(2:i, 100 * mpred.ssm[1:(i - 1), 1], col = grey(0.65), lwd = 2, lty = 2)
      lines(2:i, 100 * mpred.ssm[1:(i - 1), 2], col = grey(0.65), lwd = 2, type = "b", pch = 16)
      lines(2:i, 100 * mpred.ssm[1:(i - 1), 3], col = grey(0.65), lwd = 2, lty = 2)
      lines(i:(i + 1), 100 * mpred.ssm[(i - 1):i, 1], col = "red", lwd = 2, lty = 2)
      lines(i:(i + 1), 100 * mpred.ssm[(i - 1):i, 2], col = "red", lwd = 2)
      lines(i:(i + 1), 100 * mpred.ssm[(i - 1):i, 3], col = "red", lwd = 2, lty = 2)
    }
    abline(h = 0, lty = 2)
    abline(v = i, col = "blue")
    legend("topright", c("Observation", "Previous Forecast", "Current Forecast"), 
           col = c("black", grey(0.65), "red"), pch = c(15, 16, 16))
  }, height = 600, width = 500)
})