# 1. Load the aggregated results from the supercomputer
agg_results <- readRDS("Simulation_Results")

# 2. Extract the columns for plotting
x_plot_cic <- -agg_results$Mean_CIC
x_plot_aic <- -agg_results$Mean_AIC2ML

# 3. Create the empty plot and plot the CIC points (Black Triangles)
plot(x_plot_cic, agg_results$Scaled_MSE,
     pch = 2, col = "black", cex = 1.2,
     xlim = range(c(x_plot_cic, x_plot_aic), na.rm = TRUE),
     ylim = range(agg_results$Scaled_MSE, na.rm = TRUE),
     xlab = "Model selection criterion value",
     ylab = expression(10^3 %*% MSE(hat(P)))
     #,main = "Estimation Performance vs Criteria"
     )

# 4. Add the AIC2ML points (Red Crosses)
points(x_plot_aic, agg_results$Scaled_MSE,
       pch = 4, col = "indianred4", cex = 1.2)

# 5. Add the trend lines
abline(lm(agg_results$Scaled_MSE ~ x_plot_cic), col = "black", lty = 1, lwd = 1.5)
abline(lm(agg_results$Scaled_MSE ~ x_plot_aic), col = "darkgrey", lty = 2, lwd = 1.5)

# 6. Add the legend
legend("topright",
       legend = c("CIC", expression(AIC[2*ML]), "CIC (trend)", expression(AIC[2*ML]~"(trend)")),
       pch = c(2, 4, NA, NA),
       col = c("black", "indianred4", "black", "darkgrey"),
       lty = c(NA, NA, 1, 2),
       lwd = c(NA, NA, 1.5, 1.5),
       pt.cex = 1.2, bty = "o")

# 2. Plotting Figure 2
plot(agg_results$Mean_pen_CIC, agg_results$Scaled_MSE,
     pch = 2, col = "black", cex = 1.2,
     xlim = range(c(agg_results$Mean_pen_CIC, agg_results$Mean_pen_AIC), na.rm = TRUE),
     ylim = range(agg_results$Scaled_MSE, na.rm = TRUE),
     xlab = expression("Bias correction term"),
     ylab = expression(10^3 %*% MSE(hat(P))),
     main = "MSE vs Bias Correction Term")

# Add the vertical lines of AIC points
points(agg_results$Mean_pen_AIC, agg_results$Scaled_MSE,
       pch = 4, col = "indianred4", cex = 1.2)

# Add the single CIC trend line
abline(lm(agg_results$Scaled_MSE ~ agg_results$Mean_pen_CIC), col = "black", lty = 1, lwd = 1.5)

legend("bottomright",
       legend = c("CIC", expression(AIC[2*ML]), "CIC (trend)"),
       pch = c(2, 4, NA),
       col = c("black", "indianred4", "black"),
       lty = c(NA, NA, 1),
       lwd = c(NA, NA, 1.5),
       pt.cex = 1.2, bty = "o")

