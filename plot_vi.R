#' Plot Variational Inference Result
#'
#' Plots the data points and the fitted Gaussian ellipsoid (2D or 3D).
#'
#' @param samples A numeric matrix of samples (rows are observations).
#' @param result The output list from \code{variational_inference()}.
#' @export
plot_vi <- function(samples, result) {
  dim <- ncol(samples)
  if (dim == 2) {
    library(ggplot2)
    library(ellipse)
    samples_df <- data.frame(x = samples[, 1], y = samples[, 2])
    conf_ellipse <- data.frame(ellipse(result$Sigma, centre = result$mu, level = 0.95))
    ggplot(samples_df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3) +
      geom_path(data = conf_ellipse, aes(x = x, y = y), color = "red", size = 1.2) +
      geom_point(aes(x = result$mu[1], y = result$mu[2]), color = "blue", size = 3) +
      theme_minimal() +
      ggtitle("2D Variational Inference Result") +
      xlab("X1") + ylab("X2")
  } else if (dim == 3) {
    library(rgl)
    library(car)
    rgl.open()
    plot3d(samples[,1], samples[,2], samples[,3], col = "blue", size = 3, alpha = 0.5,
           xlab = "X1", ylab = "X2", zlab = "X3")
    ellipsoid <- ellipse3d(result$Sigma, centre = result$mu, level = 0.95)
    shade3d(ellipsoid, col = "red", alpha = 0.3)
    points3d(result$mu[1], result$mu[2], result$mu[3], col = "black", size = 10)
    title3d("3D Variational Inference Result", "", "", "", "")
    rglwidget()
  } else {
    stop("Plotting only supported for 2D or 3D.")
  }
}