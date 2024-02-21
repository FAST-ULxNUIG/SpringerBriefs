create_functional_boxplot_gg <- function(time_grid = 0:100, fd_obj,
                                         ribbon_fill = "#cd2fbc",
                                         fence_colour = "#50abea",
                                         outlier_colour = "red",
                                         outlier_linetype = 2,
                                         linewidth = 1,
                                         outlier_linewidth = 0.75) {
  # WRAPPER function to make ggplot figure from fda::fboxplot
  require(fda)
  require(ggplot2)
  require(data.table)
  
  fda_boxplot <- fda::boxplot.fd(x = fd_obj)
  fd_obj_on_grid <- eval.fd(time_grid, fdobj = fd_obj)
  
  # median
  median_curve_ind <- fda_boxplot$medcurve
  median_curve <- fd_obj_on_grid[,median_curve_ind]
  
  # quartiles (middle 50 envelope)
  n <- ncol(fd_obj_on_grid)
  n_over_2 <- ceiling(n/2) # use in case of odd number
  
  depth <- fda_boxplot$depth
  top_n_over_2_inds <- seq_len(n)[order(-depth)][seq_len(n_over_2)]
  top_n_over_2_fd_grid <- fd_obj_on_grid[, top_n_over_2_inds]
  
  max_50_envelope <- apply(top_n_over_2_fd_grid, MARGIN = 1, max)
  min_50_envelope <- apply(top_n_over_2_fd_grid, MARGIN = 1, min)
  
  # Interquartile Range #a nd Fences
  IQR <- max_50_envelope - min_50_envelope
  min_envelope_cutoff <- min_50_envelope - 1.5 * IQR
  max_envelope_cutoff <- max_50_envelope + 1.5 * IQR
  
  is_not_outlier <- apply(fd_obj_on_grid, 
                          2, 
                          function(x) {
    all((x > min_envelope_cutoff) & (x < max_envelope_cutoff))
  })
  
  lower_fence <- apply(fd_obj_on_grid[, is_not_outlier],
                                         1,
                                         min)
  upper_fence <- apply(fd_obj_on_grid[, is_not_outlier],
                                         1,
                                         max)
  
  if(sum(!is_not_outlier) > 0) {
    n_outlier <- sum(!is_not_outlier)
    outlier_df <- cbind(time_grid, fd_obj_on_grid[, !is_not_outlier])
    colnames(outlier_df)[-1] <- paste0("outlier_", seq_len(n_outlier))
    outlier_df <- as.data.table(outlier_df)
    outlier_df <- melt.data.table(outlier_df, 
                    id.vars = "time_grid",
                    measure.vars = paste0("outlier_", seq_len(n_outlier)),
                    variable.name = "curve_ind",
                    variable.factor = TRUE,
                    value.factor = FALSE, 
                    value.name = "curve_value"
                    )
  }

  
  boxplot_df <- data.frame(
    time_grid = time_grid,
    median_curve = median_curve,
    max_50_envelope = max_50_envelope,
    min_50_envelope = min_50_envelope,
    lower_fence = lower_fence,
    upper_fence = upper_fence
  )
  
  middle_time_point <- time_grid[floor(median(seq_along(time_grid)))]
  middle_time_df <- boxplot_df[boxplot_df$time_grid == middle_time_point,]
  middle_time_df <- data.frame(time_grid = rep(middle_time_point, 2), 
                               middle = c(middle_time_df$lower_fence, middle_time_df$upper_fence))
  
  if(sum(!is_not_outlier) > 0) {
    p <- ggplot(data = boxplot_df) +
      aes(x = time_grid) +
      geom_line(data = middle_time_df, aes(y = middle), linewidth=linewidth, 
                colour = fence_colour) +
      geom_ribbon(aes(ymin = min_50_envelope, ymax = max_50_envelope), 
                  fill = ribbon_fill, 
                  linewidth=linewidth, 
                  colour = fence_colour) +
      geom_line(aes(y = median_curve), linewidth = linewidth) +
      geom_line(aes(y = lower_fence),
                linewidth=linewidth, 
                colour = fence_colour) +
      geom_line(aes(y = upper_fence),
                linewidth=linewidth, 
                colour = fence_colour) + 
      geom_line(data = outlier_df, 
                aes(group = curve_ind, y = curve_value), 
                lty = outlier_linetype,
                linewidth = outlier_linewidth,
                colour = outlier_colour)
  } else if(sum(!is_not_outlier) == 0) {
    p <- ggplot(data = boxplot_df) +
      aes(x = time_grid) +
      geom_line(data = middle_time_df, aes(y = middle), linewidth=linewidth, 
                colour = fence_colour) +
      geom_ribbon(aes(ymin = min_50_envelope, ymax = max_50_envelope), 
                  fill = ribbon_fill, 
                  linewidth=linewidth, 
                  colour = fence_colour) +
      geom_line(aes(y = median_curve), linewidth = linewidth) +
      geom_line(aes(y = lower_fence),
                linewidth=linewidth, 
                colour = fence_colour) +
      geom_line(aes(y = upper_fence),
                linewidth=linewidth, 
                colour = fence_colour)
  }
  
  if(sum(!is_not_outlier) == 0) outlier_df <- data.frame(time_grid)
  list(ggplot = p, outlier_df = outlier_df, boxplot_df = boxplot_df)
  
}




