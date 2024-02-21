# ------------------------------------------------------------------------#
# Adapt example shown by ?refund::fosr2s to become example for function on
# scalar 2-step regression. It uses the canadian weather data set and
# shows Effect of latitude on daily mean temperatures.
# 
# We will sample every seven days from the daily temperature because this 
# is simply a demonstration and if all 365 points are used,
# the results of all the points cover eachother on the plot
# ------------------------------------------------------------------------#

# Packages used: ----------------------------------------------------------
library(refund)     # CRAN v0.1-24
library(fda)        # CRAN v5.5.1
library(ggplot2)    # [github::tidyverse/ggplot2] v3.3.5.9000
library(tikzDevice) # CRAN v0.12.3.1 
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{upgreek}"))
# Set up a reduced, weekly data set ---------------------------------------

plots_path <- here::here("chapter-04", "figures")

# Reduced indices - sample every 365 days
weekly_ind <- seq(1, 365, by = 7)
weekly_temp <- CanadianWeather$dailyAv[weekly_ind,, 1]
# ^ 1st slice of CanadianWeather$dailyAv gives the daily temperature
weekly_argvals <- day.5[weekly_ind]
matplot(x = weekly_argvals,
        y = weekly_temp,
        type = "b",
        lty = 1,
        pch = 1)




# Set up data for a regression --------------------------------------------

tempmat <- t(weekly_temp)
                    # Create design matrix. First is a column of 1's for intercept. 
# binded to a column representing the centered lattitude:
centered_lat <- scale(CanadianWeather$coord[ , 1],
                         center = TRUE, # centred!
                         scale = FALSE)
latmat <- cbind(1, centered_lat)


# Fit model: --------------------------------------------------------------
fzmod <- fosr2s(Y = tempmat,
                X = latmat,
                argvals = weekly_argvals,
                basistype="bspline", # define basis for smoothed response
                nbasis = 25,
                norder = 4,
                pen.order = 2)

par(mfrow = c(1, 2))
for(j in 1:2) {
  plot(fzmod$raw.coef[, j], pch = "β")
  lines(fzmod$est.func[, j])
}



# See can we do this manually ---------------------------------------------

# fit pointwise regressions:
raw_est <- purrr::map_dfr(
  .x = purrr::array_tree(array = weekly_temp, margin = 1),
  .f = ~coef(lm(.x ~ centered_lat)))
names(raw_est) <- c("beta_0", "beta_1")

# Now smooth:
# Set up basis function for smoothing:
bspl_basis_fd <- create.bspline.basis(
  rangeval = c(weekly_argvals[1], weekly_argvals[length(weekly_argvals)]),
  nbasis = 25, norder = 4)
# get basis penalty as we will smooth using mgcv::gam() function and will need this:
bspl_penmat <- getbasispenalty(basisobj = bspl_basis_fd, Lfdobj = 2)
# and evaluation of basis functions at eval points:
X <- eval.basis(evalarg = weekly_argvals, basisobj = bspl_basis_fd)
# and a finer set of grid points to evaluate the smoothed functions on:
coef_argvals <- seq(min(weekly_argvals), max(weekly_argvals), length.out = 100)

# estimate basis function coefficients using gam function and use REML estimation
                    # create an fda object using the estimate coefficients
# and evaluate them on our fine grid:
smoothed_est <- purrr::map_dfc(
  .x = raw_est,
  .f = ~{
    gam_fit <- mgcv::gam(formula = .x ~ X - 1,
                       paraPen = list(X = list(bspl_penmat)),
                       method = "REML")
    fd_obj <- fd(coef = gam_fit$coefficients, basisobj = bspl_basis_fd)
    as.vector(eval.fd(evalarg = coef_argvals, fdobj = fd_obj))
  })




# Compare manual and fosr2s function results ------------------------------

par(mfrow = c(2, 2))
for(j in 1:2){
  plot(fzmod$argvals, fzmod$raw.coef[, j])
  points(weekly_argvals, raw_est[[j]], col = "red")
  
  lines(fzmod$argvals, fzmod$est.func[, j])
  lines(coef_argvals, smoothed_est[[j]], col = "red")
}

# Identical





# Make plot for publication -----------------------------------------------
raw_est$t <- weekly_argvals
smoothed_est$t <- coef_argvals
raw_est$shape <- "Pointwise coefficient estimate"
smoothed_est$col <- "Smoothed"


p <- ggplot(data = raw_est) +
  aes(x = t, y = beta_1) +
  # scale_shape_manual(values = "β") +
  # scale_shape_manual("h") +
  geom_rug(sides = "b") +
  geom_line(col = "grey4", linewidth = 0.5, lty = 1, alpha = 0.75) +
  geom_point(col = "white", size = 3) +
  geom_line(data = smoothed_est, 
            mapping = aes(colour = col),
            lwd = 1.5) +
  scale_color_manual(values = "red4", labels = "Smoothed coefficient function $\\widehat{\\beta} (t)$") +
  ylab("$\\widehat{\\beta} (t)$") +
  xlab("$t$") +
  theme_bw() +
  # geom_point(mapping = aes(shape = shape), size = 4) +
  geom_text(mapping = aes(label = "$\\upbeta$")) +
  theme(legend.title = element_blank(),
        axis.text = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        # legend.position = c(0.2, 0.875),
        legend.position = "none",
        legend.background = element_rect(colour = "lightgrey"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 13)
        )

p

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "fosr-2s.tex"),
     width = 1 * doc_width_inches, 
     height = (5/7) * doc_width_inches,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file.path(plots_path, "fosr-2s.tex"))
# ggsave(filename = here::here("Figures/fosr-2-step.png"), device = "png", width = 7, height = 5, dpi = "retina")
           
