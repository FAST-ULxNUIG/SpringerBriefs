# ======================================================================= #
# Plot comparing the basis function expansions of Fourier and BSplines #
# Demo on the gait data set for smoothing #
# ======================================================================= #

# 1) Import Packages ------------------------------------------------------
library(fda)        # CRAN v5.5.1 
library(tidyverse)  # CRAN v1.3.1 
library(ggpubr)     # CRAN v0.4.0 
library(Polychrome) # CRAN v1.5.1 
library(tikzDevice) # CRAN v0.12.3.1

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-02", "figures")
source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12))



# 3) Data  ----------------------------------------------------------

# use the gait data set:
data('gait')
dimnames(gait)
dim(gait)
# 20 x 29 x 2 array
# rows are time points
gait_time <- as.numeric(dimnames(gait)[[1]])
# choose a single observation:
gait_trial <- gait[, 1 , 2]

# 4) Set up bases ---------------------------------------------------------
# 9 (nbasis = 9) cubic (norder = 4) bspline basis functions 
# defined on the interval of the gait data 
# with equally spaced knots
bspline_9 <- create.bspline.basis(rangeval = range(gait_time), nbasis = 9, norder = 4)
plot(bspline_9, knots = F, lty = 1)

# 9 Fourier basis functions:
fourier_9 <- create.fourier.basis(rangeval = range(gait_time), nbasis = 9)
plot(fourier_9, lty = 1)


# 5) Smoothing ------------------------------------------------------------

# smooth on bspline basis (least-squares)
smooth_bs <- smooth.basis(argvals = gait_time, y = gait_trial, fdParobj = bspline_9)
smooth_f <- smooth.basis(argvals = gait_time, y = gait_trial, fdParobj = fourier_9)



# 6) Evaluate Basis and Data Wrangling --------------------------------------

# sequence of time points to evaluate data at 
evaluation_sequence <- seq(min(gait_time), max(gait_time), length.out = 501)

# evaluate the bases at these points
bspl_eval <- eval.basis(evalarg = evaluation_sequence, basisobj = bspline_9)
fourier_eval <- eval.basis(evalarg = evaluation_sequence, basisobj = fourier_9)

# weight them according to their coefficients (from smoothing)
bs_weighted <- bspl_eval %*%  diag(as.vector(smooth_bs$fd$coefs)) %>% as.data.frame
f_weighted <- fourier_eval %*%  diag(as.vector(smooth_f$fd$coefs)) %>% as.data.frame

# name columns of matrices
colnames(bs_weighted) <- colnames(f_weighted) <- paste0("fun_", 1:9)

# calculate the sum of all weighted basis functions
# this is just the fitted function in each case (i.e. the weighted sum)
bs_weighted$total <- rowSums(bs_weighted)
f_weighted$total <- rowSums(f_weighted)

# data format: wide to long
bs_lng <- data.frame(time = evaluation_sequence, bs_weighted) %>%
  gather(-time, key = "fun", value = "val_bs")

f_lng <- data.frame(time = evaluation_sequence, f_weighted) %>%
  gather(-time, key = "fun", value = "val_f")

# and combine:
basis_df <- inner_join(bs_lng, f_lng, by = c("time", "fun")) %>%
  gather(- time, - fun, key = "basis_type", value = "val") %>%
  mutate(basis_type = case_when(basis_type == "val_bs" ~ "B-spline",
                                basis_type == "val_f" ~ "Fourier"))

# and get data frame of original data to plot points
original_data <- data.frame(time = gait_time, val = gait_trial)

# 7) Plot --------------------------------------------------------------------

# Get Palette for colors of basis functions
set.seed(1996) # randomised, so set seed:
my_pallette_39 <- createPalette(N = 39, seedcolors = c("#ff0000", "#00ff00", "#0000ff")) %>%
  as.vector

p3 <- basis_df %>%
  ggplot() +
  aes(x = time, y = val) +
  facet_wrap( ~ basis_type) +
  geom_line(mapping = aes(group = fun, color = fun),
            data = . %>% filter(!(fun %in% c("total")))) +
  geom_line(data = . %>% filter((fun %in% c("total"))), linewidth = 0.8) +
  geom_point(data = original_data) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = my_pallette_39) +
  labs(x = "Normalised Time (Proportion of Gait Cycle)",
       y = "Angle ($^{\\circ}$)") + # axis labels
  theme(legend.position = "none") # don't need a legend

# View plot
# Save Plots 
# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "basis-expansion-plot.tex"),
     width = 1 * doc_width_inches, 
     standAlone = TRUE,
     height = 0.5 *  doc_width_inches)
p3
dev.off()

tinytex::lualatex(file.path(plots_path, "basis-expansion-plot.tex"))



