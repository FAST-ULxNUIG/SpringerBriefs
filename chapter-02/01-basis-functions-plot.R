# ======================================================================= #
# Plot comparing the basis function expansions of Fourier and BSplines #
# Demo on the gait data set for smoothing #
# ======================================================================= #

# 1) Import Packages ------------------------------------------------------
library(fda) # CRAN v5.1.9
library(tidyverse) # CRAN v1.3.0
library(ggpubr) # CRAN v0.4.0
library(Polychrome) # CRAN v1.2.6
library(tikzDevice)

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-02", "figures")
source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12))


# -------------------------------------------------------------------------
gait_time <- as.numeric(dimnames(gait)[[1]])
bspline_9 <- create.bspline.basis(rangeval = range(gait_time), nbasis = 9, norder = 4)

# 9 Fourier basis functions:
fourier_9 <- create.fourier.basis(rangeval = range(gait_time), nbasis = 9)

# sequence of time points to evaluate data at 
evaluation_sequence <- seq(min(gait_time), max(gait_time), length.out = 501)

# evaluate the bases at these points
bspl_eval <- eval.basis(evalarg = evaluation_sequence, basisobj = bspline_9)
fourier_eval <- eval.basis(evalarg = evaluation_sequence, basisobj = fourier_9)

bs_eval_df<- bspl_eval %>% as.data.frame
f_eval_df<- fourier_eval %>% as.data.frame

colnames(bs_eval_df) <- colnames(f_eval_df) <- paste0("fun_", 1:9)



# calculate the sum of all weighted basis functions
# this is just the fitted function in each case (i.e. the weighted sum)
bs_eval_df$total <- rowSums(bs_eval_df)
f_eval_df$total <- rowSums(f_eval_df)

# data format: wide to long
bs_lng <- data.frame(time = evaluation_sequence, bs_eval_df) %>%
  gather(-time, key = "fun", value = "val_bs")

f_lng <- data.frame(time = evaluation_sequence, f_eval_df) %>%
  gather(-time, key = "fun", value = "val_f")

# and combine:
basis_df <- inner_join(bs_lng, f_lng, by = c("time", "fun")) %>%
  gather(- time, - fun, key = "basis_type", value = "val") %>%
  mutate(basis_type = case_when(basis_type == "val_bs" ~ "B-spline",
                                basis_type == "val_f" ~ "Fourier"))



# 7) Plot --------------------------------------------------------------------

# Get Palette for colors of basis functions
set.seed(1996) # randomised, so set seed:
my_pallette_39 <- createPalette(N = 39, seedcolors = c("#ff0000", "#00ff00", "#0000ff")) %>%
  as.vector

p3 <- basis_df %>%
  ggplot() +
  aes(x = time, y = val) +
  facet_wrap( ~ basis_type, scales ="free_y") +
  geom_line(mapping = aes(group = fun, color = fun),
            data = . %>% filter(!(fun %in% c("total")))) +
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

tikz(file.path(plots_path, "basis-functions-plot.tex"),
     width = 1.05 * doc_width_inches, 
     height = 0.5 *  doc_width_inches,
     standAlone = TRUE)
p3
dev.off()

tinytex::lualatex(file.path(plots_path, "basis-functions-plot.tex"))
