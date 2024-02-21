# ======================================================================= #
                    # Create annotated plot of choosing the smoothing parameter with gcv
# ======================================================================= #


# 1) Import Packages ------------------------------------------------------
library(fda)        # CRAN v5.1.9
library(tidyverse)  # CRAN v1.3.0
library(ggpubr)     # CRAN v0.4.0
library(ggforce)    # CRAN v0.3.2
library(ggtext)     # CRAN v0.1.2 
library(scales)     # CRAN v1.1.1
library(tikzDevice) # CRAN v0.12.3.1 


# 2) Settings -------------------------------------------------------------
# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-02", "figures")
theme_set(theme_bw())
source(file.path(functions_path, "theme_gunning.R"))
theme_update(axis.text = element_text(size = 11),
             axis.title = element_text(size = 11),
             panel.grid.minor = element_blank(),
             strip.text = element_text(size = 11),
             plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
             plot.subtitle = element_text(size = 11, hjust = 0.5),
             axis.ticks = element_blank(),
             strip.background = element_rect(size = 0.5))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937



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


# 4) Smoothing ------------------------------------------------------------
# Do a gcv search.
N_obs <- length(gait_trial)
lambdas <- 10 ^ seq(-10, -2, by = 1)
N_eval <- 501
evaluation_sequence <- seq(min(gait_time), max(gait_time), length.out = N_eval)
rich_basis <- create.bspline.basis(rangeval = range(gait_time), nbasis = 20)
lambda_store <- matrix(data = NA, nrow = N_eval, ncol = length(lambdas))

# set up names of clumn 
colnames(lambda_store) <- paste0("$\\lambda = 10^{", log10(lambdas), "}$")


gcv_store <- numeric(length = length(lambdas))
for(i in 1:length(lambdas)){
  smooth_lam <- smooth.basis(argvals = gait_time, y = gait_trial,
                             fdParobj = fdPar(fdobj = rich_basis,
                                              Lfdobj = 2, lambda = lambdas[i]))
  lambda_store[, i] <- eval.fd(fdobj = smooth_lam$fd, evalarg = evaluation_sequence)
  gcv_store[i] <- smooth_lam$gcv
}


# wrangle the data
lambda_df <- lambda_store %>%
  as.data.frame 
lambda_df$t <- evaluation_sequence
original_data <- data.frame(t = gait_time, x_true = gait_trial)


# 5) Create the 2 Plots ------------------------------------------------------


## Smoothes curves with varying lambda -------------------------------------
g2 <- lambda_df %>%
  gather(-t, key = "lam_val", value = "xhat") %>%
  mutate(lam_val = factor(lam_val, levels = paste0("$\\lambda = 10^{", log10(lambdas), "}$"))) %>%
  ggplot(aes(x = t, y = xhat, group = lam_val, color = lam_val)) +
  geom_point(data = original_data, inherit.aes = F, aes(x = t, y = x_true),
             size = 1, alpha = 0.6) + 
  facet_wrap( ~ lam_val) + #, labeller = "label_parsed") +
  geom_line(lwd = 0.75) +
  labs(y = "Angle ($^{\\circ}$)",
       x = "Normalised Time (Proportion of Gait Cycle)") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1),labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(legend.position = "none",
        legend.title = element_blank())

g2


rich_lambdas <- 10 ^ seq(-10, -2, by = 0.1)
rich_gcv_store <- numeric(length = length(rich_lambdas))
for(i in 1:length(rich_lambdas)){
  smooth_lam <- smooth.basis(argvals = gait_time, y = gait_trial,
                             fdParobj = fdPar(fdobj = rich_basis,
                                              Lfdobj = 2, lambda = rich_lambdas[i]))
  rich_gcv_store[i] <- smooth_lam$gcv
}


rich_lam_df <- data.frame(lam = rich_lambdas, gcv = rich_gcv_store)
lam_df <- data.frame(lam = lambdas, gcv = gcv_store, lam_val = factor(paste0("λ==10^", log10(lambdas)), levels = paste0("λ==10^", log10(lambdas))))

acceptable_range <- rich_lam_df %>%
  filter(gcv <= 2.7)
acceptable_range$label_rect <- "grey"

textbox_df <- data.frame(lab = "<span style='font-size:14pt; color:black'> The results of the grid search <br> suggest that a suitable value <br> of $\\lambda$ should lie in this region </span>", y = 15, x = 10^-8, color = "black")
textbox_df2 <- data.frame(lab = "<span style='font-size:14pt; color:black'> After examining the fits in the <br>middle row we might opt for<br> the smoothest fit at **<span style='color:#00B9E3'>$\\lambda = 10$<sup>-5</sup></span>**  <br>  or refine our search around <br> this area </span>", y = 28, x = 10^-6, color = c("black"))


## GCV plot --------------------------------------------------------------

g1 <- ggplot(data = rich_lam_df) +
  aes(x = lam, y = gcv) +
  geom_line() +
  geom_mark_rect(data = acceptable_range,
                 aes(fill = label_rect),
                 expand = c(0.01), fill = "darkgrey",
                 color = "grey", linetype = 2, lwd = 1) +
  geom_mark_rect(aes(fill = lam_val, filter = lam_val == "λ==10^-6"), data = lam_df, col = "#00BFC4",  fill = "#00BFC4", expand = 0.03) +
  geom_point(data = lam_df, aes(color = lam_val), size = 2)+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-2, 40), expand = c(0, 0)) +
  labs(title = "Generalised Cross Validation",
       y = "GCV$(\\lambda)$", x = "Smoothing Parameter $\\lambda$") +
  scale_x_log10(expand = c(0.02, 0),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_mark_rect(data = acceptable_range,
                 aes(fill = label_rect), fill = "darkgrey", expand = c(0.025),
                 color = "grey", linetype = 2, lwd = 1) 
#+
  # geom_richtext(inherit.aes = F, aes(x =x, y = y, label = lab),
  #               data = textbox_df, fill = "lightgrey", alpha = 0.3, col = "darkgrey",
  #               label.padding = unit(0.7, "line")) +
  # annotate("segment", x = 10^-6.5, xend = 10^-6.5, y = 11.5, yend = 3.5, col = "darkgrey", lwd = 0.75)
  # geom_richtext(inherit.aes = F, aes(x =x, y = y, label = lab),
  #               data = textbox_df2, fill = "white", alpha = 1, col = "#00B9E3", size = 4,
  #               label.padding = unit(0.65, "line")) +
  # annotate("segment", x = 10^-5, xend = 10^-5, y = 23, yend = 3.5, col = "#00B9E3", lwd = 0.75, alpha = 0.75) +
  # #annotation_logticks(sides = "b") +


g1



## Combine ---------------------------------------------------------------
(comb_gcv <- ggarrange(g2, g1, align = "hv"))



tikz(file.path(plots_path, "gcv-plot.tex"),
     width = 1.5 * doc_width_inches, 
     height = 0.7 *  doc_width_inches,
     standAlone = TRUE)
comb_gcv
dev.off()

tinytex::lualatex(file.path(plots_path, "gcv-plot.tex"))
