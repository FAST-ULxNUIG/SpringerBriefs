# ======================================================================= #
# Preliminary plot of the reduction to discrete variables for the book #
# ======================================================================= #


# 1) Import Packages ------------------------------------------------------
library(fda)                                # CRAN v5.5.1 
library(tidyverse)                          # CRAN v1.3.1 
library(ggpubr)                             # CRAN v0.4.0 
library(Polychrome)                         # CRAN v1.5.1 
library(tikzDevice)                         # CRAN v0.12.3.1
library(data.table)                         # CRAN v1.14.2 

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-01", "figures")


# Plot Settings -----------------------------------------------------------
source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_update(strip.text = element_text(size = 10),
             axis.text = element_text(size = 9),
             axis.title = element_text(size = 10),
             plot.title = element_text(size = 11))


# Read in Data ------------------------------------------------------------
gait_lng <- readRDS(file = here::here("chapter-01/data/gait_lng.rds"))

# swtich to data_table as it's easier to group by
gait_lng <- as.data.table(gait_lng)

gait_lng_knee <- gait_lng[joint == "Knee Angle"]
gait_lng_knee

# get maximum over time points for each curve (basically filter each subject and retain only point where angle == max(angle)!)
gait_lng_max <- gait_lng_knee[, .(time, angle, is_max = angle == max(angle)), by = subject]
gait_lng_max <- gait_lng_max[is_max == TRUE]

p1 <- ggplot(data = gait_lng_knee) + # call plot, specify data
  aes(x = time, y = angle) + # map time to the x axis, angle on the y
  geom_line(mapping = aes(group = subject, color = subject), alpha = 1, lwd = 0.6) + # a different line and color for each subject
  scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = my_pallette_39) +
  labs(x = "Normalised Time (Proportion of Gait Cycle)",
       y = "Angle ($^{\\circ}$)",
       title = "\\textbf{(a)} Functional Data") + # axis labels
  theme(legend.position = "none") # don't need a legend
p1


p2 <- ggplot(data = gait_lng_knee) + # call plot, specify data
  aes(x = time, y = angle) + # map time to the x axis, angle on the y
  geom_line(mapping = aes(group = subject, color = subject), alpha = 0.075, lwd = 0.6) + # a different line and color for each subject
  geom_point(aes(color = subject), data = gait_lng_max, size=0.5) +
  geom_point(aes(x = mean(gait_lng_max$time), y = mean(gait_lng_max$angle)), size = 20, pch = 1) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = my_pallette_39) +
  labs(x = "Normalised Time (Proportion of Gait Cycle)",
       y = "Angle ($^{\\circ}$)",
       title = "\\textbf{(b)} Single Discrete Variable") + # axis labels
  theme(legend.position = "none") # don't need a legend
p2


# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "discrete-variable-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, standAlone = TRUE)
ggarrange(p1, p2, nrow=1)
dev.off()

tinytex::lualatex(file.path(plots_path, "discrete-variable-plot.tex"))


