# ======================================================================= #
# Preliminary plot of the children's gait data for the book #
# ======================================================================= #

# Note: This code was originally used in EG's PhD thesis and was made publicly available:
# https://github.com/FAST-ULxNUIG/thesis-chapt-1-2/blob/main/code/analysis/chapter-01/gait-data-plot.R

# 1) Import Packages ------------------------------------------------------
library(fda)                                # CRAN v5.5.1 
library(tidyverse)                          # CRAN v1.3.1 
library(ggpubr)                             # CRAN v0.4.0 
library(Polychrome)                         # CRAN v1.5.1 
library(tikzDevice)                         # CRAN v0.12.3.1

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-01", "figures")

source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_update(strip.text = element_text(size = 11),
             axis.text = element_text(size = 10),
             axis.title = element_text(size = 11))


# 3) Data Wrangling ----------------------------------------------------------

# use the gait data set:
data('gait', package = "fda") # from fda package
dimnames(gait)
dim(gait)
# 20 x 29 x 2 array
# rows are time points
gait_time <- as.numeric(dimnames(gait)[[1]])
# subjects are columns
gait_subject_names <- dimnames(gait)[[2]]
# joints (hip, knee) are in the array 'slices'
hip_df <- data.frame(time = gait_time, gait[,, 1])
knee_df <- data.frame(time = gait_time, gait[,, 2])
head(hip_df)
head(knee_df)
# reshape the data to plot them:
hip_lng <- hip_df %>%
  gather(- time, key = "subject", value = "hip_angle")
head(hip_lng)

knee_lng <- knee_df %>%
  gather(- time, key = "subject", value = "knee_angle")
head(knee_lng)

# now the data are in 'long' form, 
# one column for time, one for angle, one for subject
# this means they are suitable for easy plotting with ggplot separately
# lets join them together to make one data frame
# we join them by the two common identifiers of the rows: dubject id and time points
# we use the inner_join function

hip_knee_df <- inner_join(x = hip_lng, y = knee_lng, by = c("time", "subject"))
head(hip_knee_df)

# we can reshape this data to make it 'longer' again
# make a column identifying the joint, and then another identifying the angle
# rather than separate angles for hip and knee

gait_lng <- hip_knee_df %>%
  gather(- time, - subject, key = "joint", value = "angle") %>%
  mutate(joint = case_when(joint == "hip_angle" ~ "Hip Angle",
                           joint == "knee_angle" ~ "Knee Angle") %>%
           as.factor)
head(gait_lng)

saveRDS(gait_lng, here::here("chapter-01/data/gait_lng.rds"))

# 4) Plot the Data -----------------------------------------------------------

# get colors
set.seed(1996)
my_pallette_39 <- createPalette(N = 39, seedcolors = c("#ff0000", "#00ff00", "#0000ff")) %>%
  as.vector

# Plot the data:
p1 <- ggplot(data = gait_lng) + # call plot, specify data
  aes(x = time, y = angle) + # map time to the x axis, angle on the y
  facet_wrap( ~ joint, scales = "free_y") + # create a sub-plot for the hip and knee angles
  geom_line(mapping = aes(group = subject, color = subject), alpha = 1, lwd = 0.6) + # a different line and colour for each subject
  scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = my_pallette_39) +
  labs(x = "Normalised Time (Proportion of Gait Cycle)",
       y = "Angle ($^{\\circ}$)") + # axis labels
  theme(legend.position = "none") # don't need a legend

p1

# Save Plots 
# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "gait-data-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, standAlone = TRUE)
p1
dev.off()

tinytex::lualatex(file.path(plots_path, "gait-data-plot.tex"))
