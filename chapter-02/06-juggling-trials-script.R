# ======================================================================= #
                    # Create a preliminary plot to show the juggling data #
# ======================================================================= #

# 1) Import Packages ------------------------------------------------------
library(fda)        # CRAN v5.1.9
library(tidyverse)  # CRAN v1.3.0
library(ggpubr)     # CRAN v0.4.0
library(registr)    # CRAN v1.0.0
library(Polychrome) # CRAN v1.2.6
library(tikzDevice) # CRAN v0.12.3.1
# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
data_path <- here::here("chapter-02", "juggling-data")
plots_path <- here::here("chapter-02", "figures")

source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_set(new = theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
             plot.subtitle = element_text(size = 12, hjust = 0.5))



# 3) Import Data ------------------------------------------------------------

# Juggling data should be downloaded and unzipped from the following link:
# https://archive.mbi.ohio-state.edu/vod/mbi-media/e-162/1400789954_juggling.zip
# the stored in a folder called Juggling within the working directory


# Import from .txt files
x <- read_table2(file.path(data_path,"Xcoord.txt"), col_names = FALSE) # x coordinates
y <- read_table2(file.path(data_path,"Ycoord.txt"), col_names = FALSE) # y coordinates
z <- read_table2(file.path(data_path,"Zcoord.txt"), col_names = FALSE) # z co-ordinates
ni <- as.vector(read_csv(file.path(data_path,"ni.txt"), col_names = FALSE))$X1 # no. of cycles in trial
seqn <- as.vector(t(read_table2(file.path(data_path,"seqn.txt"), col_names = FALSE)[1, ])) # no. of observations


# 4) Data Cleaning/ Wrangling ------------------------------------------------------
# currently zeros are filled in for NA (after trials end)
# fix this by finding last point in each trial and adding NA for entries after
x_t <- x
y_t <- y
z_t <- z
for(i in 1:length(ni)){
 if(ni[i] < nrow(x_t)){
   na_inds <- (ni[i]+1) : nrow(x_t)
   x_t[na_inds, i] <- y_t[na_inds, i] <- z_t[na_inds, i] <- NA
 }
  }

# add column names
colnames(x_t) <- colnames(y_t) <- colnames(z_t) <- paste0 ("trial", 1:ncol(x_t))
# time points to evaluate at
time_seq <- seq(0, nrow(x_t) * 5, length.out = nrow(x_t)) # converting frames to millisec (multiply by 5)

# Data: wide to long for each coord
x_lng <- data.frame(time_seq, x_t) %>%
  gather(-time_seq, key = "trial", value = "$x(t)$")

y_lng <- data.frame(time_seq, y_t) %>%
  gather(-time_seq, key = "trial", value = "$y(t)$")

z_lng <- data.frame(time_seq, z_t) %>%
  gather(-time_seq, key = "trial", value = "$z(t)$")

# Joint 3 long data sets to create another wide data containing xyz
xyz_wide <- inner_join(x_lng, y_lng, by = c("time_seq", "trial")) %>%
  inner_join(y = z_lng, by = c("time_seq", "trial"))
# wide to long again
xyz_lng <- xyz_wide %>%
  gather(-time_seq, - trial,  key = "coord", value = "fun_val")


# 5) Plot -----------------------------------------------------------------

set.seed(1996) # reproducible color palette with a different color for each trial
my_pallette_10 <- createPalette(N = 10, seedcolors = c("#ff0000", "#00ff00")) %>% as.vector()

# make plot
juggling_trials <- xyz_lng %>% 
  ggplot() +
  facet_wrap(~ coord, nrow = 3, scales = "free_y") +
  aes(x = time_seq, y = fun_val, group = trial, color = trial) +
  geom_line() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = my_pallette_10) +
  theme(legend.position = "none") +
  labs(x = "Time ($t$) in milliseconds", y = "Spatial Co-ordinate")

# preview
juggling_trials


tikz(file.path(plots_path, "juggling_trials.tex"),
     width = 7.87402, 
     height = 4.2949173228,
     standAlone = TRUE)
juggling_trials
dev.off()

tinytex::lualatex(file.path(plots_path, "juggling_trials.tex"))

# # save
# ggsave(plot = juggling_trials, filename = file.path(plots_path, "juggling_trials.png"),
#        device = "png", units = "cm", width = 20, height = (20/11) * 6)
