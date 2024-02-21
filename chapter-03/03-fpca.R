# ======================================================================== #
# Plot for visualizations provided by FPCA.

# Plot 1: Reconstruction with the FPCs as basis functions

# Plot 2: Showing the FPC functions, them added and subtracted to the mean
# and the highest and lowest scoring curves.

# Need to already have run GaitRecMeanVar_Script.R

# ======================================================================== #

# 1) Load Packages -----------------------------------------------------------
library(readr)        # CRAN v1.3.1
library(tidyverse)    # CRAN v1.3.0
library(fda)          # CRAN v5.1.9
library(RColorBrewer) # CRAN v1.1-2
library(ggtext)       # [github::wilkelab/ggtext] v0.1.0
library(ggpubr)       # CRAN v0.4.0
library(Polychrome)   # CRAN v1.2.6
library(tikzDevice)   # CRAN v0.12.3.1

# 2)  Plot settings -----------------------------------------------------------
theme_set(new = theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5),
             plot.subtitle = element_text(size = 14, hjust = 0.5))

source(here::here("functions", "theme_gunning.R"))
theme_gunning()

# 3) Read in data ----------------------------------------------------------------------

data.path <- here::here("chapter-03", "data")

# Datasets, already converted to functions:
average_profiles_reg <- readRDS(
  file = file.path(data.path, "GaitRecHealthyReg.rds"))

selected_sessions <- readRDS(
  file = file.path(data.path, "SelectedSessions.rds"))



# 2) Perform FPCA -------------------------------------------------------------

average_vGRF_pca <- pca.fd(fdobj = average_profiles_reg$regfd,
                           nharm = 3, 
                           centerfns = T)

# Plot FPCA
par(mfrow = c(1, 1))
plot.pca.fd(average_vGRF_pca)
# Nice, simple patterns found by the FPCs


# 3) Plotting FPCA Results ----------------------------------------------------------------


## i) Plot 1 ------------------------------------------------------------------

# Plot 1 will be different reconstructions using
# 1) Just the Mean
# 2) mean + FPC1
# 3) mean + FPC1 + FPC2 
# 4) mean + FPC1 + FPC2 + FPC3

eval_points <- 0:100

# evaluate the mean function:
mean_fun <- mean.fd(x = average_profiles_reg$regfd) %>% eval.fd(evalarg = eval_points)

# evaluate the fpc functions 'harmonics':
fpc_funs <- eval.fd(evalarg = eval_points, fdobj = average_vGRF_pca$harmonics)

# get the fpca scores
fpc_scores <- average_vGRF_pca$scores

dim(fpc_scores)
dim(fpc_funs)

# we are going to create three data frames:
# 1) the mean curve (same for each subject) [represent_1]
# 2) mean + fpc score 1 * fpc1 for each subject [represent_2]
# 3) mean + fpc score 1 * fpc1 + fpc score 2 * fpc2 for each subject [represent_3]


### Represent_1 -------------------------------------------------------------

n_subject <- length(selected_sessions)
# rep the mean function n_subj times in a matrix
represent_1 <- matrix(mean_fun, nrow = length(eval_points), ncol = n_subject)
# join it to a data frame and add time points
represent_1_wide <- data.frame(time_seq = eval_points, represent_1)
head(represent_1_wide) # peak at it
colnames(represent_1_wide)[-1] <- selected_sessions # name cols as session IDs

# reshape the data to long form for plotting
represent_1_long <- represent_1_wide %>% 
  gather(- time_seq, key = "session", value = "mean")

head(represent_1_long) # peak


### Represent 2 -------------------------------------------------------------

scores_fpc1 <- fpc_scores[, 1] 
eval_fpc1 <- fpc_funs[, 1]
represent_2 <- represent_1 + matrix(eval_fpc1, nrow = length(eval_points), ncol = n_subject) %*% diag(scores_fpc1)
# same procedure for reshaping
# join it to a data frame and add time points
represent_2_wide <- data.frame(time_seq = eval_points, represent_2)
head(represent_2_wide) # peak at it
colnames(represent_2_wide)[-1] <- selected_sessions # name cols as session IDs

# reshape the data to long form for plotting
represent_2_long <- represent_2_wide %>% 
  gather(- time_seq, key = "session", value = "mean + fpc1")

head(represent_2_long) # peak


### Represent 3 -------------------------------------------------------------

scores_fpc2 <- fpc_scores[, 2] 
eval_fpc2 <- fpc_funs[, 2]
represent_3 <- represent_2 + matrix(eval_fpc2, nrow = length(eval_points), ncol = n_subject) %*% diag(scores_fpc2)
# same procedure for reshaping
# join it to a data frame and add time points
represent_3_wide <- data.frame(time_seq = eval_points, represent_3)
head(represent_3_wide) # peak at it
colnames(represent_3_wide)[-1] <- selected_sessions # name cols as session IDs

# reshape the data to long form for plotting
represent_3_long <- represent_3_wide %>% 
  gather(- time_seq, key = "session", value = "mean + fpc1 + fpc2")

head(represent_3_long) # peak


### Represent 4 -------------------------------------------------------------

scores_fpc3 <- fpc_scores[, 3] 
eval_fpc3 <- fpc_funs[, 3]
represent_4 <- represent_3 + matrix(eval_fpc3, nrow = length(eval_points), ncol = n_subject) %*% diag(scores_fpc3)
# same procedure for reshaping
# join it to a data frame and add time points
represent_4_wide <- data.frame(time_seq = eval_points, represent_4)
head(represent_4_wide) # peak at it
colnames(represent_4_wide)[-1] <- selected_sessions # name cols as session IDs

# reshape the data to long form for plotting
represent_4_long <- represent_4_wide %>% 
  gather(- time_seq, key = "session", value = "mean + fpc1 + fpc2 + fpc3")

head(represent_4_long) # peak

par(mfrow = c(1, 1))
plot(average_vGRF_pca$values, type= "b")
plot(100 * cumsum(average_vGRF_pca$values) / sum(average_vGRF_pca$values), type = "b", ylab = "% Variance Explained")

### Join up representations -------------------------------------------------
# joining them all up makes the data 'wide' agin
combined_representations_wide <- inner_join(represent_1_long, represent_2_long,by = c("time_seq", "session")) %>%
  inner_join(represent_3_long, by = c("time_seq", "session"))  %>%
  inner_join(represent_4_long, by = c("time_seq", "session"))

combined_representations_long <- combined_representations_wide %>%
  gather(-time_seq, - session, key = "representation", val = "fun_val")

head(combined_representations_long)


### Plot --------------------------------------------------------------------

set.seed(1996)
# randomly select colors:
my_pallette_nsub <- createPalette(N = n_subject,
                                seedcolors = c("#ff0000", "#00ff00")) %>% as.vector()

# # we want formulae in strip text
# lookup_labels <- c(expression(italic(bar(x)(t))),
#                    expression(bar(italic(x))~italic((t)) + italic(f[i1]) ~ italic(xi[1](t))),
#                    expression(bar(italic(x))~italic((t)) + italic(f[i1]) ~ italic(xi[1](t)) + italic(f[i2]) ~ italic(xi[2](t))),
#                    expression(bar(italic(x))~italic((t)) + italic(f[i1]) ~ italic(xi[1](t)) + italic(f[i2]) ~ italic(xi[2](t)) + italic(f[i3]) ~ italic(xi[3](t))))

lookup_labels <- c("$\\bar{x}(t)$",
                   "$\\bar{x}(t) + f_{i1} \\xi_1(t)$",
                   "$\\bar{x}(t) + f_{i1} \\xi_1(t) + f_{i2} \\xi_2(t)$",
                   "$\\bar{x}(t) + f_{i1} \\xi_1(t) + f_{i2} \\xi_2(t) + f_{i3} \\xi_3(t)$")


                      # create Data Set to Add Variance Labels:
variance_labels <- round(cumsum(average_vGRF_pca$varprop), 2)
var_prop_rep <- factor(x = c("mean + fpc1",
                             "mean + fpc1 + fpc2",
                             "mean + fpc1 + fpc2 + fpc3"),
                       levels = c("mean + fpc1",
                                  "mean + fpc1 + fpc2",
                                  "mean + fpc1 + fpc2 + fpc3"),
                       labels = lookup_labels[-1])
var_df <- data.frame(label = paste0("Variance Explained = $", 100 * variance_labels, "\\%$"), representation = var_prop_rep)





combined_representations_long %>%
  mutate(representation = factor(representation, levels = c('mean', "mean + fpc1",
                                                            "mean + fpc1 + fpc2",
                                                            "mean + fpc1 + fpc2 + fpc3"),
                                    labels = lookup_labels)) %>%
  
  ggplot() +
  facet_wrap( ~ representation) +
  aes(x = time_seq, y = fun_val) +
  geom_line(aes(group = session, color = session)) +
  scale_color_manual(values = my_pallette_nsub) +
  #geom_line(data = mean_df, col = "black", lwd = 0.75) +
  theme(legend.position  = "none",
        strip.text = element_text(family = "serif", size = 14)) +
  geom_label(mapping = aes(x = 50, y = 0.3, label = label), fill = alpha("lightblue", 0.2), size = (5/14) * 11,
             data = var_df, label.padding = unit(0.4, "lines"), family = "sans") +
  theme(strip.background = element_rect(fill = alpha("lightblue", 0.22)),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 12)) +
  labs(x = "Normalised Time ($\\%$ of Stance)",
       y = "vGRF (Normalised to BW)") -> fpc_plot1

fpc_plot1

plots_path <- here::here("chapter-03", "figures")
tikz(file = file.path(plots_path, "KL_rep.tex"), width = 18 * 0.393701, height = 17 * 0.393701, standAlone = TRUE)
fpc_plot1
dev.off()
tinytex::lualatex(file.path(plots_path, "KL_rep.tex"))


## i) Plot 2 ------------------------------------------------------------------

# here we are going to plot three plots together:
# 1) Just FPC1
# 2) Mean +- C * FPC1 'legacy' fda plot
# 3) The highest and lowest scorers on FPC1.


### Data wrangling ----------------------------------------------------------
high_low_scorer_df <- data.frame(time_seq = eval_points, 
           low_scorer = eval.fd(evalarg = eval_points, average_profiles_reg$regfd[which.min(scores_fpc1)]),
           high_scorer = eval.fd(evalarg = eval_points, average_profiles_reg$regfd[which.max(scores_fpc1)])) %>%
  gather(-time_seq, key = "fun", val = "fun_val")
high_low_scorer_df$plot_type <- "high_low"

plot_labels <- c("FPC1", "Mean ± C × FPC1", "**<span style='color:#00B0F6'>Highest</span>** and **<span style='color:#F8766D'>Lowest</span>** FPC1 Scores")


### Plotting ----------------------------------------------------------------
data.frame(time_seq = eval_points, mean_fun = mean_fun, fpc1 = eval_fpc1) %>%
  mutate(mean_plus = mean + 2 * sd(scores_fpc1) * fpc1,
         mean_minus = mean - 2 * sd(scores_fpc1) * fpc1) %>%
  gather(-time_seq, key = "fun", val = "fun_val") %>% 
  mutate(plot_type = ifelse(fun == "fpc1", "pc_fun", "add_sub")) %>%
  rbind(high_low_scorer_df) %>%
  mutate(plot_type = factor(plot_type, levels = c("pc_fun", "add_sub", "high_low"),
                            labels = plot_labels)) %>%
  ggplot() +
  facet_wrap( ~ plot_type) +
  aes(x = time_seq, y = fun_val, group = fun, color = fun) +
  geom_line(data = . %>% filter(! fun %in% c("mean_minus", "mean_plus")), lwd = 0.8) +
  geom_point(data = . %>% filter(fun %in% c("mean_minus", "mean_plus")) %>% filter(time_seq %in% seq(0, 101, by =2)),
             aes(shape = fun, size = fun)) +
  scale_size_manual(values = c(2.5, 2.5)) +
  scale_shape_manual(values = c("-", "+")) +
  scale_color_manual(values = c("black", "black", "black", "black", "#F8766D", "#00B0F6")) +
  theme(tex = element_text(family = "serif")) +
  theme(legend.position = "none",
        strip.text = element_markdown(size = 10),
        strip.background = element_rect(fill = alpha("grey", 0.2)),
        axis.title = element_text(size = 11)) +
  labs(x = "Normalised Time (% of Stance)",
       y = "vGRF (Normalised to BW)") +
  geom_hline(yintercept = 0, lty = 3) -> fpc_plot2

fpc_plot2

# save
ggsave(plot = fpc_plot2,
       filename = here::here("chapter-03",
                             "figures",
                             "fpc1_plot.pdf"),
       device = "pdf",
       units = "cm",
       width = 19,
       height = 7)




