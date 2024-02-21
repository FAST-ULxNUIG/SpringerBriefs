# ======================================================================= #
# Plot the mean and covariance of the vgrf data for a group of
# healthy subjects in the gait data
# ======================================================================= #


# 1) Load Packages -----------------------------------------------------------
library(readr) # CRAN v1.3.1
library(tidyverse) # CRAN v1.3.0
library(fda) # CRAN v5.1.9
library(RColorBrewer) # CRAN v1.1-2
library(ggtext) # [github::wilkelab/ggtext] v0.1.0
library(ggpubr) # CRAN v0.4.0
library(tikzDevice)
library(Polychrome) # CRAN v1.2.6
# conflicted::conflict_prefer(name = "select", winner = "dplyr")
# conflicted::conflict_prefer(name = "filter", winner = "dplyr")

# 2) Plot settings --------------------------------------------------------

source(here::here("functions", "theme_gunning.R"))
theme_gunning()
theme_update(axis.text = element_text(size = 12),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5),
             plot.subtitle = element_text(size = 14, hjust = 0.5))

# Note: Skip to line 240 and read in the RDS files (commented out)
# if you just want to make the plots.

# 3) Read in Data ------------------------------------------------------------
# Exploratory analysis of the GaitRec data:
# https://www.nature.com/articles/s41597-020-0481-z

# Download:
# Force = vertical
# Leg: Left
# Processing: Smoothed and Normalized
# This is a big data set, it might take a while:
GRF_F_V_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063191")
# And associated MetaDeta, e.g., session and subject information:
GRF_metadata <- read_csv("https://ndownloader.figshare.com/files/22062960")
#.. you might have to wait.. it Will download!



# 4) Exploratory data analysis -----------------------------------------------
# Aim: become familiar with the structure of the meta data
# to pick out a nice subset to demonstrate exploratory techniques.

head(GRF_metadata)
# We can get all info from here: https://www.nature.com/articles/s41597-020-0481-z/tables/4
#Identifiers
#SUBJECT_ID 	integer 	— 	Unique identifier of a subject
#SESSION_ID 	integer 	— 	Unique identifier of a session
#Labels
#CLASS_LABEL 	string 	— 	Annotated class labels
#CLASS_LABEL_DETAILED 	string 	— 	Annotated class labels for subclasses
#Subject Metadata
#SEX 	binary 	— 	female=0, male=1
#AGE 	integer 	years 	Age at recording date
#HEIGHT 	integer 	centimeter 	Body height in centimeters
#BODY_WEIGHT 	double 	𝑘𝑔𝑚𝑠2
#Body weight in Newton
#BODY_MASS 	double 	kg 	Body mass
#SHOE_SIZE 	double 	EU 	Shoe size in the Continental European System
#AFFECTED_SIDE 	integer — left=0, right=1, both=2
#Trial Metadata
#SHOD_CONDITION 	integer 	— 	barefoot & socks=0, normal shoe=1, orthopedic shoe=2
#ORTHOPEDIC_INSOLE 	binary 	— 	without insole=0, with insole=1
#SPEED 	integer 	— 	slow=1, self-selected=2, fast=3 walking speed
#READMISSION 	integer 	— 	indicates the number of re-admission=0 … n
#SESSION_TYPE 	integer 	— 	initial measurement=1, control measurement=2, initial measurement after readmission=3
#SESSION_DATE 	string 	— 	date of recording session in the format “DD-MM-YYYY”
#Train-Test Split Information
#TRAIN 	binary 	— 	is part (=1) or is not part (=0) of TRAIN
#TRAIN_BALANCED 	binary 	— 	is part (=1) or is not part (=0) of TRAIN_BALANCED
# TEST 	binary 	— 	is part (=1) or is not part (=0) of TEST

# Skim the meta data:
skimr::skim(GRF_metadata)

# How many different subjects?
GRF_metadata  %>%
  summarise(n_distinct(SUBJECT_ID))
#  2295 (!!!! big)

# How many different subjects in the healthy class (HC)?
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  summarise(n_distinct(SUBJECT_ID))
# 211

# How many sessions do each of them have?
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>%
  summarise(n_sesh = n_distinct(SESSION_ID)) %>%
  summarise(min(n_sesh), max(n_sesh), median(n_sesh))

# can look on a bar chart too..
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>%
  summarise(n_sesh = n_distinct(SESSION_ID)) %>%
  ggplot(aes(x = n_sesh)) +
  geom_bar()

# Are these all on the same day and under separate conditions?
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>%
  summarise(n_sesh = n_distinct(SESSION_DATE)) %>%
  summarise(min(n_sesh), max(n_sesh), median(n_sesh))
# Answer, a clear yes:.


# do the SESSION TYPE's differ?
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  summarise( min(SESSION_TYPE), max(SESSION_TYPE))
# No, they are all coded 1



# 5) Select a subset of the sessions for Exploratory Analysis  ---------------------

# Here is the subset we will take
# try and make it  as standard as possible
# eg. healthy controls, self-selected speed, wearing normal shoes.....
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>% # we are going to select one session per subject
  # filter(SESSION_TYPE == 1) %>% 
  filter(SPEED == 2) %>% # we'll go with self selected speed
  # filter(SESSION_DATE == min(SESSION_DATE)) %>% # make sure it's their first session 
  filter(SHOD_CONDITION == 1) %>% # all shod in normal shoe
  ungroup() %>% 
  summarise(n_distinct(SESSION_ID)) # check how man subjects this is
# 208
# Nice sample size for demonstration.

# Make sure it is one session per subject
# i.e we get all 1's here.
GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>% 
  # filter(SESSION_TYPE == 1) %>%
  filter(SPEED == 2) %>%
  #filter(SESSION_DATE == min(SESSION_DATE)) %>% 
  filter(SHOD_CONDITION == 1) %>%
  count() %>%
  pull
# all good!!!

selected_sessions <- GRF_metadata %>% filter(CLASS_LABEL == "HC") %>%
  group_by(SUBJECT_ID) %>% 
  # filter(SESSION_TYPE == 1) %>%
  filter(SPEED == 2) %>%
  #filter(SESSION_DATE == min(SESSION_DATE)) %>% 
  filter(SHOD_CONDITION == 1) %>%
  pull(SESSION_ID) # Choose the session ID variable.

# make sure every element is a unique session ID
stopifnot(length(unique(selected_sessions)) == length(selected_sessions)) 

# save selected sessions as an RDS obj for later.
saveRDS(object = selected_sessions, file = here::here("chapter-03", "data", "SelectedSessions.rds"))


# 6) Subset the GRF data -----------------------------------------------------
# This is processed data (normalised to [0,1])
head(GRF_F_V_PRO_left)

# We will subset it to look at only the session we have selected above
# Let's see how many trials per session/subject: 
GRF_F_V_PRO_left %>%
  filter(SESSION_ID %in% selected_sessions) %>%
  group_by(SESSION_ID) %>%
  summarise(num_trials = n_distinct(TRIAL_ID)) %>%
  summarise(min(num_trials), max(num_trials), median(num_trials))

# min max median
#  4  16 10

# the data are already smoothed ('filtered')
# think 35 basis functions and a tiny smoothing parameter will do.
# we will use cubic spline, penalize 2nd derivative with a tiny penalty
bsplines_35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 35, norder = 4)


# 7) Register and Average and represent the data as an FD object ----------

# we want an iid sample for mean, covar and pca.
# we will average all trials within a session to have 1 per subj
# therefore, a light registration within each subject will help this average

##### N.B. ##########
# ^ the above is probably not necessary. and probably shouldn't be taken to be done in practice.

# Set up basis and smoothing for warping functions:
wbasis <- create.bspline.basis(rangeval = c(0, 100), nbasis = 5, norder = 4)
wcoef <- matrix(0, 5, 1)
wb_obj <- fd(basisobj =  wbasis, coef = wcoef)
w_par <- fdPar(fdobj = wb_obj, Lfdobj = 2,lambda = 10^-4) 
# I've chosen quite a heavy penalty on the warping functions so to only do light warping


# Going to loop through every subject, and get an average, with some light warping.
# store the coefficients of the average for each subject in this matrix
coefficients_average <- matrix(data = NA, nrow = bsplines_35$nbasis, ncol = length(selected_sessions))

# Start loop
for(sesh in 1:length(selected_sessions)){ 
  
  sesh_id <- selected_sessions[sesh] # loop through each session
  GRF_data <- GRF_F_V_PRO_left %>% filter(SESSION_ID == sesh_id) # subset GRF data to contain session
  n_trials <- nrow(GRF_data) # number of trials
  trial_id <- GRF_data$TRIAL_ID # ID of the trials
  grf_fd <- smooth.basis(argvals = 0:100, # smooth the data with the 35 basis functions
                         y = t(GRF_data[, -c(1:3)]), # and a tiny roughness penalty (10^-10)
                         fdParobj = fdPar(fdobj = bsplines_35, Lfdobj = 2, lambda = 10^-10))
  # register the trials before computing an average, using the continous reg (Ramsay and Silverman (2005)
  registered_funs <- register.fd(yfd = grf_fd$fd, WfdParobj = w_par, crit = 2)
  
  # Make a plot to see if registration performed ok:
  #jpeg(filename = paste0("Ouputs/Diagnostic_Plots/", sesh, ".jpg"), width = 900, height = 1000)
  par(mfrow = c(3, 1), cex = 1) # set up for side-by-side plot
  plot(grf_fd) # plot unregistered data
  lines(mean.fd(grf_fd$fd), type = "l", lwd = 2) # overlay thick black mean
  title("Unregistered")
  plot(registered_funs$regfd) # plot registered data 
  lines(mean.fd(registered_funs$regfd), lwd = 2, col = "black") # overlay thick black mean
  title("Registered")
  plot(mean.fd(registered_funs$regfd)) # compare registered and unregistered average trial
  lines(mean.fd(grf_fd$fd), col = "red")
  legend("bottom", c("registered", "unregistered"), col = c("black", "red"), lty = c(1, 1))
  title("Comparison")
  #dev.off()
  
  # Store the coefficients of the average of the trials for this session.
  coefficients_average[, sesh] <- mean.fd(registered_funs$regfd)$coefs
}

# create an fd object from the averaged profiles:
average_profiles <- fd(coef = coefficients_average, basisobj = bsplines_35)
average_profiles$fdnames$reps <- selected_sessions # add session names
average_profiles$fdnames$funs <- "grf"

# register these profiles, with an even smoother warps
w_par <- fdPar(fdobj = wb_obj, Lfdobj = 2,lambda = 10^-2) # lambda = 10^-2 for smoother warping functions 
# registration (again)
average_profiles_reg <- register.fd(yfd = average_profiles, WfdParobj = w_par, crit = 2)
# save this for later work before plotting.
# RDS is good for reproducibility.
saveRDS(object = average_profiles_reg,
        file = here::here("chapter-03", "data", "GaitRecHealthyReg.rds"))

# ***** Can skip to here to recreate plots ******
# Uncomment and run:
#average_profiles_reg <- readRDS(file = "Outputs/GaitRecHealthyReg.rds")
#selected_sessions <- readRDS(file = "Outputs/SelectedSessions.rds")

# Let's have a look at the data pre and post:
par(ask= F, mfrow = c(2, 1))
plot(average_profiles)
lines(mean.fd(average_profiles), lwd = 4)
plot(average_profiles_reg$regfd)
lines(mean.fd(average_profiles_reg$regfd), lwd = 4)
plot(mean.fd(average_profiles_reg$regfd), lwd = 1)
lines(mean.fd(average_profiles), lwd = 1, col = "red")
plot(average_profiles_reg$Wfd)

# Save pre and post comparisons that can be manually used to check registration.
#for(i in 1:length(selected_sessions)){
#  jpeg(paste0("Outputs/average", selected_sessions[i], ".jpg"), width = 600, height = 500)
#  lines(average_profiles_reg$regfd[i], lwd = 2, col = "red")
#  plot(average_profiles[i], lwd = 2)
#  title(paste0("Session ID", selected_sessions[i]))
#  dev.off()
#}



# 7) Plots of the mean and Covariance Structure ---------------------------


# 7.1 Covariance ----------------------------------------------------------

co_var_reg <- eval.bifd(var.fd(fdobj1 = average_profiles_reg$regfd), sevalarg = 0:100, tevalarg = 0:100 )
co_var_unreg <- eval.bifd(var.fd(fdobj1 = average_profiles), sevalarg = 0:100, tevalarg = 0:100)

# plot covariance surfaces
filled.contour(x = 0:100, y = 0:100, z = co_var_unreg)
filled.contour(x = 0:100, y = 0:100, z = co_var_reg)

co_var <- co_var_reg
# we can see registration has 'focused' the variance more
# useful for didactic purposes in chapter (practically more useful?)
# we'll use the registered data for our plots

co_var <- as.data.frame(co_var)
colnames(co_var) <- 0:100
co_var$s <- 0:100


# set up function to inerpolate a color palette:
getPalette <- colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))

# Plot of Covariance Function:
# function to label breaks for the contour plot.
breaks_labeller <- function(string_break){
  cleaned_string <- str_remove(string = string_break, pattern =  "]") %>%
    str_remove(pattern = "\\(") %>%
    str_remove(pattern = fixed(" "))
  
  inf_string <- str_extract(string = cleaned_string,
                            pattern = "[^,]*") %>%
    as.numeric %>%
    round(digits = 2)
  
  sup_string <- str_extract(string = cleaned_string,
                            pattern =  "[^,]*$") %>%
    as.numeric %>%
    round(digits = 2)
  
  rounded_lab = paste0("(", inf_string, ", ", sup_string, "]")
  return(rounded_lab)
}
# test
breaks_labeller(c("(-0.4263, -0.3654]", "(-0.4255, -0.4554]"))
# plot contour
(co_var_plot <- co_var %>%
  gather(-s, key = "t", value = "var") %>%
  mutate(t = as.numeric(t)) %>%
  ggplot() + # note rescale x 100 so axis labels have 2 rather than 4 decimal points
  geom_contour_filled(mapping = aes(x = s, y = t, z = 100*var), bins = 15, alpha = 0.9) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = getPalette(n = 15), labels = breaks_labeller) +
  theme(axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.box.margin = margin(l=10, r=10),
        legend.background = element_rect(color = "black")) +
  labs(title = "Sample Covariance Function $\\hat{v}(s, t)$",
       x = "$s$", 
       y = "$t$",
       fill = "$\\hat{v}(s, t)$ ($\\times 10^2$)"))




# 7.2 Mean Function -------------------------------------------------------

# Evaluated the registered data:
# We will plot in the background
averaged_reg_eval <- eval.fd(evalarg = 0:100, fdobj = average_profiles_reg$regfd)
# shape it as a data frame to plot:
rownames(averaged_reg_eval) <- 0:100
colnames(averaged_reg_eval) <- selected_sessions
# we are going to plot the meamn, store this in df:
averaged_reg_mean <- apply(X = averaged_reg_eval, MARGIN = 1, mean)
av_mean_df <- data.frame("t_val" = 0:100, "mean_val" = averaged_reg_mean)

# make plot:
mean_plot <- averaged_reg_eval %>%
  as.data.frame %>%
  rownames_to_column(var = "t_val") %>%
  mutate(t_val = as.numeric(as.character(t_val))) %>%
  gather(-t_val, key = "subject", value = "grf_val") %>%
  ggplot() +
  theme(plot.title = element_text(face = "bold"),
        legend.box.margin = margin(l=10, r=10),
        legend.background = element_rect(color = "black")) + 
  scale_y_continuous(expand = c(0, 0)) +
  aes(x = t_val, y = grf_val) +
  geom_line(mapping = aes(group = subject), color = "darkgrey", alpha = 0.5) +
  geom_line(data = av_mean_df, aes(y = mean_val), lwd = 1.5) +
  labs(y = "vGRF (Normalised to BW)",
       x = "Normalised Time ($\\%$ of Stance)",
       title = "Sample Mean Function $\\bar{x}(t)$")

mean_plot


# 7.3 Arrange and Save Plots ----------------------------------------------

# Arrange plots side by side:
mean_covar_plot <- ggarrange(mean_plot, co_var_plot, 
                             ncol = 2, 
                             widths = c(0.5, 0.5), align = "h",
                             common.legend = TRUE,
                             legend = "right")
mean_covar_plot


plots_path <- here::here("chapter-03", "figures")
tikz(file = file.path(plots_path, "MeanCoVar.tex"), width = 10.6299 * 0.85, height = 4.33071 * 0.85, standAlone = TRUE)
mean_covar_plot
dev.off()
tinytex::lualatex(file.path(plots_path, "MeanCoVar.tex"))

