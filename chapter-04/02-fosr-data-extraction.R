# ======================================================================= #
# Extract data set for the GaitRec data set to fit a function-on-scalar 
# linear model. This script just extracts the data set that we will fit 
# the model in the next script FLM_Model_Fitting_Script.R
# ======================================================================= #

# Load Packages -----------------------------------------------------------
library(readr) # CRAN v1.3.1
library(fda) # CRAN v5.1.9
library(tidyverse) # CRAN v1.3.0


# 1) Read in Data ------------------------------------------------------------
# Use the GaitRec data:
# https://www.nature.com/articles/s41597-020-0481-z

# Download:
# Force = vertical (F & V)
# Leg: Left (left)
# Processing: Smoothed and Normalized (PRO)
# This is a big data set, it might take a while:
GRF_F_V_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063191")
# And associated MetaDeta, e.g., session and subject information:
GRF_metadata <- read_csv("https://ndownloader.figshare.com/files/22062960")
#.. you might have to wait.. it Will download!



# 2) Data Extraction ---------------------------------------------------------

# Going to do a simple FLM 
#
# Class Labels as Scalars
#
# We will use the pre-specified test_train split in the GaitRec data set.
#
# Use balanced train data 
#
# But only left side of injured subjects. (so un-balanced vs control, who ahve no injured legs)
#
# Check that it's one session subject per subject when we do this
# even though we know it is
GRF_metadata %>% filter(CLASS_LABEL == "A") %>% # Ankle
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  group_by(SUBJECT_ID) %>%
  summarise(n = n()) %>%
  summarise(max(n), min(n))

GRF_metadata %>% filter(CLASS_LABEL == "K") %>% # knee
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  group_by(SUBJECT_ID) %>%
  summarise(n = n()) %>%
  summarise(max(n), min(n))

GRF_metadata %>% filter(CLASS_LABEL == "H") %>% # Hip
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  group_by(SUBJECT_ID) %>%
  summarise(n = n()) %>%
  summarise(max(n), min(n))

GRF_metadata %>% filter(CLASS_LABEL == "C") %>% # calcaneus
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  group_by(SUBJECT_ID) %>%
  summarise(n = n()) %>%
  summarise(max(n), min(n))

GRF_metadata %>% filter(CLASS_LABEL == "HC") %>% # healthy control
  filter(TRAIN_BALANCED == 1) %>%
  group_by(SUBJECT_ID) %>%
  summarise(n = n()) %>%
  summarise(max(n), min(n))



# Extract the session ID's for for the data set that
# we are going to use

Ankle_Sessions <- GRF_metadata %>% filter(CLASS_LABEL == "A") %>% # A = Ankle
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  select(SESSION_ID, CLASS_LABEL)

Knee_Sessions <- GRF_metadata %>% filter(CLASS_LABEL == "K") %>% # K = Knee
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  select(SESSION_ID, CLASS_LABEL)

Hip_Sessions <- GRF_metadata %>% filter(CLASS_LABEL == "H") %>% # H = hip
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  select(SESSION_ID, CLASS_LABEL)

Calcaneus_Sessions <- GRF_metadata %>% filter(CLASS_LABEL == "C") %>% # C = calcaneus
  filter(TRAIN_BALANCED == 1) %>%
  filter(AFFECTED_SIDE == 0) %>% # LEFT SIDE INJURED
  select(SESSION_ID, CLASS_LABEL)

Healthy_Sessions <- GRF_metadata %>% filter(CLASS_LABEL == "HC") %>% # HC = Helathy Control
  filter(TRAIN_BALANCED == 1) %>% 
  select(SESSION_ID, CLASS_LABEL) # No affected side because HC


Sessions_df <- bind_rows(Healthy_Sessions, Ankle_Sessions, Knee_Sessions, Hip_Sessions, Calcaneus_Sessions)


# 3) Data Processing/ Cleaning -----------------------------------------------

# We will subset it to look at only the session we have selected above
# Let's see how many trials per session/subject: 
GRF_F_V_PRO_left %>%
  filter(SESSION_ID %in% Sessions_df$SESSION_ID) %>%
  group_by(SESSION_ID) %>%
  summarise(num_trials = n_distinct(TRIAL_ID)) %>%
  summarise(min(num_trials), max(num_trials), median(num_trials))
# between 5 and 15 trials per subject, we will average them 
# before fitting our linear model
# (we want an iid sample so to stick to assumptions)

# the data are already smoothed ('filtered')
# think 35 basis functions and a tiny smoothing parameter will do.
# we will use cubic spline, penalize 2nd derivative with a tiny penalty
bsplines_35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 35, norder = 4)


# we will average all trials within a session to have 1 per subject/session
# therefore, a light registration within each subject will help this average
# smoothing applied to warping functions
# Set up basis and smoothing for warping functions:
# only 5 basis functions, not too much flexibility and will
# add a reasonable penalty to their smoothness.
wbasis <- create.bspline.basis(rangeval = c(0, 100), nbasis = 5, norder = 4)
wcoef <- matrix(0, 5, 1)
wb_obj <- fd(basisobj =  wbasis, coef = wcoef)
w_par <- fdPar(fdobj = wb_obj, Lfdobj = 2,lambda = 10^-4) 
# I've chosen a reasonably heavy penalty on the warping functions so to only do light warping


# Going to loop through every subject, and get an average, with some light warping.
# store the coefficients of the average for each subject in this matrix
coefficients_average <- matrix(data = NA, nrow = bsplines_35$nbasis, ncol = nrow(Sessions_df))

# Start loop
# we are going top loop through and check each plot
# to verify registration has worked before moving 
# to next subject

#par(mfrow = c(3, 1), cex = 1, ask = T) # set up for side-by-side plot 
# alternatively, we could save each plot and go through after
# (see commented code within loop)
for(sesh in 1:nrow(Sessions_df)){ 
  
  # See what number we are on
  print(paste("Subject", sesh, "of", nrow(Sessions_df)))
  sesh_id <- Sessions_df$SESSION_ID[sesh] # loop through each session
  GRF_data <- GRF_F_V_PRO_left %>% filter(SESSION_ID == sesh_id) # subset GRF data to contain session
  n_trials <- nrow(GRF_data) # number of trials
  trial_id <- GRF_data$TRIAL_ID # ID of the trials
  grf_fd <- smooth.basis(argvals = 0:100, # smooth the data with the 35 basis functions
                         y = t(GRF_data[, -c(1:3)]), # and a tiny roughness penalty (10^-10)
                         fdParobj = fdPar(fdobj = bsplines_35, Lfdobj = 2, lambda = 10^-10))
  # register the trials before computing an average, using the continous reg (Ramsay and Silverman (2005)
  registered_funs <- register.fd(yfd = grf_fd$fd, WfdParobj = w_par, crit = 2)
  
  #file_name <- paste0("Outputs/Diagnostic_Plots/", sesh, ".jpeg")
  # Make a plot to see if registration performed ok:
  #jpeg(filename = file_name, width = 900, height = 1000)
  
  #  Comment out if you don't want to see each plot #
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
  #  End comment out #
  
  # Store the coefficients of the average of the trials for this session.
  coefficients_average[, sesh] <- mean.fd(registered_funs$regfd)$coefs
}

# on inspecting each plot, registration barely seems to have changed each average curve
# maybe, not necessary, but did no harm.


# create an fd object from the averaged profiles:
combined_average_profiles <- fd(coef = coefficients_average, basisobj = bsplines_35)
combined_average_profiles$reps <- Sessions_df$SESSION_ID # add session names
combined_average_profiles$fdnames$funs <- "grf"

# we could plot each group mean, (commented out below but done in next scrip)

#par(mfrow = c(1, 1), ask = F)
#plot(mean.fd(combined_average_profiles[which(Sessions_df$CLASS_LABEL == "HC")]))
#lines(mean.fd(combined_average_profiles[which(Sessions_df$CLASS_LABEL == "A")]), col = "red")
#lines(mean.fd(combined_average_profiles[which(Sessions_df$CLASS_LABEL == "K")]), col = "blue")
#lines(mean.fd(combined_average_profiles[which(Sessions_df$CLASS_LABEL == "H")]), col = "green")
#lines(mean.fd(combined_average_profiles[which(Sessions_df$CLASS_LABEL == "C")]), col = "orange")
par(ask = F)

# Registration of all the curves could be done now to a single mean.
# in this case, we chose not to, as curves from differing groups had slighty 
# different shapes
# in contrast to, say, just a group of healthy controls

