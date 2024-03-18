# =======================================================================#
# create a combination of different plots which show examples of multivariate
# functional data
# We use the gait data, the juggling data and a recent open source data set on running
# =======================================================================#

# 1) Import Packages ------------------------------------------------------
library(fda)              # CRAN v5.5.1
library(tidyverse)        # CRAN v1.3.1
library(ggpubr)           # CRAN v0.4.0
library(Polychrome)       # CRAN v1.5.1
library(rdryad)           # CRAN v1.0.0    
library(scatterplot3d)    # CRAN v0.3-41 
library(tikzDevice)       # CRAN v0.12.3.1

# 2) Settings -------------------------------------------------------------
source(here::here("functions", "theme_gunning.R"))
theme_gunning()
theme_update(axis.text = element_text(size = 10),
             plot.title = element_text(face = "bold", size = 13.5),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12))
plots_path <- here::here("chapter-05", "figures")

# 3) Gait Data: Read in and Wrangling ----------------------------------------------------------

# use the gait data set:
data('gait')
dimnames(gait)
dim(gait)
# 20 x 29 x 2 array
# rows are time points
gait_time <- as.numeric(dimnames(gait)[[1]])
# subjects are columns
gait_subject_names <- dimnames(gait)[[2]]



# 4) Gait Data: Smoothing ---------------------------------------------------------------

# use a periodic (Fourier basis)
# penalty on Second Derivative.
# choose lambda using a grid search and gcv

# Set up basis:
fourier_basis11 <- create.fourier.basis(rangeval = range(gait_time), nbasis = 11)

# Choose smoothing paramter values for grid search:
lambdas <- 10^seq(-12, 0, by = 0.1)
gcvs <- sapply(lambdas, function(y){
  smth <- smooth.basis(argvals = gait_time, y = gait,
               fdParobj = fdPar(fdobj = fourier_basis11, Lfdobj = 2, lambda = y))
  return(mean(smth$gcv))
})
# check grid search:
par(mfrow = c(1, 1))
plot(log10(lambdas), gcvs)
best_lambda <- lambdas[which.min(gcvs)]

print(paste("best value of lambda is", best_lambda))
# Smooth data using best lambda:
smooth_gait <- smooth.basis(argvals = gait_time, y = gait,
                            fdParobj = fdPar(fdobj = fourier_basis11,
                                             Lfdobj = 2, lambda = best_lambda))

# Plot data par(mfrow = c(2, 1))
plot(smooth_gait)



# 5) Gait Data: More Wrangling -------------------------------------------------------

# evaluate the data and put into a data frame:
eval_seq <- seq(min(gait_time), max(gait_time), length.out = 50)
eval_smooth_gait <- eval.fd(eval_seq, fdobj = smooth_gait$fd)
colnames(eval_smooth_gait) <- smooth_gait$fd$fdnames$reps

# joints (hip, knee) are in the array 'slices'
hip_df <- data.frame(time = eval_seq, eval_smooth_gait[,, 1])
knee_df <- data.frame(time = eval_seq, eval_smooth_gait[,, 2])
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




# 6) Gait Data: Plot -----------------------------------------------------------------

# make plot 
hip_vs_knee <- hip_knee_df %>%
  filter(subject %in% paste0("boy", 1:3)) %>%
  ggplot() +
  aes( x = hip_angle, y = knee_angle, group = subject, color = subject) +
  geom_path(lwd = 1) +
  geom_point() +
  labs(x = "Hip Angle ($^{\\circ}$)",
       title = "(b) Children's Gait",
       y = "Knee Angle ($^{\\circ}$)") + # axis labels
  theme(legend.position = "none") # don't need a legend

hip_vs_knee






# 7) Running Data: Import -------------------------------------------------

# Running data is from the paper:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244361#sec020
# Data set uploaded to:
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.1c59zw3sk
# we use the 'rdryad' package to download it directly:
all_angles <- read_csv(file = dryad_download("10.5061/dryad.1c59zw3sk")[[1]])

# the data can be used to show knee angles in running.


# 8) Running Data: Plot of the ankle angle and Knee ----------------------------------------
running_hip_ankle <- all_angles %>%
  filter(Treatment == "Body weight", Replicate == "BW1", subject_id %in% paste0("ID0", 5:8)) %>%
  ggplot() +
  aes(x = hip_angle, y = ankle_angle, group = subject_id, color = subject_id) +
  geom_path(lwd = 1) +
  geom_point() +
  labs(y = "Ankle Angle ($^{\\circ}$)", x = "Hip Angle ($^{\\circ}$)", title = "(c) Running") +
  theme(legend.position = "none")

running_hip_ankle



# 9) Juggling data   -------------------------------------------------------

## 9.1  Read in DataFrom Text Files -------------------------------------------------------
# Import from .txt files
juggling_data_path <- here::here("chapter-02", "juggling-data")
x <- read_table2(file.path(juggling_data_path,"Xcoord.txt"), col_names = FALSE) # x coordinates
y <- read_table2(file.path(juggling_data_path,"Ycoord.txt"), col_names = FALSE) # y coordinates
z <- read_table2(file.path(juggling_data_path,"Zcoord.txt"), col_names = FALSE) # z co-ordinates
ni <- as.vector(read_csv(file.path(juggling_data_path,"ni.txt"), col_names = FALSE))$X1 # no. of cycles in trial
seqn <- as.vector(t(read_table2(file.path(juggling_data_path,"seqn.txt"), col_names = FALSE)[1, ])) # no. of observations


## 9.2 Data Wrangling ------------------------------------------------------
# currently zeros are filled in for NA
x_t <- x
y_t <- y
z_t <- z
for(i in 1:length(ni)){
  if(ni[i] < nrow(x_t)){
    na_inds <- (ni[i]+1) : nrow(x_t)
    x_t[na_inds, i] <- y_t[na_inds, i] <- z_t[na_inds, i] <- NA
  }
}

# read in feature times to color the curves
Vfeaturetimes <- rmatio::read.mat(filename = file.path(juggling_data_path, "Vfeaturetimes.mat"))

# Fix clumsiness of reading Matlab files into R
Vfeaturetimes <- Vfeaturetimes$Vfeaturetimes
#  The tangential velocity features for each cycle as follows:
#  1. a minor local maximum in tangential velocity positioned between 
#     the minimum and the primary maximum (launch)
#  2. a major local maximum in tangential velocity (drop), and
#  3. a primary local minimum in tangential velocity (handoff)

dimnames(Vfeaturetimes)[[1]] <- c("launch", "drop", "handoff")

dimnames(Vfeaturetimes)

dim(Vfeaturetimes)
record_1_launches <- Vfeaturetimes[1, ,1]
# last record is not an actual timing
record_1_launches <- head(record_1_launches, -1)
ni_rec1 <- ni[1]
ni
time <- 5 * (0:(ni_rec1-1)) / 1000
color_cut <- cut(time, record_1_launches)


## 9.3 Plot ----------------------------------------------------------------

# open3d()
# plot3d(x_t$X1, y_t$X1, z_t$X1, type = "l", xlab = "x", ylab = "y", zlab = "z", col = as.numeric(color_cut), size = 5)
# plot3d(x_t$X1, y_t$X1, z_t$X1, add = T, col = as.numeric(color_cut), size = 6)
# par3d(windowRect = c(20, 30, 800, 800), cex = 1.5)
# save current angle viewpoint
#(view <- par3d("userMatrix"))
#par3d(userMatrix = view)


xlim <- range(unlist(x_t), na.rm = TRUE)
ylim <- range(unlist(y_t), na.rm = TRUE)
zlim <- range(unlist(z_t), na.rm = TRUE)

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "juggling-3d-plot.tex"),
     width = 0.5 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, 
     standAlone = TRUE)
par(mfrow=c(1, 1))
sp <- scatterplot3d(x = x_t[, 1, drop = TRUE],
                    y = y_t[, 1, drop = TRUE],
                    z = z_t[, 1, drop = TRUE],
                    type = "l",
                    cex.axis = 0.5,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    main = "(a) Juggling",
                    zlab = "$z(t)$",
                    xlab = "$x(t)$",
                    angle = 10,
                    pch = 20,
                    cex.main = 1)
# par("usr")
text(y = -2.5,  x = 7, "$y(t)$", xpd=TRUE, srt = 10, cex = 1)


for(j in seq_len(10)) {
  sp$points3d(x = x_t[, j, drop = TRUE],
              y = y_t[, j, drop = TRUE],
              z = z_t[, j, drop = TRUE],
              type = "l",
              col = j)
}

dev.off()
tinytex::lualatex(file.path(plots_path, "juggling-3d-plot.tex"))



tikz(file.path(plots_path, "gait-and-running-2d-plots.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, 
     standAlone = TRUE)
ggarrange(hip_vs_knee, running_hip_ankle, nrow = 1, ncol = 2)
dev.off()

tinytex::lualatex(file.path(plots_path, "gait-and-running-2d-plots.tex"))

