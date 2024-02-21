# ======================================================================= #
# Create two figures for exploratory data analysis section #
# 1) Normal boxplot with annotation (for comparison)
# 2) Functional boxplot
# ======================================================================= #

# Load neccesary packages and data ----------------------------------------

# Packages:
library(fda) # CRAN v5.1.5
library(tidyverse) # CRAN v1.3.1
library(funData) # CRAN v1.3-5

data.path <- here::here("chapter-03", "data")

# Datasets, already converted to functions:
average_profiles_reg <- readRDS(
  file = file.path(data.path, "GaitRecHealthyReg.RDS"))

selected_sessions <- readRDS(
  file = file.path(data.path, "SelectedSessions.RDS"))



# Manual data analysis ----------------------------------------------------

# Evaluate functional data on a grid of points for plotting.
average_session_reg_eval <- eval.fd(
  evalarg = 0:100,
  fdobj = average_profiles_reg$regfd)

# we will want to make a function boxplot of the data.
# can be done for discrete data using fbplot function
fb_plot_output <- fbplot(fit = average_session_reg_eval,
                         ylim = c(0, 1.5))
# note that we 'save' the plot with an assignment statement
# because some useful output (the depth values) is also printed with the plot
fb_plot_output

# can do a short manual recreation of the plot to gain intuition

# plot median curve to get us started
# which one is it
fb_plot_output$medcurve
plot(0:100,
     average_session_reg_eval[, fb_plot_output$medcurve],
     ylim = c(0, 1.5), # these limits make no sense for force data but just to compare output
     ylab = "vgrf",
     xlab = "normalised time",
     type = "l")

# now let's try to make 50% reason
# we need top n/2 curves ranked by depth and then
# at each value t get the minumum and maximum value

# first get total number of curves"
n <- ncol(average_session_reg_eval)
n_over_2 <- ceiling(n/2) # use in case of odd number

# create data frame with both curve values and depth values
data_frame_w_depth <- data.frame(
  t(average_session_reg_eval),
  depth = fb_plot_output$depth
)

names(data_frame_w_depth)[1:101] <- paste(0:100)

# use tidyverse tools to 
# order by depth in descending order
# and choose first n/2
n_over_2_df <- data_frame_w_depth %>%
  arrange(- depth) %>%
  slice(1:n_over_2) %>%
  select( - depth)

# create enevelopes
max_50_envelope <- apply(n_over_2_df, MARGIN = 2, max) %>%
  unlist()

min_50_envelope <- apply(n_over_2_df, MARGIN = 2, min) %>%
  unlist()

# and shade

polygon(c(0:100,
          rev(0:100)),
        c(min_50_envelope,
          rev(max_50_envelope)),
        col="deeppink",
        border = "blue")
# need to overlay median curve again

lines(0:100,
      lwd = 2,
      average_session_reg_eval[, fb_plot_output$medcurve])

# now for minimum and maximum fences
# get iqr width at each point t

iqr <- max_50_envelope - min_50_envelope


med_curve <- unlist(average_session_reg_eval[, fb_plot_output$medcurve])

# The minimum and maximum fence cut-offs are the min and max 50% envelopes
# extended by 1.5 times the IQR
min_envelope_cutoff <- min_50_envelope - 1.5 * iqr

max_envelope_cutoff <- max_50_envelope + 1.5 * iqr

# Now. we check does each curve check this point at any point t in the domain
# set up vector to store if it dies, yes or no
curve_exclude <- logical(length = ncol(average_session_reg_eval))

# loop through each curve (column of the matrix)
# can't think of a better way to do this.

for(i in 1:ncol(average_session_reg_eval)){
  
  # curve is a column of data frame - subset it
  curve <- average_session_reg_eval[, i]
  
  # fill in true if breaks the upper or lower fence and false otherwise
  # (the TRUE values will be outlier curves and removed from our calculations for min and max fence)
  curve_exclude[i] <- any(curve < min_envelope_cutoff | curve > max_envelope_cutoff)
}

# calculate the enevelopes, the se are max and min curves
# with outliers excluded
max_envelope <- apply(average_session_reg_eval[, !curve_exclude], 1, max)
min_envelope <- apply(average_session_reg_eval[, !curve_exclude], 1, min)

# draw them on
lines(0:100, min_envelope, col = "blue")
lines(0:100, max_envelope, col = "blue")


# finally, add outliers
matlines(0:100, average_session_reg_eval[, curve_exclude],
        col = 'red', type = "l", lty = 2)

#... looks like a good recreation

# could also use boxplot.fd to plot the fd objects directly
boxplot.fd(x = average_profiles_reg$regfd)

# ^ for our figure, we use fbplot() as it gives a little more freedom with
# (see below)

# for fancy, annotated figures... see below:
# The only reason that I've written them into functions
# is to save me running line by line (or highlighting chunks)

# Normal boxplot ----------------------------------------------------------
# with annotations

# choose the functional data set at a single point to create scalar data
chosen_time_point <- 50
pointwise_vec <- average_session_reg_eval[chosen_time_point, ]


library(tikzDevice)

cex_val <- 1.2

create_normal_boxplot <- function(){
  
  tikz(file  = here::here("chapter-03", "figures","normal-boxplot.tex"),
       width = 550 * (1/95),
       height = 600* (1/95),
       standAlone = TRUE)
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.8, 4.1, 2.1))
  q1 <- quantile(pointwise_vec, 0.25)
  med <- median(pointwise_vec)
  q3 <- quantile(pointwise_vec, 0.75)
  
  mybp <- boxplot(pointwise_vec,
                  cex=cex_val,
                  medcolor = "black",
                  medlwd = 4,
                  col = "#cd2fbc",
                  border = "#50abea",
                  outcol = "red",
                  outwex = 10, 
                  whiskcol = "darkgrey",
                  whisklwd = 2,
                  ylab = "vGRF (Normalised to BW) at $t=50\\%$",
                  cex.lab = 1.3 * cex_val,
                  staplelwd = 3,
                  boxlwd = 3,
                  boxcol= "#50abea",
                  outpch = 16,
                  xlim = c(0.5, 1.8),
                  ylim = c(0.4, 1),
                  cex.main = 1.5 * cex_val,
                  main = "Boxplot for Scalar Data")
  
  arrows(x0 = 0.75, x1 = 0.75,
         y0 = q1,
         y1 = q3,
         cex=cex_val,
         col = "#cd2fbc",
         length =  0.05)
  arrows(x0 = 0.75, x1 = 0.75,
         y1 = q1,
         y0 = q3,
         cex=cex_val,
         col = "#cd2fbc",
         length =  0.05)
  
  text(x = 0.6,
       y = med,
       cex=cex_val,
       col = "#cd2fbc",
       "Interquartile\nRange (IQR)")
  
  arrows(
    y0 = median(pointwise_vec),
    x1 = 1.25,
    cex=cex_val,
    x0 = 1.4,
    lwd = 2,
    length =  0.05
  )
  
  text(x = 1.5,
       y = med,
       cex=cex_val,
       "Median")
  
  arrows(y0 = q1,
         x1 = 1.25,
         x0 = 1.4,
         cex=cex_val,
         lwd = 2,
         length =  0.05,
         col = "#50abea")
  
  text(x = 1.6,
       y = q1,
       cex=cex_val,
       "Q1 First Quartile",
       col = "#50abea")
  
  arrows(y0 = q3,
         x1 = 1.25,
         x0 = 1.4,
         cex=cex_val,
         lwd = 2,
         length =  0.05,
         col = "#50abea")
  
  text(x = 1.6,
       y = q3,
       cex=cex_val,
       "Q3 Third Quartile",
       col = "#50abea")
  
  lower_fence <- c(mybp$stats)[1]
  upper_fence <- c(mybp$stats)[5]
  
  arrows(y0 = lower_fence,
         x1 = 1.15,
         x0 = 1.35,
         lwd = 2,
         cex=cex_val,
         length =  0.05,
         col = "#50abea")
  
  text(x = 1.6, y = lower_fence,
       cex=cex_val,
       "Smallest data value\n $\\geq$ Q1 - (1.5) IQR",
       col = "#50abea")
  
  arrows(y0 = upper_fence,
         x1 = 1.15,
         x0 = 1.35,
         lwd = 2,
         cex=cex_val,
         length =  0.05,
         col = "#50abea")
  
  text(x = 1.6, y = upper_fence,
       "Largest data value\n $\\leq$ Q3 + (1.5) IQR",
       col = "#50abea", cex=cex_val)
  
  outliers <- mybp$out
  
  arrows(x1 = 1.1,
         x0 = 1.35,
         cex=cex_val,
         y1 = mean(outliers),
         y0 = 0.5,
         lwd = 2,
         length =  0.05,
         col = "red")
  
  text(y = 0.5, 
       x = 1.575,
       cex=cex_val,
       "Outlying values",
       col = "red")
  
  points(x = 1,
         y = mean(outliers),
         pch = 1,
         cex = 6,
         col = 'red')
  
  dev.off()
  tinytex::lualatex(here::here("chapter-03", "figures","normal-boxplot.tex"))
}


create_normal_boxplot()


# Functional Boxplot ------------------------------------------------------


create_functional_boxplot <- function(){
  png(filename = here::here(
    "chapter-03", "figures", "functional-boxplot.png"),
    pointsize = 15,
    width = 1000,
    height = 480)
  par(xpd = T, family = "serif")
  layout(mat = matrix(c(1, 2, 3), nrow = 1),
         widths = c(0.4, 0.1, 0.4))
  fbplot(average_session_reg_eval,
         0:100,
         ylim = c(0, 1.5),
         main = "",
         xlab = "Normalised Time (% of Stance)",
         ylab = "vGRF (Normalised to BW)",
         xlim = c(0, 100),
         xaxs="i",
         yaxs="i",
         cex.lab = 1.5
  )
  title(main = "Functional Boxplot",
        cex.main = 1.5)
  plot(x = 1,
       y = 1,
       cex = 1.5,
       type = "n",
       yaxt = "n",
       xaxt = "n",
       bty = "n",
       xlab = "",
       ylab = "")
  
  arrows(x0 = 1.35,
         x1 = 1.9,
         y0 = 1.05,
         y1 = 1.04,
         col = "#cd2fbc",
         length = 0.1
  )
  
  text(x = 0.5,
       y = 1.05,
       "Central\n50% Region",
       cex = 1.5,
       col = "#cd2fbc")
  
  arrows(x0 = 1.055,
         x1 = 1.85,
         y0 = 1.125,
         y1 = 1.06,
         col = "black",
         length = 0.1,
  )
  
  text(x = 0.695,
       y = 1.175,
       "Functional\nMedian",
       cex = 1.5,
       col = "black")
  
  fbplot(average_session_reg_eval,
         ylim = c(0.2, 1.2),
         xlim = c(45, 55),
         xaxs="i",
         yaxs="i",
         yaxt = "n",
         xaxt = "n",
         bty = "n", 
         xlab = "",
         ylab = "")
  
  lower_text_height <- 0.325
  upper_text_height <- 1.1
  upper_arrow_start <- 1.08
  
  text(x = 46,
       y = lower_text_height,
       cex = 1.5,
       "Minimum Envelope",
       col = "#50abea",
       xpd = T)
  
  arrows(x0 = 45,
         x1 = 45,
         y0 = 0.35,
         y1 = 0.575,
         col = "#50abea",
         length = 0.1
  )
  
  
  text(x = 46,
       y = upper_text_height,
       cex = 1.5,
       "Maximum Envelope",
       col = "#50abea",
       xpd = T)
  
  arrows(x0 = 45,
         x1 = 45,
         y1 = 0.9175,
         y0 = upper_arrow_start,
         col = "#50abea",
         length = 0.1
  )
  
  text(x = 50,
       y = lower_text_height,
       cex = 1.5,
       "Outlying curves",
       col = "red",
       xpd = T)
  
  arrows(x0 = 50,
         x1 = 50,
         y0 = 0.35,
         y1 = 0.5,
         col = "red",
         length = 0.1
  )
  
  
  text(x = 52,
       y = 1.1,
       cex = 1.5,
       "Upper 50% Envelope",
       col = "#50abea",
       xpd = T)
  
  arrows(x0 = 52,
         x1 = 52,
         y1 = 0.825,
         y0 = upper_arrow_start,
         col = "#50abea",
         length = 0.1
  )
  
  
  text(x = 54,
       y = lower_text_height - 0.025,
       cex = 1.5,
       "Lower\n50% envelope",
       col = "#50abea",
       xpd = T)
  
  arrows(x0 = 54,
         x1 = 54,
         y0 = 0.35,
         y1 = 0.65,
         col = "#50abea",
         length = 0.1
  )
  
  dev.off()
}

create_functional_boxplot()

