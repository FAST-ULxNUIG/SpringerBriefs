library(readr)
x <- read_table2("juggling/Xcoord.txt", col_names = FALSE) # x coordinates
y <- read_table2("juggling/Ycoord.txt", col_names = FALSE) # y coordinates
z <- read_table2("juggling/Zcoord.txt", col_names = FALSE) # z co-ordinates
ni <- as.vector(read_csv("juggling/ni.txt", col_names = FALSE)) # no. of cycles in trial
seqn <- as.vector(t(read_table2("juggling/seqn.txt", col_names = FALSE)[1, ])) # no. of observations

time <- 5 * ( 0:(nrow((x)-1)) ) # calculate time (ms), data sampled at 200Hz = once every 5 miliseconds

juggling <- array(data = NA, dim = c(nrow(x), ncol(x), 3))
juggling[,, 1] <- as.matrix(x)
juggling[,, 2] <- as.matrix(y)
juggling[,, 3] <- as.matrix(z)



par(mfrow = c(3, 1), cex = 1.05)
matplot(time[1:1500], x[1:1500, ], type = "l", lty = 1, xlab = "Time (milliseconds)", ylab = "x(t)")
matplot(time[1:1500], y[1:1500, ], type = "l", lty = 1, xlab = "Time (milliseconds)", ylab = "y(t)")
matplot(time[1:1500], z[1:1500, ], type = "l", lty = 1, xlab = "Time (milliseconds)", ylab = "z(t)")

