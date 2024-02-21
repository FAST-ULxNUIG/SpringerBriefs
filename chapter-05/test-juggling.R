library(fda)

x <- readr::read_table2("chapter-02/juggling-data/Xcoord.txt", col_names = FALSE) # x coordinates
y <- readr::read_table2("chapter-02/juggling-data//Ycoord.txt", col_names = FALSE) # y coordinates
z <- readr::read_table2("chapter-02/juggling-data//Zcoord.txt", col_names = FALSE) # z co-ordinates

cols <- rainbow(ncol(x)) 
t <- (0:(nrow(x)-1)) * 5 # time arguments (approx. ms)
par(mfrow = c(3, 1))
matplot(t, as.matrix(x), type = "l", lty = 1, xlab = expression(t), ylab = expression(x(t)), col = cols)
grid(lty = 1)
title(expression(x(t)))
matplot(t, as.matrix(y), type = "l", lty = 1, xlab = expression(t), ylab = expression(y(t)), col = cols)
grid(lty = 1)
title(expression(y(t)))
matplot(t, as.matrix(z), type = "l", lty = 1, xlab = expression(t), ylab = expression(z(t)), col = cols)
grid(lty = 1)
title(expression(z(t)))


log10_lambda_grid <- seq(0, 5, by = 0.5)
order_6_basis <- create.bspline.basis(rangeval = range(t), nbasis = length(t), norder = 6, )
gcv_vec <- vector(mode = "numeric", length = length(log10_lambda_grid))
  
for(i in seq_along(log10_lambda_grid)) {
  gcv_vec[i] <- smooth.basis(argvals = t, y = as.matrix(z)[,1], fdParobj = fdPar(fdobj = order_6_basis, Lfdobj = 4, lambda = 10^log10_lambda_grid))
}
z1_fd <- Data2fd(argvals = t, )
plot(z1_fd, )

z_deriv <- deriv.fd(z1_fd)
plot(eval.fd(t, z_deriv), eval.fd(t, z1_fd), type = "l")
