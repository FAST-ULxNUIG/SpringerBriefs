library(fda)
library(scatterplot3d)
library(tikzDevice)


plots_path <- here::here("chapter-05", "figures")

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

x <- readr::read_table2("chapter-02/juggling-data/Xcoord.txt", col_names = FALSE) # x coordinates
y <- readr::read_table2("chapter-02/juggling-data//Ycoord.txt", col_names = FALSE) # y coordinates
z <- readr::read_table2("chapter-02/juggling-data//Zcoord.txt", col_names = FALSE) # z co-ordinates



cols <- rainbow(ncol(x)) 
t <- (0:(nrow(x)-1)) # time arguments (approx. ms)
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


z_use <- z[,1, drop = TRUE]
z_use <- z_use[z_use!=0]
stopifnot(!any(is.na(z_use)))
t <- t[seq_len(length(z_use))]
plot(t, z_use)


z_fd <- Data2fd(y = z_use, argvals = t, basisobj = create.bspline.basis(norder = 6, rangeval = range(t), nbasis = round(length(t)/2)), lambda = 10000)

plot_df <- data.frame(t = t,
           )

z <- c(eval.fd(t, z_fd))
Dz <- c(eval.fd(t, z_fd, Lfdobj = 1))
D2z <- c(eval.fd(t, z_fd, Lfdobj = 2))


plot(z, D2z, type = "l")

xlim <- range(unlist(z), na.rm = TRUE)
ylim <- range(unlist(Dz), na.rm = TRUE)
zlim <- range(unlist(D2z), na.rm = TRUE)


tikz(file.path(plots_path, "juggling-phase-plane-plot.tex"),
     width = 0.5 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, 
     standAlone = TRUE)
par(mfrow = c(1, 1))
sp <- scatterplot3d(x = z,
                    y = Dz,
                    z = D2z,
                    type = "l",
                    cex.axis = 0.5,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    main = "Juggling",
                    zlab = "$D^2 z(t)$",
                    xlab = "$z(t)$",
                    angle = 10,
                    pch = 20,
                    cex.main = 1)
text(y = -2.25,  x = 7, "$Dz(t)$", xpd=TRUE, srt = 30, cex = 1)
dev.off()

tinytex::lualatex(file.path(plots_path, "juggling-phase-plane-plot.tex"))



