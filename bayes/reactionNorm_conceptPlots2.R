#conceptual plot of continuous and binary reaction norm sampling
setwd("~/Documents/Projects/toleranceCurves/bayes/")
source("tolerance_functions.R")
library(scales)

scale.kumara2 <- function(x, a, b, c, d, e1){
  xs <- (x - d)/e1
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  return(list(xs = xs, mod.fit=mod.fit))
}


#-----------------interaction scenerio----------------

xseq <- seq(0, 1, length.out = 500)
a <- c(4, 4)
b <- c(4, 4)
c <- c(1.4, 1)
d <- c(7, 1)
e <- c(17, 15)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.9, 1.1)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)

pdf("figures/conceptual_ReNorm.pdf")

par(mfrow=c(2,1))
par(mar=c(1,4,2,2)*1.1)
plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="", ylab="Trait or Fitness Value", main="")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

liwd <- 5
lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(12, 15, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

#pty <- 22
pty <- 21
pcx <- 2

alph <- 0.7

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = alpha("black", alph), cex=pcx)

#-------------false parallel scenerio-----------------

xseq <- seq(0, 1, length.out = 500)
a <- c(4, 4)
b <- c(4, 4)
c <- c(1, 0.7)
d <- c(5.8, 3)
e <- c(12.8, 30)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.8, 1.2)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)

par(mar=c(5,4,1,2)*1.1)
plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="Environment", ylab="Trait or Fitness Value")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(11, 16, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = alpha("black", alph), cex=pcx)


dev.off()


