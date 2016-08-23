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

par(mfrow=c(2,2))
par(mar=c(2,4,4,2)*1.1)
plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.5,
     xlab="", ylab="Trait or Fitness Value", main="Continuous (Multi Env.)")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

liwd <- 5
lines(s1$xs, s1$mod.fit, lwd=liwd)
lines(s2$xs, s2$mod.fit, lwd=liwd, col="grey60")

xseq <- seq(12, 15, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

#pty <- 22
pty <- 21
pcx <- 2


points(xseq[1], s1$mod.fit[1], bg="black", pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg="black", pch=pty, cex=pcx)
points(xseq[1], s2$mod.fit[1], pch=pty, bg="grey60", cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1), pch=pty, bg="grey60", cex=pcx)

#lines(xseq, s1$mod.fit, lwd=8)
#lines(xseq, s2$mod.fit, lwd=8, col="grey60")


yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.8,1.2)
plot(NA, NA, xlim=range(xseq)*c(0.95, 1.05), ylim=yrng, axes=F, 
     cex.lab=1.5, xlab="", ylab="",  main="Discrete (Two Env.)")
box()
axis(side = 1, at = range(xseq), labels = F)
axis(side = 2, labels = F)

alph <- 0.7
segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col="grey60", lwd=liwd)

lines(xseq, s1$mod.fit, col=alpha("black", alph), lty=2)
lines(xseq, s2$mod.fit, col=alpha("black", alph), lty=2)

points(xseq[1], s1$mod.fit[1], bg="black", pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg="black", pch=pty, cex=pcx)
points(xseq[1], s2$mod.fit[1], pch=pty, bg="grey60", cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1), pch=pty, bg="grey60", cex=pcx)

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

par(mar=c(5,4,2,2)*1.1)
plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.5,
     xlab="Environment", ylab="Trait or Fitness Value")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

lines(s1$xs, s1$mod.fit, lwd=liwd)
lines(s2$xs, s2$mod.fit, lwd=liwd, col="grey60")

xseq <- seq(11, 16, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e1 = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e1 = e[2])

points(xseq[1], s1$mod.fit[1], bg="black", pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg="black", pch=pty, cex=pcx)
points(xseq[1], s2$mod.fit[1], pch=pty, bg="grey60", cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg="grey60", cex=pcx)

#lines(xseq, s1$mod.fit, lwd=8)
#lines(xseq, s2$mod.fit, lwd=8, col="grey60")

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0, 1.2)

plot(NA, NA, xlim=range(xseq)*c(0.95, 1.05), ylim=yrng, axes=F, 
     xlab="Environment", ylab="", cex.lab=1.5)
box()
axis(side = 1, at = range(xseq), labels = F)
axis(side = 2, labels = F)


segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col="grey60", lwd=liwd)

lines(xseq, s1$mod.fit, col=alpha("black", alph), lty=2)
lines(xseq, s2$mod.fit, col=alpha("black", alph), lty=2)

points(xseq[1], s1$mod.fit[1], bg="black", pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg="black", pch=pty, cex=pcx)
points(xseq[1], s2$mod.fit[1], pch=pty, bg="grey60", cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg="grey60", cex=pcx)


dev.off()


