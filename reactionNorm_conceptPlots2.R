
my_alpha <- function(col, alpha){
  hexed <- col2rgb(col)
  return(
    rgb(hexed[1,1], hexed[2,1], hexed[3,1], alpha = (alpha)*255,
        maxColorValue = 255)
  )
}

stretch.kumara <- function(x, a, b, c){
  return(c*((a*b*x^(a-1) ) * (1-x^a)^(b-1)))
}


scale.kumara2 <- function(x, a, b, c, d, e){
  xs <- (x - d)/(e - d)
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  return(list(xs = xs, mod.fit=mod.fit))
}


unscale.kumara <- function(x, a, b, c, d, e){
  xs <- x*(e - d) + d
  mod.fit <- c*((a*b*x^(a-1) ) * (1-x^a)^(b-1))
  return(list(xs = xs, mod.fit=mod.fit))
}


#conceptual plot of continuous and binary reaction norm sampling

#-----------------interaction scenerio----------------
pdf("figures/fig1.pdf")
#png(filename = "analyses_and_viz/conceptual_ReNorm.png")
par(mfrow=c(2,2))
par(mar=c(2,4,4,1)*1.1)

xseq <- seq(0, 1, length.out = 500)
a <- c(3, 4)
b <- c(4, 4)
c <- c(1.1, 0.8)
d <- c(4, 0)
e <- c(14.5, 10)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.9, 1.1)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)

plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="", ylab="Trait or Fitness Value", main="")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

liwd <- 5
lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(6.25, 9.2, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

#pty <- 22
pty <- 21
pcx <- 2

alph <- 0.7

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=my_alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=my_alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = my_alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = my_alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = my_alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = my_alpha("black", alph), cex=pcx)


#-------------false parallel scenerio-----------------
par(mar=c(2,2,4,3)*1.1)

xseq <- seq(0, 1, length.out = 500)
a <- c(4, 4)
b <- c(4, 4)
c <- c(1, 0.7)
d <- c(6, 3)
e <- c(14, 25)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.8, 1.2)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)

plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="", ylab="")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(9, 12.55, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=my_alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=my_alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = my_alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = my_alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = my_alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = my_alpha("black", alph), cex=pcx)



#-------------specialist generalist-----------------
par(mar=c(5,4,1,1)*1.1)

xseq <- seq(0, 1, length.out = 500)
a <- c(4, 4)
b <- c(4, 4)
c <- c(1, 0.65)
d <- c(3.75, 1)
e <- c(7.5, 9)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.8, 1.2)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)


plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="Environment", ylab="Trait or Fitness Value")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(4.85, 6.25, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=my_alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=my_alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = my_alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = my_alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = my_alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = my_alpha("black", alph), cex=pcx)


#-------------always lower always higher parallel-----------------
par(mar=c(5,2,1,3)*1.1)
xseq <- seq(0, 1, length.out = 500)
a <- c(4, 4)
b <- c(4, 4)
c <- c(1, 0.6)
d <- c(2, 4)
e <- c(10, 9)

s1 <- unscale.kumara(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- unscale.kumara(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

yrng <- range(c(s1$mod.fit, s2$mod.fit))*c(0.8, 1.2)
xrng <- range(c(s1$xs, s2$xs))*c(0.9, 1.1)


plot(NA, NA, xlim=xrng, ylim=yrng, axes=F, cex.lab=1.2,
     xlab="Environment", ylab="")
box()
axis(side = 1, labels = F)
axis(side = 2, labels = F)

lines(s1$xs, s1$mod.fit, lwd=4, lty=2)
lines(s2$xs, s2$mod.fit, lwd=4, lty=2, col="grey60")

xseq <- seq(5.5, 6.75, length.out = 500)
s1 <- scale.kumara2(x = xseq, a = a[1], b = b[1], c = c[1], d = d[1], e = e[1])
s2 <- scale.kumara2(x = xseq, a = a[2], b = b[2], c = c[2], d = d[2], e = e[2])

segments(x0 = xseq[1], y0 = s1$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s1$mod.fit, 1), col=my_alpha("grey60", alph), lwd=liwd)

segments(x0 = xseq[1], y0 = s2$mod.fit[1],
         x1 = tail(xseq, 1), y1 = tail(s2$mod.fit, 1), col=my_alpha("black", alph), lwd=liwd)

points(xseq[1], s1$mod.fit[1], bg = my_alpha("grey60", alph), pch=pty, cex=pcx)
points(tail(xseq, 1), tail(s1$mod.fit, 1), bg = my_alpha("grey60", alph), pch=pty, cex=pcx)

points(xseq[1], s2$mod.fit[1], pch=pty, bg = my_alpha("black", alph), cex=pcx)
points(tail(xseq, 1), tail(s2$mod.fit, 1),pch=pty, bg = my_alpha("black", alph), cex=pcx)


dev.off()
