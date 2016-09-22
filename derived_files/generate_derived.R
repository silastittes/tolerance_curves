#GENERATE DERIVED PARAMETERS AND ANALYSES FOR LASTHENIA TOLERANCE CURVE DATA

source("load_data.R")
emery <- load_emery()

#Approximate location of water treatment optima
library(parallel)
xseq_big <- seq(0,1, length.out=10000)
draw_list <- as.list(1:(ndraws*1))
sppint_list <- as.list(1:max(emery$sppint))
opt <- mclapply(sppint_list, function(sp){
  sapply(draw_list, function(pr){
    opt.kumara(xs = xseq_big,
               a = posts$a[pr,sp],
               b = posts$b[pr,sp],
               c = posts$c[pr,sp],
               d = posts$d[pr,sp],
               e1 = posts$e1[pr,sp]
    )
    library(parallel)
    xseq_big <- seq(0,1, length.out=10000)
    draw_list <- as.list(1:(ndraws*1))
    sppint_list <- as.list(1:max(emery$sppint))
    opt <- mclapply(sppint_list, function(sp){
      sapply(draw_list, function(pr){
        opt.kumara(xs = xseq_big,
                   a = posts$a[pr,sp],
                   b = posts$b[pr,sp],
                   c = posts$c[pr,sp],
                   d = posts$d[pr,sp],
                   e1 = posts$e1[pr,sp]
        )
      }
      )
    }
    )
    
    names(opt) <- unique(emery$Species)
    maximadf <- data.frame(opt)
    write.table(x = maximadf, file = "derived_files/maxima_draws.txt")
    
    mden <- apply(maximadf, 2, density)
    ylimit <- max(sapply(mden, function(x) max(x$y)))
    xlimit <- range(sapply(mden, function(x) range(x$x)))
    plot(NA,NA, xlim=c(0,ylimit), ylim=xlimit, xlab="", ylab="", main="")
    sapply(mden, FUN = function(x) lines(x = x$y, y = x$x) )
  }
  )
}
)

names(opt) <- unique(emery$Species)
maximadf <- data.frame(opt)
write.table(x = maximadf, file = "derived_files/maxima_draws.txt")

mden <- apply(maximadf, 2, density)
ylimit <- max(sapply(mden, function(x) max(x$y)))
xlimit <- range(sapply(mden, function(x) range(x$x)))
plot(NA,NA, xlim=c(0,ylimit), ylim=xlimit, xlab="", ylab="", main="")
sapply(mden, FUN = function(x) lines(x = x$y, y = x$x) )
