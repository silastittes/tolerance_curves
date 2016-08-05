library(parallel)

ndraw <- length(posts_Normal$lp__)*1

xseq <- seq(0,1, length.out=500)

plot.kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  list(x, mod.fit)
}

opt.kumara <- function(xs, a, b, c, d, e1){
  x <- xs * e1 + d
  mod.fit <- c*((a*b*xs^(a-1) ) * (1-xs^a)^(b-1))
  x[which.max(mod.fit)]
}

######################
######################
##### WORKING!!! #####
######################
######################

xseq <- seq(0,1, length.out=500)
draw_list <- as.list(1:(ndraws*1))
sppint_list <- as.list(1:max(emery$sppint))
post_plot_draws <- lapply(sppint_list, function(sp)
            lapply(draw_list, function(pr)
                plot.kumara(xs = xseq,
                            a = posts_Normal$a[pr,sp],
                            b = posts_Normal$b[pr,sp],
                            c = posts_Normal$c[pr,sp],
                            d = posts_Normal$d[pr,sp],
                            e1 = posts_Normal$e1[pr,sp]
                            )
                   )
               )
names(post_plot_draws) <- unique(emery$Species)

spp <- "debilis"
xr <- range(sapply(post_plot_draws[[spp]], function(x) x[[1]]))
yr <- range(sapply(post_plot_draws[[spp]], function(x) x[[2]]))
plot(NA,NA, xlim=xr, ylim=yr)
lapply(post_plot_draws[[spp]], function(x) lines(x[[1]], x[[2]]))


xseq_big <- seq(0,1, length.out=10000)
draw_list <- as.list(1:(ndraws*1))
sppint_list <- as.list(1:max(emery$sppint))
opt <- mclapply(sppint_list, function(sp)
  sapply(draw_list, function(pr)
    opt.kumara(xs = xseq_big,
               a = posts_Normal$a[pr,sp],
               b = posts_Normal$b[pr,sp],
               c = posts_Normal$c[pr,sp],
               d = posts_Normal$d[pr,sp],
               e1 = posts_Normal$e1[pr,sp]
    )
  )
)

names(opt) <- unique(emery$Species)
maximadf <- data.frame(opt)
write.table(x = maximadf, file = "bayes/maxima_draws.txt")
