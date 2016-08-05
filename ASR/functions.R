require(phytools)
require(geiger)
require(geomorph)

# logit function
# first parameter is the glm intercept
# second parameter the glm slope
logit_fx <- function(theta,x)
{
  b <- theta[1]
  m <- theta[2]
  exp(m*x+b)/(1+exp(m*x+b))
}

# inverse logit function
logit_inv <- function(theta,y)
{
  b <- theta[1]
  m <- theta[2]
  ((log(-y/(y-1))-b)/m)
}

# transform logit glm parameters to uncorrelated slope parameter
logit_slope_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  m/4
}

# transform logit glm parameters to uncorrelated EC50 parameter
logit_ec50_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  -b/m
}

# back-calculates glm intercept parameter
logit_b_func <- function(slope,ec50)
{
  m <- slope*4
  b <- -ec50 * m
  return(b)
}

# back-calculates glm slope parameter
logit_m_func <- function(slope,ec50)
{
  m <- slope*4
  return(m)
}

# simulates dataset (plant wilting) for N species over an environmental gradient (drought)
# with x_length levels and n replicates per level
f=function(seed=2,N=30,x_length=30,n=20,startree=FALSE,lambda=1,measurement_error_in_data=FALSE,extra_trait=FALSE,covariation=.8,covary_with_slope=FALSE)
{
  x <- seq(0,15,length=x_length) # environmental gradient
  set.seed(seed)
  tree <- pbtree(n=N)
  simtree <- rescale(tree,"lambda",lambda=lambda)
  if(startree) simtree <- starTree(tree$tip.label,rep(1,length(tree$tip.label)))
  if(extra_trait=="univariate" & !covary_with_slope)
  {
    logit_slope <- fastBM(simtree,.5,bounds=c(.2,1),internal=TRUE)
    logit_ec50 <- sim.corrs(simtree,vcv = matrix(c(1,covariation,covariation,1),2,2),anc = 5,internal = TRUE)
    extra_trait_val <- logit_ec50[1:length(tree$tip.label),2]
    logit_ec50 <- logit_ec50[,1]
    if(any(logit_ec50<2) | any(logit_ec50>20))
    {
      logit_ec50 <- scale_between(X = logit_ec50,lower = 2,upper = 20,start = 5)
    }
  } else if(extra_trait=="univariate" & covary_with_slope)
  {
    logit_ec50 <- fastBM(simtree,5,bounds=c(2,20),internal=TRUE)
    logit_slope <- sim.corrs(simtree,vcv = matrix(c(1,covariation,covariation,1),2,2),anc = .5,internal = TRUE)
    extra_trait_val <- logit_slope[1:length(tree$tip.label),2]
    logit_slope <- logit_slope[,1]
    if(any(logit_ec50<.2) | any(logit_ec50>1))
    {
      logit_slope <- scale_between(X = logit_slope,lower = .2,upper = 1,start = .5)
    }
  } else
  {
    logit_ec50 <- fastBM(simtree,5,bounds=c(2,20),internal=TRUE)
    logit_slope <- fastBM(simtree,.5,bounds=c(.2,1),internal=TRUE)
  }
  logit_b <- logit_b_func(slope=logit_slope,ec50=logit_ec50)
  logit_b_anc <- logit_b[(N+1):length(logit_b)]
  logit_b <- logit_b[1:N]
  logit_m <- logit_m_func(slope=logit_slope,ec50=logit_ec50)
  logit_m_anc <- logit_m[(N+1):length(logit_m)]
  logit_m <- logit_m[1:N]
  logit_coefs <- cbind(logit_b,logit_m)
  logit_sim_coefs <- logit_coefs
  logit_y <- matrix(nrow=length(x),ncol=N)
  for(i in 1:N)
  {
    if(measurement_error_in_data)
    {
      logit_y[,i] <- (rbinom(length(x),n,logit_fx(c(logit_b[i],logit_m[i]),x))/n)
    } else
    {
      logit_y[,i] <- logit_fx(c(logit_b[i],logit_m[i]),x)
    }
    temp_y <- logit_y[,i]*n
    if(measurement_error_in_data)
    {
      logit_sim_coefs[i,] <- suppressWarnings(glm(cbind(temp_y,n-temp_y)~x,family=binomial("logit"))$coefficients)
      rownames(logit_sim_coefs) <- rownames(logit_coefs)
    }
    colnames(logit_y) <- tree$tip.label
  }
  solveC <- solve(vcv(tree))
  D <- dist.nodes(tree)
  N <- length(tree$tip.label)
  Nnode <- tree$Nnode
  MRCA <- mrca(tree, full = TRUE)
  M <- D[as.character(N + 1), MRCA]
  dim(M) <- rep(sqrt(length(M)), 2)
  varAY <- M[(N+1):(N+Nnode), 1:N]
  logit_anc_coefs <- pgls_curve(tree = tree,tip_coefficients = logit_coefs,
             solveC = solveC,varAY = varAY)
  
  ret <- (list(tree=tree,
              wilting=list(
                observed=list(
                  water_deficiency=x,
                  raw_data=logit_y*n,
                  replicates=n,
                  tip_coefficients=logit_sim_coefs,
                  anc_coefficients=logit_anc_coefs),
                true=list(
                  tip_coefficients=cbind(logit_b,logit_m),
                  anc_coefficients=cbind(logit_b_anc,logit_m_anc)
                )
              )))
  if(extra_trait=="univariate")
  {
    ret[["extra_trait"]] <- extra_trait_val
  }
  ret
}

# wrapper for physignal
multivar_phylosig <- function(x,tree,ndims,iter)
{
  x <- t(apply(x,1,function(X) logit_inv(X,seq(.01,.99,length=ndims))))
  physignal_no_plot(tree,x,iter)
}

# wrapper for procD.pgls
multivar_pgls <- function(x,univariate_trait,tree,ndims,iter,error_type)
{
  x1 <- t(apply(x,1,function(X) logit_inv(X,seq(.01,.99,length=ndims))))
  if(missing(univariate_trait))
  {
    if(error_type==1)
    {
      univariate_trait <- fastBM(tree)
    } else if(error_type==2)
    {
      univariate_trait <- apply(x1,1,mean) + rnorm(length(tree$tip.label),0,.1)
    }
  }
  list(procd=procD.pgls(x1~univariate_trait,tree,iter),univariate_trait=univariate_trait)
}

pgls_curve <- function(tree,tip_coefficients,solveC,varAY,vals_only=FALSE)
{
  if(missing(solveC))
  {
    solveC <- solve(vcv(tree))
  }
  if(missing(varAY))
  {
    D <- dist.nodes(tree)
    N <- length(tree$tip.label)
    Nnode <- tree$Nnode
    MRCA <- mrca(tree, full = TRUE)
    M <- D[as.character(N + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    varAY <- M[(N+1):(N+Nnode), 1:N]
  }
  N <- length(tree$tip.label)
  Nnode <- tree$Nnode
  X <- matrix(1,nrow(solveC),1)
  y <- seq(.01,.99,length=100)
  Yinv <- t(apply(tip_coefficients,1,function(X) logit_inv(X,y)))
  root <- as.double(solve(t(X)%*%solveC%*%X)%*%(t(X)%*%solveC%*%Yinv))
  root <- (matrix(rep(root,N),length(root),N))
  anc_vals <- t(varAY%*%solveC%*%(Yinv-t(root)))+root[,1:Nnode]
  if(vals_only) return(t(anc_vals))
  ret <- t(apply(anc_vals,2,function(X) coef(glm(y~X,family=quasibinomial("logit")))))
  rownames(ret) <- (N+1):(N+Nnode)
  ret
}

# this function incorrectly reconstructs ancestral states of sigmoid curves using f(x) instead of f-1(y)
# for simulations only
pgls_curve_wrong <- function(tree,tip_coefficients,solveC,varAY,x_eval)
{
  if(missing(solveC))
  {
    solveC <- solve(vcv(tree))
  }
  if(missing(varAY))
  {
    D <- dist.nodes(tree)
    N <- length(tree$tip.label)
    Nnode <- tree$Nnode
    MRCA <- mrca(tree, full = TRUE)
    M <- D[as.character(N + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    varAY <- M[(N+1):(N+Nnode), 1:N]
  }
  N <- length(tree$tip.label)
  Nnode <- tree$Nnode
  X <- matrix(1,nrow(solveC),1)
  Y <- t(apply(tip_coefficients,1,function(X) logit_fx(X,x_eval)))
  root <- as.double(solve(t(X)%*%solveC%*%X)%*%(t(X)%*%solveC%*%Y))
  root <- (matrix(rep(root,N),length(root),N))
  anc_vals <- t(varAY%*%solveC%*%(Y-t(root)))+root[,1:Nnode]
  return(t(anc_vals))
}

# faster version of physignal from geomorph package by avoiding multiple calls to solve()
# also avoids plotting
physignal_no_plot <- function (phy, A, iter = 249, method = c("Kmult", "SSC")) 
{
  method <- match.arg(method)
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) {
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")
    }
    x <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(rownames(A))) {
      stop("Data matrix does not include taxa names as dimnames for rows.")
    }
    x <- A
  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  N <- length(phy$tip.label)
  if (N != dim(x)[1]) {
    stop("Number of taxa in data matrix and tree are not not equal.")
  }
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x))))) == 
        T) {
    stop("Names do not match between tree and data matrix.")
  }
  x <- x[phy$tip.label, ]
  if (is.null(dim(x)) == TRUE) {
    x <- matrix(x, dimnames = list(names(x)))
  }
  if (method == "Kmult") {
    Kmult <- function(x, phy,C,solveC,D.mat,K.denom) {
      x <- as.matrix(x)
      N <- length(phy$tip.label)
      ones <- array(1, N)
      C <- C[row.names(x), row.names(x)]
      solveC <- solveC[row.names(x), row.names(x)]
      D.mat <- D.mat[row.names(x), row.names(x)]
      a.obs <- colSums(solveC) %*% x/sum(solveC)
      distmat <- as.matrix(dist(rbind(as.matrix(x), a.obs)))
      MSEobs.d <- sum(distmat[(1:N), (N + 1)]^2)
      dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - (ones %*% a.obs))), 0)))
      MSE.d <- sum(dist.adj[(1:N), (N + 1)]^2)
      K.stat <- (MSEobs.d/MSE.d)/K.denom
      return(K.stat)
    }
    C <- vcv.phylo(phy)
    solveC <- solve(C)
    eigC <- eigen(C)
    D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% 
                     t(eigC$vectors))
    colnames(D.mat) <- rownames(D.mat) <- colnames(C)
    ones <- array(1, N)
    K.denom <- (sum(diag(C)) - N * solve(t(ones) %*% 
                                           solveC %*% ones))/(N - 1)
    K.obs <- Kmult(x, phy,C,solveC,D.mat,K.denom)
    P.val <- 1
    K.val <- rep(0, iter)
    for (i in 1:iter) {
      x.r <- as.matrix(x[sample(nrow(x)), ])
      rownames(x.r) <- rownames(x)
      K.rand <- Kmult(x.r, phy,C,solveC,D.mat,K.denom)
      P.val <- ifelse(K.rand >= K.obs, P.val + 1, P.val)
      K.val[i] <- K.rand
    }
    P.val <- P.val/(iter + 1)
    K.val[iter + 1] = K.obs
    if (dim(x)[2] > 1) {
      #plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = K.obs, pvalue = P.val))
  }
  if (method == "SSC") {
    anc.states <- NULL
    options(warn = -1)
    for (i in 1:ncol(x)) {
      tmp <- as.vector(fastAnc(phy, x[, i]))
      anc.states <- cbind(anc.states, tmp)
    }
    colnames(anc.states) <- NULL
    dist.mat <- as.matrix(dist(rbind(as.matrix(x), as.matrix(anc.states)))^2)
    SSC.o <- 0
    for (i in 1:nrow(phy$edge)) {
      SSC.o <- SSC.o + dist.mat[phy$edge[i, 1], phy$edge[i, 
                                                         2]]
    }
    P.val <- 1
    SSC.val <- rep(0, iter)
    for (ii in 1:iter) {
      x.r <- x[sample(nrow(x)), ]
      if (is.null(dim(x.r)) == TRUE) {
        x.r <- matrix(x.r)
      }
      row.names(x.r) <- row.names(x)
      anc.states.r <- NULL
      options(warn = -1)
      for (i in 1:ncol(x.r)) {
        tmp <- as.vector(fastAnc(phy, x.r[, i]))
        anc.states.r <- cbind(anc.states.r, tmp)
      }
      colnames(anc.states.r) <- NULL
      dist.mat.r <- as.matrix(dist(rbind(as.matrix(x.r), 
                                         as.matrix(anc.states.r)))^2)
      SSC.r <- 0
      for (i in 1:nrow(phy$edge)) {
        SSC.r <- SSC.r + dist.mat.r[phy$edge[i, 1], phy$edge[i, 
                                                             2]]
      }
      P.val <- ifelse(SSC.r <= SSC.o, P.val + 1, P.val)
      SSC.val[ii] <- SSC.r
    }
    P.val <- P.val/(iter + 1)
    SSC.val[iter + 1] = SSC.o
    if (dim(x)[2] > 1) {
      #plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = SSC.o, pvalue = P.val))
  }
}

# for constraining trait values from sim.corrs output
scale_between = function(X,lower,upper,start)
{
  lower <- max(start - (start - min(X)),lower)
  upper <- min(start + (max(X) - start),upper)
  ((upper - lower) * (X-min(X)))/(max(X)-min(X)) + lower
}

# for dynamic time warping
convert <- function(scaled1,scaled2,original_x)
{
  if(missing(scaled2))
  {
    scaled2 <- scaled1
  }
  min <- min(c(scaled1,scaled2))
  max <- max(c(scaled1,scaled2))
  a <- min(original_x)
  b <- max(original_x)
  new_x1 <- (((b-a)*(scaled1-min)) / (max - min)) + a
  new_x2 <- (((b-a)*(scaled2-min)) / (max - min)) + a
  return(list(new_x1,new_x2))
}

# calculates Felsenstein's ancestral state reconstruction values
ace_hand <- function(x,Y,tree,gpr_fit=FALSE)
{
  if(gpr_fit) x <- (x-min(x))/(max(x)-min(x))
  Y <- Y[tree$tip.label,]
  Y <- rbind(Y,matrix(0,nrow=tree$Nnode,ncol=ncol(Y)))
  v1_lookup <- v2_lookup <- n1 <- n2 <- d1 <- d2 <- double(1)
  corrected_d <- tree$edge.length
  for(i in (tree$Nnode+length(tree$tip.label)):(length(tree$tip.label)+1))
  {
    v1_lookup <- which(tree$edge[,1]==i)[1]
    v2_lookup <- which(tree$edge[,1]==i)[2]
    n1 <- tree$edge[v1_lookup,2] # node 1 number
    n2 <- tree$edge[v2_lookup,2] # node 2 number
    d1 <- (1-(corrected_d[v1_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    d2 <- (1-(corrected_d[v2_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    if(gpr_fit)
    {
      gp1 <- GP_fit(x,Y[n1,])
      gp2 <- GP_fit(x,Y[n2,])
      time_warp <- dtw(predict(gp1,x)$Y_hat,predict(gp2,x)$Y_hat)
    } else
    {
      gp1 <- coef(glm(Y[n1,]~x,family=quasibinomial("logit")))
      gp2 <- coef(glm(Y[n2,]~x,family=quasibinomial("logit")))
      time_warp <- dtw(logit_fx(gp1,x),logit_fx(gp2,x))
    }
    
    new_x1 <- (((max(x)-min(x))*(time_warp$index1-min(time_warp$index1))) / (max(time_warp$index1) - min(time_warp$index1))) + min(x)
    new_x2 <- (((max(x)-min(x))*(time_warp$index2-min(time_warp$index2))) / (max(time_warp$index2) - min(time_warp$index2))) + min(x)
    if(gpr_fit)
    {
      new_y1 <- predict(gp1,new_x1)$Y_hat
      new_y2 <- predict(gp2,new_x2)$Y_hat
    } else
    {
      new_y1 <- logit_fx(gp1,new_x1)
      new_y2 <- logit_fx(gp2,new_x2)
    }
    anc_x <- (d1*new_x1) + (d2*new_x2)
    anc_y <- (d1*new_y1) + (d2*new_y2)
    if(gpr_fit)
    {
      anc_gp <- GP_fit(anc_x,anc_y)
      Y[i,] <- predict(anc_gp,x)$Y_hat
    } else
    {
      anc_gp <- coef(glm(anc_y~anc_x,family=quasibinomial("logit")))
      Y[i,] <- logit_fx(anc_gp,x)
    }
    corrected_d[which(tree$edge[,2]==i)] <- corrected_d[which(tree$edge[,2]==i)]+((corrected_d[v1_lookup]*corrected_d[v2_lookup])/(corrected_d[v1_lookup]+corrected_d[v2_lookup]))
  }
  anc_y <- Y[(length(tree$tip.label)+1):(tree$Nnode+length(tree$tip.label)),]
  return(anc_y)
}

# calculates ancestral state reconstruction with pairwise dynamic time warping and 
# the fast PIC method (adapted from fastAnc in phytools)
#
# NOTE: generalized time warping is not currently implemented in R, but pairwise 
# alignment is implemented in the dtw package. DTW allows ancestral state
# reconstruction of any function using the fast PIC method (as in fastAnc),
# but does not align landmarks across species for distance-based multivariate
# methods such as K-mult (Adams 2014a) or D-PGLS (Adams 2014b).
# To perform alignment across species, generalized time warping may be performed
# using the current GTW MATLAB implementation: http://www.f-zhou.com/ta_code.html
# 
fast_anc_hand <- function(x,Y,tree,root_only=TRUE,gpr_fit=FALSE)
{
  if(is.null(rownames(Y)))
  {
    rownames(Y) <- tree$tip.label
    warning("No taxa labels attached to coefficients. Assuming order matches tips on tree.")
  } else
  {
    temp <- match(tree$tip.label,rownames(Y))
    Y <- Y[temp,]
  }
  
  M <- tree$Nnode
  N <- length(tree$tip.label)
  a <- tree
  node_vals <- matrix(0,M,ncol(Y))
  # calculate root state
  for(i in 1:M + N)
  {
    a <- multi2di(root(tree,node=(i)))
    node_vals[i-N,] <- ace_hand(x,Y,tree=a,gpr_fit = gpr_fit)[1,]
    if(root_only)
    {
      return(node_vals[1,])
    }
  }
  rownames(node_vals) <- 1:M + N
  return(node_vals)
}
