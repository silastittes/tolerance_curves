#setwd("")
source("functions.R")

simK <- function(N,ndims,lambda,startree,niters=1000)
{
  x_eval <- seq(0,15,length=ndims)
  pmat <- Kmat <- matrix(NA,nrow = niters,length(x_eval)+1)
  j <- 0
  counter <- 0
  while(counter<niters)
  {
    j <- j+1
    print(j)
    data <- try(f(seed = j,N = N,startree = startree,lambda = lambda),silent=TRUE)
    if(class(data)=="try-error") next
    univar <- vector("list",length(x_eval))
    mvar <- try(multivar_phylosig(data$wilting$true$tip_coefficients,data$tree,length(x_eval),999),silent=TRUE)
    if(class(mvar)=="try-error") next
    for(i in 1:length(x_eval)) univar[[i]] <- try(phylosig(data$tree,apply(data$wilting$true$tip_coefficients,1,function(X) logit_fx(X,x_eval[i])),test = TRUE,nsim = 249),silent=TRUE)
    if(any(sapply(univar,function(X) class(X)=="try-error"))) next
    Kmat[counter+1,1] <- mvar[[1]]
    Kmat[counter+1,1:length(x_eval)+1] <- sapply(univar,function(X) X[[1]])
    pmat[counter+1,1] <- mvar[[2]]
    pmat[counter+1,1:length(x_eval)+1] <- sapply(univar,function(X) X[[2]])
    colnames(pmat) <- colnames(Kmat) <- c("Kmult",x_eval)
    counter <- counter+1
    Pfilename <- paste("Kmult_P_N",N,"_dims",ndims,"_","type",abs(as.integer(startree)-2),"_lambda",lambda,".csv",sep="")
    Kfilename <- paste("Kmult_K_N",N,"_dims",ndims,"_","type",abs(as.integer(startree)-2),"_lambda",lambda,".csv",sep="")
    write.table(pmat,Pfilename,sep=",",row.names = FALSE)
    write.table(Kmat,Kfilename,sep=",",row.names = FALSE)
  }
}

N <- c(32,128) # number of species
ndims <- c(2,5,10,25,50,100,250) # number of dimensions
lambda <- c(0,.05,.1,.2,.4,.6,.8,1)
startree <- c(FALSE)
niters <- 1000 # number of simulations

for(loop_N in 1:length(N))
{
  for(loop_ndims in 1:length(ndims))
  {
    for(loop_lambda in 1:length(lambda))
    {
      for(loop_startree in 1:length(startree))
      {
        if((!startree) | (startree & lambda==1))
        {
          simK(N=N[loop_N],ndims=ndims[loop_ndims],lambda=lambda[loop_lambda],startree = startree[loop_startree],niters = niters)
        }
      }
    }
  }
}