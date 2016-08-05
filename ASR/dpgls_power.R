#setwd("")
source("functions.R")

simDPGLS <- function(N,ndims,covariation,niters=1000)
{
  x_eval <- seq(0,15,length=ndims)
  pmat <- R2mat <- matrix(NA,nrow = niters,length(x_eval)+1)
  j <- 0
  counter <- 0
  while(counter<niters)
  {
    j <- j+1
    print(j)
    data <- try(f(seed = j,N = N,extra_trait = "univariate",covariation = covariation,covary_with_slope = TRUE),silent=TRUE)
    extra_trait <- data$extra_trait
    extra_trait_pic <- pic(extra_trait,data$tree)
    if(class(data)=="try-error") next
    univar <- vector("list",length(x_eval))
    mvar <- try(multivar_pgls(x = data$wilting$true$tip_coefficients,univariate_trait = data$extra_trait,tree = data$tree,length(x_eval),999),silent=TRUE)
    if(class(mvar)=="try-error") next
    for(i in 1:length(x_eval)) 
    {
      univar_trait <- apply(data$wilting$true$tip_coefficients,1,function(X) logit_fx(X,x_eval[i]))
      uni_pic <- try(pic(univar_trait,data$tree),silent=TRUE)
      univar_model <- try(summary(lm(uni_pic~extra_trait_pic-1)),silent=TRUE)
      univar[[i]] <- list(R2=univar_model$r.squared,p=univar_model$coefficients[1,4])
    }
    if(any(sapply(univar,function(X) class(X)=="try-error"))) next
    R2mat[counter+1,1] <- mvar[[1]]$Rsq[1]
    R2mat[counter+1,1:length(x_eval)+1] <- sapply(univar,function(X) X[[1]])
    pmat[counter+1,1] <- mvar[[1]]$P.val[1]
    pmat[counter+1,1:length(x_eval)+1] <- sapply(univar,function(X) X[[2]])
    colnames(pmat) <- colnames(R2mat) <- c("DPGLS",x_eval)
    counter <- counter+1
    Pfilename <- paste("DPGLSslope_P_N",N,"_dims",ndims,"_covariation",covariation,".csv",sep="")
    R2filename <- paste("DPGLSslope_R2_N",N,"_dims",ndims,"_covariation",covariation,".csv",sep="")
    write.table(pmat,Pfilename,sep=",",row.names = FALSE)
    write.table(R2mat,R2filename,sep=",",row.names = FALSE)
  }
}