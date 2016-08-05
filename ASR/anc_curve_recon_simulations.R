#setwd("")
require(minqa)
require(inline)
require(caTools)
require(MASS)
source("functions.R")
as.real <- as.numeric

# NOTE: MUST first source the following files from https://github.com/fpgpr/DemoScripts (also linked via Dryad)
# source("ICA_functions_SHOWCASE.R")
# source("cxxFunctions_OU_Generic_SHOWCASE.cxx")

# The following functions (reorder.by.max() through Reconstructed_Curve()
#    are adapted from Hadjipantelis et al. 2013)

reorder.by.max <- function(M){
  #returns M with columns rearranged by absolute max function values
  M[,order(sapply(1:ncol(M), function(x) which(abs(M[,x])==max(abs(M[,x])))))]
}

flip.ICs <- function(ICs, flip.vectors = TRUE, zero.intercept = TRUE){
  #apply flip.vector to all columns in an IC matrix 
  ICs <- t(ICs)
  flipped <-sapply(1:ncol(ICs), 
                   function(x) flip.vector(ICs[,x],flip=flip.vectors, zero=zero.intercept))
  t(flipped)
}

flip.vector <- function(x, flip = TRUE, zero = TRUE){
  #flips a vector if the maximum value is negative
  # aligns lhs to zero
  if(flip) if(max(x) != max(abs(x))) x <- -x
  if(zero) x <- x + (0-x[1])
  x
}

GetThetas <-function(Distances, Data){
  #This is not a generic function, it assumes that Data is 128-by-1 and Distances 128-by-128
  Ksc = Distances;
  Q1  = Data;  
  
  N <- ncol(Distances)
  
  ff<-c();
  for (u in 1:100){
    SubSample <- sort(sample.int(n=N,size=100));  Ksc100 <- Ksc[SubSample,SubSample];  Q1_100 <- Q1[SubSample];
    S0 <-  uobyqa( log(runif(3)) ,fLogLik_General_Only_F,  X=Ksc100,Y=Q1_100)  
    ff<- c(ff, exp( S0$par), S0$fval)
  }
  GG <- matrix(ff, nrow=4)
  #Sanity Check for length scales
  San <- which(GG[2,] > max(Ksc));
  GG[2,San] <- max(Ksc)
  San <- which(GG[1,]  <GG[3,]);
  GG[2,San] <- mean(GG[2,-San])
  
  ThetasMean= rowMeans(GG)[1:3] 
  ThetasMedian= c( median(GG[1,]), median(GG[2,]), median(GG[3,]) )  
  S<-uobyqa( log(runif(3)) ,fLogLik_General_Only_F,  X=Ksc , Y=Q1);
  ThetasStraight <- exp(S$par); 
  
  return( list(ThetasMean, ThetasStraight, ThetasMedian))
}

Full_Curve_Estimate <- function(Y , X, thetas, new_X , Basis , Sample_Mean){
  if(missing(Sample_Mean)) Sample_Mean <- rep(0,nrow(Basis))
  N <- nrow(Y)
  P1 = Predictions_ForAspecificNode (Y= Y[1:N,1], X= X, Theta= thetas[1,], new_X= new_X)
  #print( Y[1:128,1]  )
  P2 = Predictions_ForAspecificNode (Y= Y[1:N,2], X= X, Theta= thetas[2,], new_X= new_X) 
  P3 = Predictions_ForAspecificNode (Y= Y[1:N,3], X= X, Theta= thetas[3,], new_X= new_X) 
  Ms= c( P1$Means, P2$Means, P3$Means );  Vs= c( P1$Vars, P2$Vars, P3$Vars );
  
  #print(Vs)
  Curve_Estimate= Reconstructed_Curve( Means= Ms, Vars= Vs, Basis =  t(Basis), Sample_Mean ) 
  
  return (Curve_Estimate)
}

Predictions_ForAspecificNode <- function(Y,X,Theta,new_X){  
  #Theta: OU hyperparameters
  #X : Pairwise distances of known points in the phylogeny
  #Y : Values at the tips    
  #new_X : Distance between the point of estimation and the known tips.
  # returns two-element list. 1st element : Means, 2nd element Variances
  N <- length(Y);#Number of tips   
  M <- dim(X)[1]
  #new_X : Pairwise distances of new point in the phylogeny with the respect to the existing ones. 
  
  K  <- matrix( rep(0, N^2) , ncol= N)#covariance matrix
  
  s_f = (  Theta[1]  ) #function's amplitude
  l =   (  Theta[2]  ) #characteristic lengh scale
  s_n = (  Theta[3]  ) #noise's amplitude 
  
  for ( i in 1:N){
    for ( j in i:N){ 
      K[i,j] = s_f^2 * exp(-(abs( X[i,j] ))/(l))  #+ s_c
    }
  } 
  K = K + t(K)
  K = K + diag(N) * (s_n^2 - s_f^2)# -s_c
  
  #Calculate the COVAR from each data point to the new one:
  K_x_xs = rep( 0,N)
  K_xs_xs = s_f^2 + s_n^2 # + s_c
  #Set Means and Vars
  Means =  0
  Vars  =  0
  
  xs_coord = new_X
  for (i in 1:N){
    K_x_xs[i] = s_f^2  * exp(  -(abs(  xs_coord[i]))/(l)  )# + s_c
  }
  Means =  K_x_xs %*% solve(K  ,Y)
  Vars  =  K_xs_xs -  K_x_xs %*% solve( K , K_x_xs)  
  return(list(Means=Means,Vars=Vars))
}
#================================================================================================

#Ancestral Node Reconstruction
Reconstructed_Curve <-function( Means= M, Vars= V, Basis = ICs, Sample_Mean = Sample_Mean ){
  #Mean Curve
  MeanCurve = t(Means) %*% Basis;
  #Standard Deviation
  Z = sqrt(rowSums(t(Basis^2) %*% diag(Vars) ))
  #3 column matrix with lower, mean and upper CI for prediction
  PredMat =  matrix(c(MeanCurve + Sample_Mean ,MeanCurve- (1.96*Z) + Sample_Mean ,  MeanCurve + (1.96*Z) + Sample_Mean ), ncol=3)
  #x11(); matplot(PredMat,type='l' ) ;
  return(PredMat)
}

pgpr <- function(tree,tip_coefficients,x_eval)
{
  XPLANE <- seq(.001,.999,length=100)
  SignalsAtTheTips <- apply(tip_coefficients,1,FUN = function(X) logit_inv(X,XPLANE))
  SignalsAtTheTips <- SignalsAtTheTips[,paste("t",1:N,sep="")]
  Ksc <- cophenetic(tree)
  thetas <- Nthetas <- matrix(0,3,3)
  MeanCurve   =  rowMeans(SignalsAtTheTips)
  PCA <- eigen(cov(t(SignalsAtTheTips)))
  ica <- cubica34(t(PCA$vectors[,1:3]))
  ica$y <- t(reorder.by.max(t(ica$y)))
  ica$Evals <- PCA$values
  ica$Evectors <- PCA$vectors
  ica$y <- flip.ICs(ica$y)
  ICA_MM <- t(ginv(ica$y)) %*% as.matrix(SignalsAtTheTips)
  ICA_sig <- ica$y
  
  Nthetas[1,] <- GetThetas(Data=ICA_MM[1,1:N], Ksc)[[1]]
  if (Nthetas[1,1] < Nthetas[1,3]) { Nthetas[1,] = c(0.00001,1, sd( ICA_MM[1,1:N])) }
  Nthetas[2,] <- GetThetas(Data=ICA_MM[2,1:N], Ksc)[[1]]
  if (Nthetas[2,1] < Nthetas[2,3]) { Nthetas[2,] = c(0.00001,1, sd( ICA_MM[2,1:N])) }
  Nthetas[3,] <- GetThetas(Data=ICA_MM[3,1:N], Ksc)[[1]]
  if (Nthetas[3,1] < Nthetas[3,3]) { Nthetas[3,] = c(0.00001,1, sd( ICA_MM[3,1:N])) }
  
  Estimated_thetas=Nthetas
  E_thetas =  Estimated_thetas; 
  E_thetas[,3]= 0.00001;
  
  NodeEstimate <- matrix(0,tree$Nnode,length(XPLANE))
  y_eval <- matrix(0,tree$Nnode,length(x_eval))
  coefs <- matrix(0,tree$Nnode,2)
  i <- 1
    node <- i
    NodeEstimate[i,] <- Full_Curve_Estimate(Y=t(ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[node+N,1:N], Basis = t(ICA_sig))[,1]
  ret_coefs <- coef(glm(XPLANE~NodeEstimate[1,],family=quasibinomial("logit")))
  ret <- logit_fx(ret_coefs,x_eval)
  if(ret[1]>.05 | ret[length(x_eval)]<.95 | abs((-ret_coefs[1]/ret_coefs[2])-5)>2) # checks for convergence failure
  {
    obj <- 1
    class(obj) <- "try-error"
    return(obj)
  }
  return(ret)
}

N <- 128
niters <- 1000
x_eval <- seq(0,15,length=100)
actual_y <- inv_curve_y <- pgpr_diff <- inv_curve_diff <- anc_y <- y_diff <- pgpr_y <- matrix(NA,nrow = niters,length(x_eval))
min_inv <- min_y <- min_pgpr <- rep(Inf,length(x_eval))
max_inv <- max_y <- max_pgpr <- rep(-Inf,length(x_eval))
inv_EC50 <- y_EC50 <- pgpr_EC50 <- vector("double",niters)
j <- 0
counter <- 0
while(counter<niters)
{
  j <- j+1
  print(j)
  data <- try(f(seed = j,N = 128),silent=TRUE)
  if(class(data)=="try-error") next
  pgpr_results <- try(pgpr(data$tree,data$wilting$true$tip_coefficients,x_eval=x_eval),silent=TRUE)
  if(class(pgpr_results)=="try-error") next
  i=1
    node <- i
  pgpr_y[counter+1,] <- pgpr_results
  actual_y[counter+1,] <- logit_fx(data$wilting$true$anc_coefficients[node,],x = x_eval)
  inv_curve_y[counter+1,] <- logit_fx(data$wilting$observed$anc_coefficients[node,],x = x_eval)
  anc_y[counter+1,] <- pgls_curve_wrong(data$tree,data$wilting$true$tip_coefficients,x_eval = x_eval)[1,]
  pgpr_diff[counter+1,] <- pgpr_y[counter+1,]-actual_y[counter+1,]
  inv_curve_diff[counter+1,] <- inv_curve_y[counter+1,]-actual_y[counter+1,]
  y_diff[counter+1,] <- anc_y[counter+1,] - actual_y[counter+1,]
  
  inv_curve_individual <- apply(inv_curve_y,2,function(X) mean(X,na.rm=TRUE))
  y_curve_individual <- apply(anc_y,2,function(X) mean(X,na.rm=TRUE))
  pgpr_y_individual <- apply(pgpr_y,2,function(X) mean(X,na.rm=TRUE))
  
  inv_curve <- apply(inv_curve_diff,2,function(X) mean(X,na.rm=TRUE))
  y_curve <- apply(y_diff,2,function(X) mean(X,na.rm=TRUE))
  pgpr_curve <- apply(pgpr_diff,2,function(X) mean(X,na.rm=TRUE))

  counter <- counter+1
  ret <- data.frame(x=x_eval,inv_curve_root=inv_curve,y_root=y_curve,pgpr_root=pgpr_curve,
                    inv_curve_individual=inv_curve_individual,y_curve_individual=y_curve_individual,pgpr_y_individual=pgpr_y_individual
                    )
  write.table(ret,"anc_curve.csv",sep=",",row.names = FALSE)
}
