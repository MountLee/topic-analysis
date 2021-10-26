# This function performs the T-SCORE methods for pLSI model.
# Input:  
#   K: The number of topics;
#   K0: The number of vertices selected through greedy step, after the kmeans,
#       we usually set it to be ceiling(1.5*K);
#   m: The number of clusters in the kmeans step. In the paper we use L to denote this quantity,
#       we usually set it to be 10*K;
#   D: The word-document matrix. The ith column of D is the empirical distribution of 
#       words in the ith document;
#   Mquantile: The quantile of the diagonal of matrix M, which is used to upper threshold the 
#       diagonal entries of M. By default it is set to be 1, which degenerates to the case
#       without pre-normalization of D.
# Output:
#   A_hat: Estimation of word-topic matrix A;
#   R: Matrix of ratios of left singular vectors;
#   V: Rows of V are the vertices among the rows of R, found through VH procedures 
#   Pi: Each row is a convex combination vectors, that is used to construct the corresponding
#       row in of R from the learned vertices in V.

norm_score <- function(K, K0, m, D, Mquantile=1, scatterplot=FALSE, kmeans_start = 200){
  library('rARPACK')
  library('nnls')
  p <- dim(D)[1]
  n <- dim(D)[2]
  M <- rowMeans(D)
  M_trunk <- pmax(M,quantile(M,Mquantile))
  
  obj <- svds(sqrt(M_trunk^(-1))*D, K)
  Xi <- obj$u
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(Xi[,2:K],2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m,kmeans_start = kmeans_start)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,0)
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- sqrt(M_trunk)*Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V, Pi=Pi))
}

vertices_est <- function(R,K0,m,kmeans_start){
  library(quadprog)
  K <- dim(R)[2] + 1
  
  #Step 2a
  obj <- kmeans(R,m,iter.max=200,nstart = kmeans_start)
  theta <- as.matrix(obj$centers)
  theta_original <- theta
  # plot(R[,1],R[,2])
  # points(theta[,1], theta[,2], col=2,lwd=4)
  
  #Step 2b'
  inner <- theta%*%t(theta)
  distance <- diag(inner)%*%t(rep(1,m)) + rep(1,m)%*%t(diag(inner)) - 2*inner
  top2 <- which(distance==max(distance),arr.ind=TRUE)[1,]
  theta0 <- as.matrix(theta[top2,])
  theta <- as.matrix(theta[-top2,])
  
  if (K0 > 2){
    for (k0 in 3:K0){
      inner <- theta%*%t(theta)
      distance <- rep(1,k0-1)%*%t(diag(inner))-2*theta0%*%t(theta)
      ave_dist <- colMeans(distance)
      index <- which(ave_dist==max(ave_dist))[1]
      theta0 <- rbind(theta0, theta[index,])
      theta <- as.matrix(theta[-index,])
    }
    theta <- theta0
  }
  
  #Step 2b
  comb <- combn(1:K0, K)
  max_values <- rep(0, dim(comb)[2])
  for (i in 1:dim(comb)[2]){
    for (j in 1:K0){
      max_values[i] <- max(simplex_dist(as.matrix(theta[j,]), as.matrix(theta[comb[,i],])), max_values[i])
    }
  }
  
  min_index <- which(max_values == min(max_values))
  
  #plot(theta[,1],theta[,2])
  #points(theta[comb[,min_index],1],theta[comb[,min_index],2],col=2,pch=2)
  
  return(list(V=theta[comb[,min_index[1]],], theta=theta_original))
}

simplex_dist <- function(theta, V){
  VV <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))%*%V
  D <- VV%*%t(VV)
  d <- VV%*%(theta-V[dim(V)[1],])
  
  A <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))
  b0 <- c(rep(0,dim(V)[1]-1),-1)
  
  
  obj <- try(solve.QP(D, d, A, b0))
  if(class(obj) == "try-error"){
    obj <- list(value = Inf)
  }
  
  return(sum((theta-V[dim(V)[1],]) ^2)+ 2*obj$value)
}

