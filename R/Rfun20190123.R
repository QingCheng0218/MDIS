propRfun = function(x, Y, Y.res){
  RhoEst <- mean(x[,1]*x[,2]);
  
  # Y <- Y-mean(Y);
  
  #------------------------------------------------#
  #Y.res=1:E(Y-Xbeta)(X1-rhoX2)(X2-rhoX1);residual-based screening procedure
  #else :  EY(X1-rhoX2)(X2-rhoX1); (Y center)The centralized response-based screening procedure
  #------------------------------------------------#
  # linear regression should contain intercept term.
  res.lm <- lm(Y~x[,1]+x[,2]+I(x[,1]*x[,2]));
  beta1 <- as.numeric(coef(res.lm)[2]);
  beta2 <- as.numeric(coef(res.lm)[3]);
  
  RES1 <- as.vector(res.lm$residuals);
  RESX <- cbind(x[,1] - RhoEst*x[,2], x[,2] - RhoEst*x[,1]);
  sigma1 <- sd(RESX[,1]*RESX[,2])*sd(RES1);
  
  #sigma2 <- (1-RhoEst^2)*sqrt(3*(beta1^2 + beta2^2 + 2*beta1*beta2*RhoEst*(2*RhoEst^2-1)));
  sigma2 <- sd(beta1*x[,1]*RESX[,1]*RESX[,2]+beta2*x[,2]*RESX[,1]*RESX[,2]);
  
  if (Y.res==1) {
    # linear regression should contain intercept term.
    RES <- as.vector(lm(Y~ x[,1]+x[,2])$residuals);
    P.ASYMSTD <- sigma1/sqrt(n);
  } else {
    RES <- as.vector(Y);
    P.ASYMSTD <- sqrt(sigma1^2 + sigma2^2)/sqrt(n);
  }
  
  
  P.STATIC <- mean(RESX[,1]*RESX[,2]*RES);
  
  
  P.zscore <- P.STATIC/P.ASYMSTD;
  
  return(P.zscore)
}
proploopR = function(X0, Y, Y.res, at){
  Zscore <- rep(0,nrow(at));
  #z.ml <- TPR <- FDR.T <- rep(0,length(L))
  #P.ZSCORE <- matrix(rep(0, nrow(at)*iterMax), ncol=iterMax)
  #no.rej <- rep(0, iterMax);
  #------------------------------------------------#
  # X0 <- apply(X0,2,scale);# standardize variable X
  Y <- Y - mean(Y)
  
  for ( j in 1:nrow(at))
  {  		 
    x <- X0[,at[j,]];
    Zscore[j] <- abs(propRfun(x, Y, Y.res));
    
    #if (abs(P.zscore[j]) >= abs(qnorm((1-pVALUE/2)))) temp <- append(temp,j);
  }
  
  Trank <- order(abs(Zscore), decreasing=T);
  return(list(Zscore = Zscore, Trank = Trank));
}

# Obtain the parameters for generate simulation data.
parafun <- function(p, r1, d1, rho, b, kappa, model, Model){
  #model(Sigma): 1:Block-diagonal(d1);2:rho^|i-j|; 3: identity matrix  4:block-diagonal(autoregressive); 5: banding matrix(d1)
  #Model(int_B): 1:upper  triangle with r1 row(col);   2:upper  triangle & rho^|i-j|; 
  #       3:"|-" matrix (row:1-r1;col:1-r1)  4:off-block-diagonal(row:r1*d1;col:r1*d1;
  #       cut-off bockdiagonal with d1(e.g:d1=10,r1=3) )
  if(model==1){
    R1 = rho^(matrix(rep(1,d1^2), d1, d1) - diag(d1));
    L = p/d1;
    S = kronecker(diag(L), R1);
  }else if(model==2){
    S <- rho^(abs(matrix(rep(c(1:p), each=p), ncol = p, byrow=T) - matrix(rep(c(1:p), each = p), ncol=p)));
  }else if(model==3){
    S <- diag(p);
  }else if(model==4){
    R1 = rho^(abs(matrix(rep(c(1:d1), each = d1), ncol = d1, byrow =  T) - matrix(rep(c(1:d1), each = d1), ncol = d1)));
    L = p/d1;
    S = kronecker(diag(L), R1);
  }else if(model==5){
    S0 <- rho^(abs(matrix(rep(c(1:p), each=p), ncol = p, byrow=T) - matrix(rep(c(1:p), each = p), ncol=p)));
    S <- as.matrix(band(S0, -d1, d1));
  }
  
  
  
  if(Model==1){
    B1 <- (kappa^(matrix(rep(1, r1^2), r1, r1) - diag(r1)) - diag(r1))/2; 
    B <- as.matrix(bdiag(B1,matrix(0, p - r1, p - r1)));
    tint <- combinations(r1, 2); 
  }else if(Model==2){
    B1 <- (kappa^(abs(matrix(rep(c(1:r1),each=r1),ncol=r1,byrow=T)
                      -matrix(rep(c(1:r1),each=r1),ncol=r1)))-diag(r1))/2;
    B <- as.matrix(bdiag(B1,matrix(0,p-r1,p-r1)));
    tint <- combinations(r1,2);
  }else if(Model==3){
    B<-matrix(rep(0,p*p),ncol=p);
    B[1:d1, 1:r1] <-kappa/2;
    B[1:r1, 1:d1] <- kappa/2;
    diag(B) <- 0;
    id.int <- which(B!=0, arr.ind=T);
    tint <- id.int[order(id.int[, 1]), ];
    tint <- tint[tint[, 1] < tint[, 2], ];
  }else if(Model==4){
    B0 <- (kappa^(matrix(rep(1, d1^2), d1, d1) - diag(d1)) - diag(d1))/2; 
    B1 <- B0;
    for(i in 1:(r1-1))
    {
      B1 <- as.matrix(bdiag(B0, B1));
      i <- i+1;
    }
    
    B2 <- matrix(rep(kappa, d1^2*r1^2), d1*r1)/2;
    B3 <- B2 - B1;
    B <- as.matrix(bdiag(B3, diag(p - r1*d1)));
    diag(B) <- 0;
    id.int <- which(B != 0, arr.ind = T);
    tint <- id.int[order(id.int[, 1]), ];
    tint <- tint[tint[, 1] < tint[, 2], ];
  }
  
  b.sit <- as.vector(tint)[-which(duplicated(as.vector(tint)))];
  beta <- rep(0,p);
  beta[b.sit] <- b;   
  
  return(list(S = S, B = B, tint = tint, beta=beta, b.sit = b.sit));
  
}


# JCIS method from 'Strong Sure Screening of Ultra-high Dimensional Categorical Data'
JCISrank <- function(x1, x2, y){
  n = length(y);
  term1 = abs(sum((x1 - mean(x1))*(x2 - mean(x2))*(y - mean(y))));
  term2 = sum((x1 - mean(x1))^2);
  term3 = sum((x2 - mean(x2))^2);
  term4 = sum((y - mean(y)));
  res = sqrt(n)*term1/sqrt(term2*term3*term4);
  return(res);
}

JCISfun <- function(X, y, at){
  p1 = dim(at)[1];
  n = dim(X)[1];
  rankres = rep(0, p1);
  for(j in 1:p1){
    d1 = at[j, 1];
    d2 = at[j, 2];
    x1 = X[, d1];
    x2 = X[, d2];
    rankres[j] = JCISrank(x1, x2, y);
  }
}

Tibshirani <- function(Y, X0, at){
  groups <- unique(Y);
  if(length(groups) > 2){
    y <- as.vector(1*(Y > median(Y)));
    groups <- unique(y);
  }else{
    y <- Y;
  }
  
  ind1 <- which(y == groups[1])
  ind2 <- which(y == groups[2])
  
  X1 <- X0[ind1,]
  X2 <- X0[ind2,]
  
  # mean center and scale X in each class.
  X1 <- t((t(X1) - apply(X1, 2, mean))/apply(X1,2,sd))
  X2 <- t((t(X2) - apply(X2, 2, mean))/apply(X2,2,sd))
  
  cors1 <- cor(X1, X1, method = "pearson")
  cors2 <- cor(X2, X2, method = "pearson")
  
  # the statistic
  # stats.M <- atanh(cors1) - atanh(cors2)
  stats.M <- cors1 - cors2;
  Zscore = abs(stats.M[at]);
  rank.i = order(Zscore, decreasing = TRUE);
  return(rank.i)
}

#~~~~~~~~~~~~~~~~~~#
#Modified FOR screening(Aug.29,2017)
#~~~~~~~~~~~~~~~~~~#
#-------------------------------Criterias for proposal---------------------------------#

#-------------------------------------------------------------------------#
# Functions to compute the criteria in the simulation.
# 1. M: to compute the minimum model size to ensure the inclusion of all active predictors. 
# 2. mqtl: to compute the 5%, 25%, 50%, 75% and 95% quantiles of the minimum model size out of 1,000 replications.
# 3. Sel.rate: to compute the proportion that every single active predictor is selected 
#    for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.


M <- function(true.v,rank.mtx) {
  # Input
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  r <- min(dim(rank.mtx)[2],length(rank.mtx))
  M <- c()
  for (j in 1:r) {
  index <- match(true.v,rank.mtx[,j]);
  id <- which(is.na(index));
  if(length(id)!=0){
  index[id] <- dim(rank.mtx)[1];}
  M[j] <- max(index);
  }
  return(M)
}

mqtl <- function(M) {
  # Input    
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  # Output
  # 5%,25%,50%,75%,95% quantiles of minimum model sizes out of 1000 replications
  quantile(M, probs =c(0.05,0.25,0.5,0.75,0.95))
}

# ifc1=0，the model size is: c1 times of true predictor number
# if c1!=0 the model size is: c1*floor(n/log(n))

Sel.rate <- function(n, c1, c2, true.v, rank.mtx) {
  # Input
  # n        :  the sample size
  # c        :  coeficient of cutoffs, for example c=2, cutoff=2[n/log(n)]
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  #
  #all.rate :   the  proportions that all active predictors are selected 
  
  if(c1==0){d <- c2*length(true.v)}else{d <- c1*floor(n/log(n))}
  rank.mtx.sel <- rank.mtx[1:d,]
  r <- min(dim(rank.mtx)[2],length(rank.mtx))
  p0 <- length(true.v)
  R <- matrix(0,p0,r)
  rate<-c()
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i]))==0) }
    rate[i]<-mean(R[i,])
  }
  all.rate <- mean(apply(R, 2, sum)==length(true.v))
  return(list(rate = rate, all.rate = all.rate))
}

# result for screening criteria.
resfun <- function(true.v, c1, c2id, rank.mtx){
  modelsize <- M(true.v,rank.mtx);
  modelqt <- mqtl(modelsize);
  
  selrate1 <- matrix(0, length(true.v), length(c2id)) ;
  selrate2 <- rep(0, length(c2id));
  for(i in 1:length(c2id)){
    c2 <- c2id[i];
    selr <- Sel.rate(n, c1, c2, true.v, rank.mtx)
    selrate1[, i] <- selr$rate; 
    selrate2[i] <- selr$all.rate
  }
  selres <- rbind(colMeans(selrate1), selrate2);
  
  return(list(modelsize = modelsize, modelqt = modelqt, selrate1 = selrate1, selrate2 = selrate2, selres = selres));
}

ispcresfun1 <- function(rank.t, rank.ispc, c2id){
  # result for big size of true interactions
  nr = length(rank.t);
  prres1 = matrix(0, nrow = nr, ncol = iterMax);
  prres2 = matrix(0, nrow = nr, ncol = iterMax);
  prres3 = matrix(0, nrow = nr, ncol = iterMax);
  
  for(iter in 1:iterMax){
    prres1[, iter] = rank.t%in%rank.ispc[1: (c2id[1]*nr), iter];
    prres2[, iter] = rank.t%in%rank.ispc[1: (c2id[2]*nr), iter];
    prres3[, iter] = rank.t%in%rank.ispc[1: (c2id[3]*nr), iter];
  }
  res = rbind(c(mean(colMeans(prres1)), mean(colSums(prres1)==nr)),
              c(mean(colMeans(prres2)), mean(colSums(prres2)==nr)),
              c(mean(colMeans(prres3)), mean(colSums(prres3)==nr)));
  return(res);
}


ispcresfun2 <- function(rank.t, rank.ispc, c1, c2id){
  # result for small size of  true interactions
  iterMax = dim(rank.ispc)[2];
  nr = length(rank.t);
  n3 = length(c2id);
  n4 = floor(n/log(n));
  prres = array(0, dim =  c(nr, iterMax, n3));
  for(j in 1:n3){
    for(iter in 1:iterMax){
      if(c1==0){
        prres[, iter, j] = rank.t%in%rank.ispc[1: (c2id[j]*nr), iter];
      }else{
        prres[, iter, j] = rank.t%in%rank.ispc[1: (c2id[j]*n4), iter];
      }

    }
  }
  
  res = rbind(apply(prres, 3, rowMeans),
              colMeans(apply(prres, 3, colSums)==nr));
  return(res);
}


