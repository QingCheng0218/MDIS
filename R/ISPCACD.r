
library(pcaPP)   # for cor.fk (x, y = NULL)  #fast algorithms to calculate Kendall's tau

 
#######################################################################
# partial correlation procedure   ACelerateD version
##########################################################################
ISPCACD = function(x,y,modelsize=NULL,standardization=TRUE, method = "pearson", bigN = 100000){#, fullreturn=FALSE){ # do not need full return
#standardize data if needed
   if (standardization) {x = scale(x); y=scale(y)}
   n=length(y);p=dim(x)[2]

#set selected model size if NULL
   if (is.null(modelsize)) {modelsize =  max(floor(n/log(n)),2)}

   partialcorr = NULL
   threshold = 0
   coordinate = NULL

   if (method=="kendall") {MainCor = abs(apply(x,2,cor.fk,y))} else {MainCor = abs(cor(x,y))}
   # rearrange x by marginal strength
   MCorder = order(MainCor,decreasing=T)
   x = x[,MCorder]; NewMainCor=MainCor[MCorder]

   for (i in 1:(p-1)){
      if (NewMainCor[i]<threshold/(1+2*threshold)) {
         if (method=="kendall"){ttcorr = apply(x[,i]*x[,c(i:p)],2,cor.fk,y)} else {ttcorr = cor(x[,i]*x[,c(i:p)],y,method=method)}
         if (max(abs(ttcorr))<=threshold/(1+2*threshold)) next
         index = which(abs(ttcorr)>threshold/(1+2*threshold)) +i-1     # find index of interactions with large marginal correlation

         if (length(index)>1)
            {tempcorr = apply(x[,index],2,mypcor,y,x[,i],method=method)}         # calculate partial correlation
         else  {tempcorr = mypcor(x[,index],y,x[,i],method=method)}

 	       if ( max(abs(tempcorr))<=threshold) next                           # skip when correlations are small
         partialcorr= c(partialcorr,tempcorr)                                    # otherwise, store it

         coordinate = rbind(coordinate,cbind(rep(i,length(index)),index))           # and corresponding positions

      } else {

         tempcorr = apply(x[,c(i:p)],2,mypcor,y,x[,i],method=method)         # calculate partial correlation
      	 #tempcorr = apply(x[,i]*x[,c(i:p)],2,cor.fk,y)                      # calculate correlation
 	       if ( max(abs(tempcorr))<=threshold) next                           # skip when correlations are small
         partialcorr= c(partialcorr,tempcorr)                                    # otherwise, store it
         coordinate = rbind(coordinate,cbind(rep(i,p+1-i),c(i:p)))           # and corresponding positions
      }
      if (length(partialcorr)>bigN) {                       #when partialcorr is long enough, use sparse vector to record it and save space
           ob= ThresholdVector(partialcorr,Size=modelsize,threshold=threshold)
           partialcorr = ob$v; position = ob$p; coordinate = coordinate[position,]
           threshold = min(abs(partialcorr))
        }
    }
    partialcorr= c(partialcorr,mypcor(x[,p],y,x[,p],method=method))              # final part, as we can not use apply to a single vector
    coordinate = rbind(coordinate,c(p,p))           # and corresponding positions
    ob= ThresholdVector(partialcorr,Size=modelsize,threshold=threshold)
    partialcorr = ob$v; position = ob$p; coordinate = coordinate[position,]
    # recall that we re-ordered features, now it is time to put it back
    coordinate = cbind(MCorder[coordinate[,1]], MCorder[coordinate[,2]])
    for (k in 1:dim(coordinate)[1]) {
       if (coordinate[k,1]>coordinate[k,2])
         {temp=coordinate[k,1]; coordinate[k,1]=coordinate[k,2];coordinate[k,2]=temp}      #make the index in correct order
    }
    return(list(interactionindex=coordinate,topcorr=partialcorr))
}
 
ISPCACD_adjust = function(x, y, MAT.idx, modelsize=NULL,standardization=TRUE, method = "pearson", bigN = 100000){#, fullreturn=FALSE){ # do not need full return
  #standardize data if needed
  if (standardization) {x = scale(x); y=scale(y)}
  n=length(y);p=dim(x)[2]
  
  #set selected model size if NULL
  # if (is.null(modelsize)) {modelsize =  max(floor(n/log(n)),2)}
  if (is.null(modelsize)) {modelsize =  p*(p+1)/2}
  
  partialcorr = NULL
  threshold = 0
  coordinate = NULL
  
  if (method=="kendall") {MainCor = abs(apply(x,2,cor.fk,y))} else {MainCor = abs(cor(x,y))}
  # rearrange x by marginal strength
  MCorder = order(MainCor,decreasing=T)
  x = x[,MCorder]; NewMainCor=MainCor[MCorder]
  
  for (i in 1:(p-1)){
    if (NewMainCor[i]<threshold/(1+2*threshold)) {
      if (method=="kendall"){ttcorr = apply(x[,i]*x[,c(i:p)],2,cor.fk,y)} else {ttcorr = cor(x[,i]*x[,c(i:p)],y,method=method)}
      if (max(abs(ttcorr))<=threshold/(1+2*threshold)) next
      index = which(abs(ttcorr)>threshold/(1+2*threshold)) +i-1     # find index of interactions with large marginal correlation
      
      if (length(index)>1)
      {tempcorr = apply(x[,index],2,mypcor,y,x[,i],method=method)}         # calculate partial correlation
      else  {tempcorr = mypcor(x[,index],y,x[,i],method=method)}
      
      if ( max(abs(tempcorr))<=threshold) next                           # skip when correlations are small
      partialcorr= c(partialcorr,tempcorr)                                    # otherwise, store it
      
      coordinate = rbind(coordinate,cbind(rep(i,length(index)),index))           # and corresponding positions
      
    } else {
      
      tempcorr = apply(x[,c(i:p)],2,mypcor,y,x[,i],method=method)         # calculate partial correlation
      #tempcorr = apply(x[,i]*x[,c(i:p)],2,cor.fk,y)                      # calculate correlation
      if ( max(abs(tempcorr))<=threshold) next                           # skip when correlations are small
      partialcorr= c(partialcorr,tempcorr)                                    # otherwise, store it
      coordinate = rbind(coordinate,cbind(rep(i,p+1-i),c(i:p)))           # and corresponding positions
    }
    if (length(partialcorr)>bigN) {                       #when partialcorr is long enough, use sparse vector to record it and save space
      ob= ThresholdVector(partialcorr,Size=modelsize,threshold=threshold)
      partialcorr = ob$v; position = ob$p; coordinate = coordinate[position,]
      threshold = min(abs(partialcorr))
    }
  }
  partialcorr= c(partialcorr,mypcor(x[,p],y,x[,p],method=method))              # final part, as we can not use apply to a single vector
  coordinate = rbind(coordinate,c(p,p))           # and corresponding positions
  ob= ThresholdVector(partialcorr,Size=modelsize,threshold=threshold)
  partialcorr = ob$v; position = ob$p; coordinate = coordinate[position,]
  # recall that we re-ordered features, now it is time to put it back
  coordinate = cbind(MCorder[coordinate[,1]], MCorder[coordinate[,2]])
  for (k in 1:dim(coordinate)[1]) {
    if (coordinate[k,1]>coordinate[k,2])
    {temp=coordinate[k,1]; coordinate[k,1]=coordinate[k,2];coordinate[k,2]=temp}      #make the index in correct order
  }
  
  # --------------------------------------------------
  # remove Xj^2
  idx <- which(coordinate[, 1]==coordinate[, 2]);
  coordinate <- coordinate[-idx, ];
  rank.i = MAT.idx[coordinate];
  # --------------------------------------------------
  
  return(list(rank.i = rank.i, interactionindex=coordinate,topcorr=partialcorr))
}


#######################################################################
# some support functions
##########################################################################
ThresholdVector = function(SV,Size=NULL,threshold=0){
   if (is.null(SV)) return(list(v=NULL,p=NULL))
   index=c(1:length(SV))
   if (threshold>0) {index = which(abs(SV)>=threshold); SV=SV[index]}         # threshold it ">="!
   if (length(SV)==0) return(list(v=NULL,p=NULL))                            # nothing left, return NULL
   if (is.null(Size)) return(list(v=SV,p=index))
   if (Size>=length(SV)) return(list(v=SV,p=index))

      fullorder=order(abs(SV),decreasing=T)
      position = fullorder[1:Size]
      value = SV[position]; position =index[position]
      return(list(v=value,p=position))
}

pcor4 = function(A, method="pearson"){    # given a n by 4 matrix, calculate the partial correlation of first two conditional on last two.
      if (method=="kendall") {Corr=cor.fk(A)}
      else {  Corr = cor(A,method=method)}
      #Inv  = solve(Corr); return(-Inv[1,2]/sqrt(Inv[1,1]*Inv[2,2]))  #this is slightly faster but may cause a problem when Corr is near singular
      if (Corr[3,4]>=1-0.000001) Corr=Corr[-4,-4]    # make it work for quadratic effects
      return(det(Corr[-1,-2])/sqrt(det(Corr[-1,-1]%*%Corr[-2,-2])))
}



mypcor = function(xj,y,xi,method="pearson"){ # calculate partial correlation between y and xi*xj conditional xi and xj
   A = cbind(y,xi*xj,xi,xj)
   return(pcor4(A,method=method))
}