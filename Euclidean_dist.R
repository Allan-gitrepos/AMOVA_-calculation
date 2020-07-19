#Function for calculating euclidean distances from DArTSeq 2row or 1row SNP data
#Author   :Allan. V
#Date     : 18 - July - 2020 
#Email    :albertoogy@gmail.com
#ORCID iD :https://orcid.org/0000-0001-7230-081X

#Note: This code may be slower than the native "dist" function 
#      for smaller data. But, for bigger data it provides a stable
#      real time output without exhausting memory.

#      Use-case: Efficient when used for more than 500 individuals or samples



euc_dis=function(data=NULL, save=T, type=NULL,callrate_ind=0.50, callrate_loc=0.80 ,repavg=0.95, mono.rm=F){
  
  cat("Loading the packages\n")
  if (!require("dartR")) install.packages("dartR")
  if (!require("dplyr")) install.packages("dplyr")
  #Reading the data
  if(is.null(data)){
    stop("Please provide a DArTSeq two row or one row SNP file")
  }
  gl=gl.read.dart(paste(data))
  cat("Filtering for call rate locus with threshold          :", callrate_loc, "\n")
  cat("Filtering for call rate individuals  with threshold   :", callrate_ind, "\n")
  cat("Filtering for reproducability of locus with threshold :", repavg, "\n")
  
  gl <- gl.filter.callrate(gl, method="ind", threshold = callrate_ind,mono.rm = F)
  gl <- gl.filter.repavg(gl, t=repavg)
  gl<-gl.filter.callrate(gl, method="loc",t=callrate_loc, mono.rm = F)
  
  
  
  cat("Data Read...........!!!!!!\n")
  
  cat("Setting indivual names\n")
  
  pop(gl)<-indNames(gl)
  
  cat("Loading the function for calculation \n")
  e.dis=function(A,B){
    if(any(is.na(A) | is.na(B))){
      i <- is.na(A) | is.na(B)
      temp=(rbind(A[!i], B[!i])) * sqrt(length(A) / length(A[!i]))
      out=sqrt(sum(((temp[1,] -temp[2,]) ^ 2),na.rm = T))
    }else{out=sqrt(sum(((A -B) ^ 2),na.rm = T))}
    return(out)
  }
  dis_ecl=function(a,ind_full){
    #g_mat is the whole matrix 
    dis=matrix(NA, nrow = length(ind_full), ncol =length(ind_full))
    colnames(dis)<-rownames(dis)<-ind_full
    for (i in 1:length(ind_full)) {
      cat(paste(i,"calculating distance between", ind_full[i], "and others"),sep="\n")
      for(j in 1:length(ind_full)){
        dis[i,j]<-e.dis(a[ind_full[i],], a[ind_full[j],])
      }
    }
    return(dis)
  }
  
  e_gl= function(gl) {
    
    #Converting to frequencies
    a=as.matrix(gl)
    x=which(a==2)#Changed to 1 #Homozys
    y=which(a==1)#Chaged to 0.5 #Hetrozys
    a[x]<-1;rm(x)
    a[y]<-0.5;rm(y)
    a=a[order(rownames(a)),]
    #Calculating pairwise euclidean distances
    ind_full=rownames(a)
    dis=dis_ecl(a, ind_full=ind_full)
    return(dis)
  }
  
  cat("Starting the pairwise distance estimation for big data")
  system.time({
    x=e_gl(gl)
    if(save){
      cat("Writing the Distance matrix to directory")
      if(is.null(type)){
        save(x,file= "Distance_mat.Rdata")
        write.csv(x, "Distance_mat.csv")
      }else if(type=="Rdata"){
        save(x,file= "Distance_mat.Rdata")
      }else if(type=="csv"){
        write.csv(x, "Distance_mat.csv")
      }else if (type=="both"){
        save(x,file= "Distance_mat.Rdata")
        write.csv(x, "Distance_mat.csv")
      }
    }
  })
  x=as.data.frame(x)
  return(x)
}

