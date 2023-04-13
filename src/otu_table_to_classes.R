
otu_table_to_classes = function(OTUs,distance=c(),euclidean=FALSE,Nclus=15){
  # This function takes an OTU table (Samples in columns) and an optional distance matrix describing
  # the beta-diversity between samples and clusters communities according to this
  # distance with partition-around-medoids. It returns a vector quantifying the 
  # Calinski-Harabasz index for each clustering threshold and the distance matrix.
  # If the distance matrix was not provided, it also computes and returns the Jensen-Shannon
  # divergence between samples.
  # dependencies: dist.JSD.R, pam.clustering.R
  library(phyloseq)
  library(clusterSim) # index G1
  library(ade4)
  # library(clusterCrit) # intCriteria
  # library(fpc) # cluster.stats
  
  
  # --- Samples should be in columns
  OTUs.relAb=apply(OTUs,2,function(x)(x/sum(x))) # Normalize the table
  OTU_table=otu_table(as.matrix(OTUs),taxa_are_rows = TRUE)
  #browser()
  if(length(distance)==0){
    mes="... Computing JSD matrix. "
    cat(mes)
    data.dist = dist.JSD(OTUs.relAb)
    euclidean=TRUE
  }else{
    data.dist = distance
    if(is.euclid(data.dist)){ # I comment it, let the user verify it
      #euclidean = TRUE
    }
  }
  if(Nclus > dim(OTUs.relAb)[2]){ # we can't have more clusters than elements
    Nclus= dim(OTUs.relAb)[2]-1
  }

  nclusters=NULL
  #cluster.stats=list()
  for (k in 1:Nclus) { 
    mes=paste("... Computing PAM clustering for",k,"clusters. \n")
    cat(mes)
    if(k==1){
      nclusters[k]=NA
      #cluster.stats[[k]]=NA
    }else{
      data.cluster_temp=pam.clustering(data.dist, k)
      if(euclidean == FALSE){ 
        nclusters[k]=index.G1(t(OTUs.relAb),data.cluster_temp) 
        #nclusters[k]=indexG1.check(t(OTUs.relAb),data.cluster_temp) 
        #intCriteria(x,data.cluster_temp,"Calinski_Harabasz")
      }else{
        nclusters[k]=index.G1(t(OTUs.relAb),data.cluster_temp, d = data.dist,
                              centrotypes = "medoids") # <<- 
        #nclusters[k]=indexG1.check(t(OTUs.relAb),data.cluster_temp, d = data.dist,
        #                      centrotypes = "medoids") 
      } 
      # cluster.stats[[k]]=cluster.stats(data.dist, data.cluster_temp) # more possible statisticss
    }
  }
  if(length(distance)==0){
    return(list(data.dist=data.dist,nclusters=nclusters)) #,cluster.stats=cluster.stats))
  }else{
    return(list(nclusters=nclusters)) #,cluster.stats=cluster.stats))
  }
}