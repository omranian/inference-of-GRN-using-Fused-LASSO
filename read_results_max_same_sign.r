
read_results <- function()
{
  
  # initialization of the network structure with respect to the overall data set
  links <- data.frame(from=c(),to=c(),connected=c(),weight=c(),stringsAsFactors=F)
  
  # initialization of the network structure with respect to each single data set
  links_cold <- data.frame(from=c(),to=c(),connected=c(),weight=c(),stringsAsFactors=F)
  links_heat <- data.frame(from=c(),to=c(),connected=c(),weight=c(),stringsAsFactors=F)
  links_ox <- data.frame(from=c(),to=c(),connected=c(),weight=c(),stringsAsFactors=F)
  links_lac <- data.frame(from=c(),to=c(),connected=c(),weight=c(),stringsAsFactors=F)
  
  # set the directory to the folder in which the results are stored
  files <- list.files("Results/obj")
  g_orf_syn <- genes$ORF
  names(g_orf_syn) <- genes$synonym
  
  # provide the data for creating the histogram of differences between the regression coefficiences obtained from each single data set
  hist_mat <- matrix(0,ncol=4)
  
  # read each file and make the networks
  for (f in files)
  {
    g_names <- as.character(g_orf_syn[which(names(g_orf_syn)==gsub("_lqa_at_each_tp.obj","",f) )])
    #load(paste("~/Heraklit-Nooshin//home/Nooshin/Network inference method (with Jeanne)/Results//obj_at_each/",f,sep=""))
    load(paste("Results/obj/",f,sep=""))
    beta <- fit$beta.opt
    
    mat <- matrix(beta,ncol=4,nrow=length(beta)/4,byrow=F)
    tmp <- names(beta)[1:(length(beta)/4)]
    rownames(mat) <- as.character(g_orf_syn[tmp])
    
    hist_mat <- rbind(hist_mat,mat)
    
    ids<- which(apply(mat,1,function(r)any(r!=0))==T)
    ids2 <- which(apply(mat[ids,],1,function(r)any(r<=0) & any(r>=0))==T)
    ids <- setdiff(ids,ids2)
    
    if (length(ids)>1)   
    {
      links <- rbind(links,data.frame(from=rownames(mat)[ids],to=rep(g_names,length(ids)),connected=as.integer(rep(1,length(ids))),weight=as.numeric(apply(mat[ids,],1,function(r)r[which.max(abs(r))])),stringsAsFactors=F))
      links_heat <- rbind(links_heat,
                          data.frame(from=rownames(mat)[ids],to=rep(g_names,length(ids)),
                                     connected=as.integer(rep(1,length(ids))),
                                     weight=as.numeric(mat[ids,1]),stringsAsFactors=F))
      links_cold <- rbind(links_cold,
                          data.frame(from=rownames(mat)[ids],to=rep(g_names,length(ids)),
                                     connected=as.integer(rep(1,length(ids))),
                                     weight=as.numeric(mat[ids,2]),stringsAsFactors=F))
      links_ox <- rbind(links_ox,
                        data.frame(from=rownames(mat)[ids],to=rep(g_names,length(ids)),
                                   connected=as.integer(rep(1,length(ids))),
                                   weight=as.numeric(mat[ids,3]),stringsAsFactors=F))
      links_lac <- rbind(links_lac,
                         data.frame(from=rownames(mat)[ids],to=rep(g_names,length(ids)),
                                    connected=as.integer(rep(1,length(ids))),
                                    weight=as.numeric(mat[ids,4]),stringsAsFactors=F))
    }
    else
      if (length(ids) == 1)
      {
        links <- rbind(links,data.frame(from=rownames(mat)[ids],to=g_names,connected=1,weight=as.numeric(mat[ids,][which.max(abs(mat[ids,]))]),stringsAsFactors=F))
        links_heat <- rbind(links_heat,
                            data.frame(from=rownames(mat)[ids],to=g_names,
                                       connected=1,weight=as.numeric(mat[ids,1]),stringsAsFactors=F))
        links_cold <- rbind(links_cold,
                            data.frame(from=rownames(mat)[ids],to=g_names,
                                       connected=1,weight=as.numeric(mat[ids,2]),stringsAsFactors=F))
        links_ox <- rbind(links_ox,
                          data.frame(from=rownames(mat)[ids],to=g_names,
                                     connected=1,weight=as.numeric(mat[ids,3]),stringsAsFactors=F))
        links_lac <- rbind(links_lac,
                           data.frame(from=rownames(mat)[ids],to=g_names,
                                      connected=1,weight=as.numeric(mat[ids,4]),stringsAsFactors=F))
      }
    rownames(links) <- NULL
    rownames(links_heat) <- NULL
    rownames(links_cold) <- NULL
    rownames(links_ox) <- NULL
    rownames(links_lac) <- NULL
    
  }
  
  hist(c(hist_mat[,2]-hist_mat[,1],hist_mat[,3]-hist_mat[,2],hist_mat[,4]-hist_mat[,3]),breaks=1000,xlab="Differences between the regression coefficiences", main="Histogram")
  
  # GRNs from all data sets 
  links[,1] <- tolower(links[,1])
  links[,2] <- tolower(links[,2])
  
  neg_link <- links[which(links$weight<0),]
  pos_link <- links[which(links$weight>0),]
  
  ROC_analysis_links <- (links[,c(1,2,4)])
  ROC_analysis_links$weight <-(abs(ROC_analysis_links$weight))/
    (max(abs(ROC_analysis_links$weight)))
  
  # GRNs from cold data set
  links_cold[,1] <- tolower(links_cold[,1])
  links_cold[,2] <- tolower(links_cold[,2])
  
  neg_link_cold <- links_cold[which(links_cold$weight<0),]
  pos_link_cold <- links_cold[which(links_cold$weight>0),]
  
  ROC_analysis_cold <- (links_cold[,c(1,2,4)])
  ROC_analysis_cold$weight <-(abs(ROC_analysis_cold$weight))/
    (max(abs(ROC_analysis_cold$weight)))
  
  # GRNs from heat data set
  links_heat[,1] <- tolower(links_heat[,1])
  links_heat[,2] <- tolower(links_heat[,2])
  
  neg_link_heat <- links_heat[which(links_heat$weight<0),]
  pos_link_heat <- links_heat[which(links_heat$weight>0),]
  
  ROC_analysis_heat <- (links_heat[,c(1,2,4)])
  ROC_analysis_heat$weight <-(abs(ROC_analysis_heat$weight))/
    (max(abs(ROC_analysis_heat$weight)))
  
  # GRNs from oxidative data set
  links_ox[,1] <- tolower(links_ox[,1])
  links_ox[,2] <- tolower(links_ox[,2])
  
  neg_link_ox <- links_ox[which(links_ox$weight<0),]
  pos_link_ox <- links_ox[which(links_ox$weight>0),]
  
  ROC_analysis_ox <- (links_ox[,c(1,2,4)])
  ROC_analysis_ox$weight <-(abs(ROC_analysis_ox$weight))/
    (max(abs(ROC_analysis_ox$weight)))
  
  # GRNs from lactose data set
  links_lac[,1] <- tolower(links_lac[,1])
  links_lac[,2] <- tolower(links_lac[,2])
  
  neg_link_lac <- links_lac[which(links_lac$weight<0),]
  pos_link_lac <- links_lac[which(links_lac$weight>0),]
  
  ROC_analysis_lac <- (links_lac[,c(1,2,4)])
  ROC_analysis_lac$weight <-(abs(ROC_analysis_lac$weight))/
    (max(abs(ROC_analysis_lac$weight)))
 
}
