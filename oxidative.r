
oxidative <- function()
{
  # timepoints
  tps <- c(10,20,30,40,50,90)
  # number of replications
  reps <- c(1,2,3)
  # conditions
  conds <- c("c", "t")
  # making a data frame which shows the replications and conditions with the given timepoints
  exp_fac <- data.frame(dataorder=seq(1,length(reps)*length(conds)*length(tps)), conditions=rep(conds,each=length(reps)*length(tps)), 
                        replicates=rep(reps,length(tps),length(reps)),
                        tps=rep(tps,length(conds),each=length(reps)))
  
  exp_fac$dataorder=1:((length(reps)*length(conds)*length(tps)))
  # sort the exp_fac first by condition and then by replicates and tps
  exp_fac=with(exp_fac,exp_fac[order(conditions,replicates,tps),])
  
  oxidative_data <- read.table("Data/oxidativestress.txt",header=T,sep="\t",row.names=1,stringsAsFactors=F)
  control_data <- read.table("Data/control.txt",header=T,sep="\t",row.names=1,stringsAsFactors=F)
  
  frep_oxidative_data <-oxidative_data[,c((0:4*3)+1)]
  srep_oxidative_data <-oxidative_data[,c((0:4*3)+2)]
  trep_oxidative_data <-oxidative_data[,c((0:4*3)+3)]
  
  domain <- c(10,20,30,40,90)
  X <- tps
  i_frep_oxidative_data <- (cubic.spline.interpol(frep_oxidative_data,domain,X))
  fid <- which(i_frep_oxidative_data[,5]<0)
  
  i_srep_oxidative_data <- (cubic.spline.interpol(srep_oxidative_data,domain,X))
  sid <- which(i_srep_oxidative_data[,5]<0)

  i_trep_oxidative_data <- (cubic.spline.interpol(trep_oxidative_data,domain,X))
  tid <- which(i_trep_oxidative_data[,5]<0)
  
  id <- unique(union(tid,union(sid,fid)))
  
  colnames(i_frep_oxidative_data) <- colnames(i_srep_oxidative_data) <- 
    colnames(i_trep_oxidative_data) <- c("t10","t20","t30","t40","t50","t90")
  
  rownames(i_frep_oxidative_data) <- rownames(i_srep_oxidative_data) <- 
    rownames(i_trep_oxidative_data) <- rownames(oxidative_data)
  
  exp_val <- cbind(control_data, oxidative_data[,1:12],
                   t50=i_frep_oxidative_data[,5],t50=i_srep_oxidative_data[,5],t50=i_trep_oxidative_data[,5],
                   oxidative_data[,13:15])
  
  trash <- (sapply(id,function(ids){id_neg<-which(exp_val[ids,which(colnames(exp_val)%in%"t50")]<0);
                          id_pos<-which(exp_val[ids,which(colnames(exp_val)%in%"t50")]>0);
                          exp_val[ids,which(colnames(exp_val)%in%"t50")][id_neg] <<- 
                                      (exp_val[ids,which(colnames(exp_val)%in%"t50")][id_pos[1]])
                          }))
  
  colnames(exp_val) <- gsub("\\.[1-3]","", colnames(exp_val))
  
  exp_val=exp_val[,exp_fac$dataorder]
  
  colnames(exp_val) <- gsub("\\.[1-3]","", colnames(exp_val))
  
  exp_val <- log2(exp_val)
  colnames(exp_val) <- gsub("\\.[1-3]","", colnames(exp_val))
  
  # E(Yg)=X*Ag ------> we have Yg, therefor we have to create X to calculate Ag
  exp_structure=factor(colnames(exp_val))
  X=model.matrix(~0+exp_structure)
  colnames(X)=levels(exp_structure)
  # this is Ag :
  lm.fit=lmFit(exp_val,X)
  
  # we have to calculate differences between control and treatment at each timepoint j, therefore:
  mc=makeContrasts('t10-c10','t20-c20','t30-c30','t40-c40','t50-c50','t90-c90',levels=X)  
  
  c.fit=contrasts.fit(lm.fit,mc)
  
  # to calculate P.values we use moderate t-statistic with eBayes:
  eb=eBayes(c.fit)
  
  bstat <- eb$lods
  
  pr_oxidative <- (exp(bstat) / (matrix(1,nrow(bstat),ncol(bstat))+exp(bstat)))
  
  write.table(pr_oxidative,"Data/pr_oxidative.txt",row.names=T,col.names=T,sep="\t")
  
  data <- sapply(paste("t",tps,sep=""),function(name,x){apply(x[,which(colnames(x)==name)],1,median)},exp_val)
  write.table(data,"Data/oxidative.txt",row.names=T,col.names=T,sep="\t")
  
  return(pr_oxidative)
  
}