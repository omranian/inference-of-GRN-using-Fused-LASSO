
lactose <- function(lac,control)
{
  # timepoints
  tps <- c(0,10,20,30,40)
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
  
  colnames(lac) <- paste("t",rep(tps,length(reps)),sep="")
  colnames(control) <- paste("c",rep(tps,length(reps)),sep="")
    
  lac_data <- lac
  control_data <- control
  
  exp_val <- cbind(control_data, lac_data)
  
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
  mc=makeContrasts('t10-c10','t20-c20','t30-c30','t40-c40',levels=X)  
  
  c.fit=contrasts.fit(lm.fit,mc)
  
  # to calculate P.values we use moderate t-statistic with eBayes:
  eb=eBayes(c.fit)
  
  bstat <- eb$lods
  
  pr_lactose <- (exp(bstat) / (matrix(1,nrow(bstat),ncol(bstat))+exp(bstat)))
  
  write.table(pr_lactose,"Data/pr_lac.txt",row.names=T,col.names=T,sep="\t")
  
  data <- sapply(paste("t",tps,sep=""),function(name,x){apply(x[,which(colnames(x)==name)],1,median)},exp_val)
  write.table(data,"Data/lac.txt",row.names=T,col.names=T,sep="\t")
  
  return(pr_lactose)
  
}