
# R version 3.1.1 (2014-07-10) -- "Sock it to Me"
# Copyright (C) 2014 The R Foundation for Statistical Computing
# Platform: x86_64-redhat-linux-gnu (64-bit)

# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.

# Natural language support but running in an English locale

# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.

# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.

main <- function()
{  
  # load all necessary libraries, functions and prerequisite data from DREAM5 and regulonDB
  #source("sources.r")

  # find probability of being differentially expressed at each specific time point
  Pr_h <- heat()
  Pr_c <- cold()
  Pr_o <- oxidative()
  Pr_l <- lactose(lac,control)
  
  rownames(Pr_h) <- genes[rownames(Pr_h),4]
  rownames(Pr_c) <- genes[rownames(Pr_c),4]
  rownames(Pr_o) <- genes[rownames(Pr_o),4]
  
  
  # timepoints
  tps <- c(10,20,30,40,50,90)
  
  # Find the similarity matrix between genes with respect to the probability of being differentially expressed (calculated separately for each data set)
  pr_heat <- apply(Pr_h,1,function(r){res <- c(); for (j in 1:nrow(Pr_h))
                                            res<-c(res,1/length(tps)* sum(abs(r-Pr_h[j,])))
                           return(res)})
  
  rownames(pr_heat) <- colnames(pr_heat) <- tolower(colnames(pr_heat))
  
  pr_cold <- apply(Pr_c,1,function(r){res <- c(); for (j in 1:nrow(Pr_c))
    res<-c(res,1/length(tps)* sum(abs(r-Pr_c[j,])))
    return(res)})
  
  rownames(pr_cold) <- colnames(pr_cold)
  
  pr_oxidative <- apply(Pr_o,1,function(r){res <- c(); for (j in 1:nrow(Pr_o))
    res<-c(res,1/length(tps)* sum(abs(r-Pr_o[j,])))
    return(res)})
  
  rownames(pr_oxidative) <- colnames(pr_oxidative)
  
  # timepoints
  tps <- c(0,10,20,30,40)
  pr_lactose <- apply(Pr_l,1,function(r){res <- c(); for (j in 1:nrow(Pr_l))
    res<-c(res,1/length(tps)* sum(abs(r-Pr_l[j,])))
                                           return(res)})
  
  rownames(pr_lactose) <- colnames(pr_lactose)

  # saving the weight matrices
  write.table(pr_heat,"Data/weight_heat.txt",row.names=T,col.names=T,sep="\t")
  write.table(pr_cold,"Data/weight_cold.txt",row.names=T,col.names=T,sep="\t")
  write.table(pr_oxidative,"Data/weight_oxidative.txt",row.names=T,col.names=T,sep="\t")
  write.table(pr_lactose,"Data/weight_lactose.txt",row.names=T,col.names=T,sep="\t")

  # saving the probability matrices
  save(pr_heat,pr_cold,pr_oxidative,pr_lactose,file="Data/probability_mats.obj")    
  
  ################################
  
  # Probablities to calc W times X
  load("Data/probability_mats.obj")
  
  # K datasets which is here 4
  load("Data/datasets.obj")
  
  ################ preparing data in the format of inputs to the penalized function
  
  ##### seclect TFs as regressors
  regressors <- as.matrix(bdiag(t(heat[tfs_to_use,]),t(cold[tfs_to_use,]),t(oxidative[tfs_to_use,]),t(lactose[tfs_to_use,])))
  rownames(regressors) <- c(colnames(heat),colnames(cold),colnames(oxidative),colnames(lactose))
  colnames(regressors) <- rep(tfs_to_use,4)
  
  ##### select genes + TFs for responses
  responses <- as.matrix(cbind(heat[genes_to_use,],cold[genes_to_use,],oxidative[genes_to_use,],lactose[genes_to_use,]))
  
  ##### Select pr mats based on the gene used in regressors and responses
  ph<-pr_heat[genes_to_use,tfs_to_use]
  pc<-pr_cold[genes_to_use,tfs_to_use]
  po<-pr_oxidative[genes_to_use,tfs_to_use]
  pl<-pr_lactose[genes_to_use,tfs_to_use]
  
  ##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Make sure that you substituted the the function in fused.lasso.modified.r with the original one through following steps:
  # 1: trace(fused.lasso,edit=T)
  # 2: substitute the function in fused.lasso.modified.r with the existing  one
  
  ##### performing the regression model on data
  fits_penalized <- (mclapply(genes_to_use,
                              function(i,res,reg){
                                print(i)
                                y<-as.numeric(res[i,]);
                                if (i %in% colnames(reg))
                                {
                                  ib<-diag(c(ph[i,-which(colnames(ph)==i)],pc[i,-which(colnames(pc)==i)],po[i,-which(colnames(po)==i)],pl[i,-which(colnames(pl)==i)]))
                                  x <-reg[,-which(colnames(reg)==i)];
                                  colnames(x) <-colnames(reg)[-which(colnames(reg)==i)]
                                }
                                else
                                {
                                  ib<-(diag(c(ph[i,],pc[i,],po[i,],pl[i,])))
                                  x <- reg;
                                  colnames(x) <-colnames(reg);
                                }
                                fit <- cv.lqa(y.train=y,x.train=x,penalty.family=fused.lasso,intercept=F,...=diag(ib),standardize=T,
                                              lambda.candidates=list(c(0.05,0.1,0.5,1,1.5),c(0.1,0.5,1,1.5,2)),n.fold=10,loss.func="dev.loss",family=gaussian())
                                mat <- matrix(fit$beta.opt,ncol=4,nrow=ncol(x)/4,byrow=F);
                                rownames(mat) <- colnames(x)[1:(ncol(x)/4)];
                                print(i)
                                write.table(mat,file=paste("Results/mat/",i,"_lqa.xls",sep=""),row.names=T,col.names=T,sep="\t");
                                save(fit,file=paste("Results/obj/",i,"_lqa_at_each_tp.obj",sep=""));
                                return(fit);
                              },responses,regressors,mc.preschedule=T,mc.set.seed=T,mc.cores=8,mc.cleanup=T))
  
  # reading the results
  read_results()
}

