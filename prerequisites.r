
# read data
load("Data/exprData_fixed.Rdata")

# separate control data set and fix lactose dataset
control <- exprData.fixed[,(grep("control",colnames(exprData.fixed)))]
control<- control[,-(grep("_0|_6|_7",colnames(control)))]
lac <- exprData.fixed[,(grep("lac",colnames(exprData.fixed)))]
lac<- lac[,-(grep("_5",colnames(lac)))]


#read gene mapping - match between gene name and probe id
genes <- read.table("Data/gene_map.txt",header=T,sep="\t",stringsAsFactors=F)
genes$X.probe.id <- gsub("\t","",genes$X.probe.id)
genes$description <- gsub("\t","",genes$description)
rownames(genes) <- genes[,1]

#read tfs obtaine from dream5 challenge
tfs <- read.table("Data/TFS_dream5.txt",header=F,sep="\t",stringsAsFactors=F)
probable_tfs <- genes[which((genes$ORF) %in% (unique((tfs$V1)))),4]
tfs_silencing <- as.numeric(as.matrix(read.table("Data/141_TFS_silencing.txt",header=F,sep="\t",stringsAsFactors=F)))
tfs_s <- tfs[tfs_silencing,]
probable_tfs_s <- genes[which(toupper(genes$ORF) %in% (unique(toupper(tfs_s)))),4]

# the gold standard network from DREAM5
names_gs <- read.table("Data/net3_gene_ids.txt",header=F,sep="\t",stringsAsFactors=F)
gs <- read.table("Data/DREAM5_NetworkInference_GoldStandard_Network3.tsv",header=F,sep="\t",stringsAsFactors=F)
gs$V1 <- names_gs$V2[gs$V1]
gs$V2 <- names_gs$V2[gs$V2]
gs_links <- gs[which(gs$V3==1),]
g_gs <- graph.data.frame(gs_links)
length(E(g_gs))
# genes in Dream5
all_gene_D5 <- read.table(file="Data/genesInDream5.txt",header=T,sep="\t",stringsAsFactors=F)

# the experimentally verified links (network) from regulonDB
reg_db <- read.table("Data/network_tf_gene_nodup_lowercase.txt",header=T,sep="\t",stringsAsFactors=F)
reg_db_genes <- unique(tolower(c(reg_db[,1],reg_db[,2])))
genes_to_use <- intersect(intersect(tolower(names_gs$V2),tolower(reg_db_genes)),tolower(genes$ORF[(which(genes$synonym%in%rownames(control)))]))
tfs <- unique(reg_db[,1])
tfs_to_use <- genes$synonym[which(tolower(genes$ORF)%in%intersect(tolower(tfs),tolower(genes_to_use)))]
genes_to_use <- genes$synonym[which(tolower(genes$ORF)%in%tolower(genes_to_use))]

# Selecting genes
gene_InNET_D5 <- read.table(file="Data/genesParticpatingInnetwrokDream5.txt",header=F,sep="\t",stringsAsFactors=F)
gene_to_use <- genes[which(toupper(genes$ORF) %in% toupper(all_gene_D5[gene_InNET_D5[,1],"Name"])),4]

genes_to_use <- unique(union(genes_to_use,gene_to_use))
tfs_to_use <- unique(union(tfs_to_use,probable_tfs_s))

