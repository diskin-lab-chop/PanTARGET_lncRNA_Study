args <- commandArgs()
print(args)

cancer = args[6]
num = args[7]

#Expression matrix
load("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/PanTarget_TotalGeneExpression_Cutoff.rda")


#TF Target Pairs
print(paste("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/",cancer,"_TFTarget_Pairs.rda",sep=""))
load(paste("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/",cancer,"_TFTarget_Pairs.rda",sep=""))

genes_anno = read.csv("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/PanTarget_GeneAnnotation_4419_HighConf_V2.csv",header=F)
rownames(genes_anno) = genes_anno$V1
noncoding = genes_anno[grep("lncRNA",genes_anno$V3),]
noncoding = noncoding[!(noncoding$V4%in%c("RMRP","RPPH1","RN7SL832P","RNU11","RNU12","RNU6ATAC35P","NR_003003","SCARNA15","SCARNA17","SCARNA2","SCARNA9")),]


cancer_tf = genes_anno[genes_anno$V4%in%all_pairs$V2,1]
cancer_tf = cancer_tf[!duplicated(cancer_tf)]

linreg = data.frame(fpkm_list[as.numeric(num)])

#Keep only the genes involved in the TFTarget Pairs for this step
linreg = linreg[rownames(linreg)%in%all_pairs$genes,]
length(linreg[,1])

#Further only keep genes that have differential expression (have highly variable expression in the data set).

DE_genes = apply(linreg,1,function(x) IQR(x))
DE_genes = log2(DE_genes)
DE_genes = DE_genes[DE_genes>1]
temp1 = linreg[rownames(linreg)%in%names(DE_genes),]
#remove any of the DE TF's - all TF's will be added next
temp1 = temp1[!(rownames(temp1)%in%cancer_tf),]
#get the matrix for all TF's - want to include even TF's that are not DE
temp2 = linreg[rownames(linreg)%in%cancer_tf,]
linreg = rbind(temp1,temp2)

length(linreg[,1])

DE_lncs = linreg[rownames(linreg)%in%rownames(noncoding),]
save(DE_lncs,file=paste(cancer,"_DE_lncRNA_Expression.rda",sep=""))

#linreg, all_pair

count = 0
for (genes in rownames(linreg)){
	break
	temp = paste("echo \'Rscript linreg_for_tftarget.R",genes,cancer,num,"\' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e /dev/null -o /dev/null",sep=" ")
	print(temp)
	#system(temp)
	count = count+1
	print(count)
	#break
}
