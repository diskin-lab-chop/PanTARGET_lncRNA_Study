library(broom)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(reshape2)

args <- commandArgs()
#print(args)
genes = args[6]
cancer = args[7]
num = args[8]

#Expression matrix
load("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/PanTarget_TotalGeneExpression_Cutoff.rda")

#TF Target Pairs - Get the pair specific to that gene
#print(paste("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/",cancer,"/",genes,"_",cancer,"_TFTarget_Pairs.rda",sep=""))
load(paste("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/",cancer,"_TFTarget_Pairs.rda",sep=""))

genes_anno = read.csv("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419/PanTarget_GeneAnnotation_4419_HighConf_V2.csv",header=F)
rownames(genes_anno) = genes_anno$V1

#Expression matrix
cancer_linreg= data.frame(fpkm_list[as.numeric(num)])
cancer_linreg$Name = genes_anno[as.character(rownames(cancer_linreg)),4]

#Get the TFs expressed in that cancer
cancer_tf = names(table(as.character(all_pairs[all_pairs$genes==genes,]$V2)))
#print(paste(genes,length(cancer_tf)))
#Get the expression matrix for the TF's
cancer_tf = rownames(cancer_linreg[cancer_linreg$Name%in%cancer_tf,])
#print(paste(genes,length(cancer_tf)))
cancer_linreg$Name = NULL


#Here we select one gene and perform linear regression for combinations of TFs
#By using tidy we are able to perform the multiple analyses all at once after melting the data frame

geneexp = data.frame(t(cancer_linreg[genes,]))
geneexp$Var1 = rownames(geneexp)
datasubset = melt(t(cancer_linreg[cancer_tf,]),id.vars = colnames(cancer_linreg))
datasubset = merge(datasubset,geneexp,by=c("Var1"))
rm(geneexp)
colnames(datasubset)[4] = "expression"
linear_models <- datasubset %>% group_by(Var2) %>% do(tidy(lm(expression ~ value, .)))
linear_models$Target_Gene = genes
linear_models$padj = p.adjust(linear_models$p.value,method="hochberg")
linear_models$Var2 = as.character(linear_models$Var2)
linear_models = linear_models[linear_models$term=="value",]
#if(empty(linear_models[linear_models$term=="value" & linear_models$padj<0.01,])!=TRUE){

if(empty(linear_models)!=TRUE){
	tfs = as.character(linear_models$Var2)
	target_gene = as.character(linear_models$Target_Gene)
	pvalue = linear_models$p.value
	padj = linear_models$padj
}

temp = data.frame(tfs,target_gene,pvalue,padj)
save(temp,file=paste("./",cancer,"/",cancer,"_TFTarget_",genes,"_LinReg.rda",sep=""))


