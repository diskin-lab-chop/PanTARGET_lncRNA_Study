#input linreg 

#load("NBL_TFTarget_SigPairs.rda")
#tftarget$combo = paste(tftarget$target_gene,tftarget$tf_name,sep="_")

#print("okay1")

library(stringr)

args <- commandArgs()
print(args)

cancer = args[6]
num = args[7]
load(args[8])

trans = read.table("PanTarget_AllTranscript_Promoters_0420.bed")
rownames(trans) = trans$V4

#all_pairs = read.table("./PanTarget_Promoter_Motifs/PanTarget_Promoter_All_Predicted_TF_51419.bed")

load("All_Predicted_TF_112119_Updated.rda")
all_pairs$V2 = toupper(all_pairs$V2)
all_pairs$V2 = gsub(".VAR.2.","",all_pairs$V2)
all_pairs$V2 = gsub(".VAR.3.","",all_pairs$V2)

#all_pairs$combo = paste(all_pairs$V1,all_pairs$V2,sep="_") 
all_pairs$genes = trans[as.character(all_pairs$V1),5]

all_pairs$genes = as.character(all_pairs$genes)
temp = all_pairs$genes[grep("lnc-",all_pairs$genes)]
temp = gsub("lnc-","",temp)
temp = gsub('.{2}$', '', temp)
temp[str_sub(temp,-1)=="-"] = str_sub(temp[str_sub(temp,-1)=="-"],1,-2)
temp = paste("lnc",temp,sep="-")
all_pairs$genes[grep("lnc-",all_pairs$genes)] = temp

#print(head(all_pairs))

genes_anno = read.csv("/home/modia/Apexa/lncRNA_Analysis/PanTarget_StringTie_PostProcessing_Part2/Annotate_StringTie_4419/Corrected_Gencode/PanTarget_GeneAnnotation_111819_HighConf.csv",header=F)
rownames(genes_anno) = genes_anno$V1
genes_anno$V4 = as.character(genes_anno$V4)
genes_anno$V1 = as.character(genes_anno$V1)

#this is needed if comparing gene names with HUGO symbol
genes_anno[grep("NM_|NR_",genes_anno$V4),4] = genes_anno[grep("NM_|NR_",genes_anno$V4),1]


#genes_anno - only used to convert the TF names, so no need to correct to lncipedia names here 

#temporary fix for lncipedia genes in genes_anno 
lncipedia = genes_anno[grep("lnc-",genes_anno$V1),]
temp = rownames(lncipedia)
temp = gsub("lnc-","",temp)
temp = gsub('.{2}$', '', temp)
temp[str_sub(temp,-1)=="-"] = str_sub(temp[str_sub(temp,-1)=="-"],1,-2)
temp = paste("lnc",temp,sep="-")
lncipedia$V1 = temp
lncipedia$V4  = lncipedia$V1
lncipedia = lncipedia[order(lncipedia$V1),]

#start = aggregate(lncipedia$V6~lncipedia$V1,data=lncipedia,FUN=min)
#stop = aggregate(lncipedia$V7~lncipedia$V1,data=lncipedia,FUN=max)

lncipedia = lncipedia[!duplicated(lncipedia$V1),]

rownames(lncipedia) = lncipedia$V1

#lncipedia[start$`lncipedia$V1`,]$V6 = start$`lncipedia$V6`
#lncipedia[stop$`lncipedia$V1`,]$V7 = stop$`lncipedia$V7`

genes_anno = genes_anno[grep("lnc-",genes_anno$V1,invert=T),]
genes_anno = rbind(genes_anno,lncipedia)
rownames(genes_anno) = genes_anno$V1

print("okay2")

#all_pairs = all_pairs[all_pairs$V1%in%tftarget$target_gene,]
#all_pairs = all_pairs[all_pairs$V2%in%tftarget$tf_name,]
#all_pairs = all_pairs[all_pairs$combo%in%tftarget$combo,]

fpkm = data.frame(fpkm_list[as.numeric(num)])
#print(head(fpkm))

all_pairs = all_pairs[all_pairs$genes%in%rownames(fpkm),]

#print(head(all_pairs))

tf_name = genes_anno[genes_anno$V1%in%rownames(fpkm),4]

#print(head(tf_name))
print("Number of Genes")
print(length(table(as.character(tf_name))))
print(length(tf_name))

#print(head(all_pairs))

all_pairs = all_pairs[all_pairs$V2%in%tf_name,]

print("Number of Exp Genes")
print(length(table(as.character(all_pairs$genes))))
print("Number of Exp TFs")
print(length(table(as.character(all_pairs$V2))))

all_pairs$combo = paste(all_pairs$genes,all_pairs$V2,sep="_")
all_pairs= all_pairs[!duplicated(all_pairs$combo),]

print(length(all_pairs[,1]))
print("okay3")

save(all_pairs,file=paste(cancer,"TFTarget_Pairs_0420.rda",sep="_"))

