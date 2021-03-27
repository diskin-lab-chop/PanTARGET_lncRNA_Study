library(reshape2)
library(dplyr)
library(data.table)

#input: lingreg,tftarget, and lncs - data frame of cancer gene expression and the lncRNA name 

args <- commandArgs()
print(args)

lncs = args[6]
cancer = args[7]
num = args[8]

#Expression matrix
load("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_112119/PanTarget_TotalGeneExpression_Cutoff_0420.rda")

#Significant TF-Target gene pairs
load(paste("../linReg_Run_Nov2019/",cancer,"_TFTarget_LinReg_Combined_0420.rda",sep=""))

#Expression matrix
cancer_linreg= data.frame(fpkm_list[as.numeric(num)])

tftarget = linreg_results

lncexp = sort(cancer_linreg[lncs,])
#0.25*95 (top and bottom 24 samples)

bottom_start = floor((length(cancer_linreg)*0.25))
top_end = floor(length(cancer_linreg)-(length(cancer_linreg)*0.25))
top_end = top_end+1
bottomlnc = colnames(lncexp[c(1:bottom_start)])
toplnc = colnames(lncexp[c(top_end:length(cancer_linreg))])


#label each gene as TF or lnc or PCG (remove malat1 from table)

#cancer.dat = cancer_linreg[,colnames(cancer_linreg)%in%c(toplnc,bottomlnc)]
cancer.dat = cancer_linreg[,colnames(cancer_linreg)%in%c(toplnc)]
cancer.dat$Name = rownames(cancer.dat)

#randomize the target gene samples to break the association of y to TF (x variable)

#removes half to reduce the calculations 
#temp = names(table(tftarget$tfs))
#tftarget_1 = tftarget[tftarget$tfs%in%temp[1:550],]

#only keep the TFs that dont target the lncRNA in question 
tflncs = tftarget[tftarget$target_gene==lncs,]
tftarget_1 = tftarget[!(tftarget$tfs%in%tflncs$tfs),]
tftarget_1 = tftarget_1[,c(1,2)]
#cancer.dat = melt(cancer.dat,id.vars = c("Name"))

#determine how correlated the tftargets are in the bottom, in the top, and the difference
colnames(tftarget_1)[1] = "Name"
cancer.TF = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.TF = melt(cancer.TF)
cancer.TF$combo = paste(cancer.TF$Name,cancer.TF$target_gene,cancer.TF$variable,sep="_")
cancer.TF = cancer.TF[order(cancer.TF$combo),]
cancer.TF$combo = NULL

colnames(tftarget_1)[1] = "tfs"
colnames(tftarget_1)[2] = "Name"
cancer.target = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.target = melt(cancer.target)
cancer.target$combo = paste(cancer.target$tfs,cancer.target$Name,cancer.target$variable,sep="_")
cancer.target = cancer.target[order(cancer.target$combo),]
cancer.target$combo = NULL


cancer.join = cbind(cancer.TF,cancer.target$value)
colnames(cancer.join)[4] = "TF.Exp"
colnames(cancer.join)[5] = "Target.Exp"
cancer.join = cancer.join %>%as.data.table()


cancer.cor.top = cancer.join[,.(cor=cor.test(TF.Exp,Target.Exp,method="spearman")$estimate),by=c("Name","target_gene")]
cancer.pval.top = cancer.join[,.(pval=cor.test(TF.Exp,Target.Exp,method="spearman")$p.value),by=c("Name","target_gene")]
cancer.cor.top = merge(cancer.cor.top,cancer.pval.top,by=c("Name","target_gene"))

warnings()
print(warnings())
rm(cancer.TF,cancer.target,cancer.join)

#-------------------calculate bottom correlation ---------------------------

cancer.dat = cancer_linreg[,colnames(cancer_linreg)%in%c(bottomlnc)]

print("Top and Bottom Length")
print(length(cancer.dat))

cancer.dat$Name = rownames(cancer.dat)

#removes half to reduce the calculations 
#temp = names(table(tftarget$tfs))
#tftarget_1 = tftarget[tftarget$tfs%in%temp[1:550],]

#only keep the TFs that dont target the lncRNA in question 
tflncs = tftarget[tftarget$target_gene==lncs,]
tftarget_1 = tftarget[!(tftarget$tfs%in%tflncs$tfs),]
tftarget_1 = tftarget_1[,c(1,2)]

#cancer.dat = melt(cancer.dat,id.vars = c("Name"))

#determine how correlated the tftargets are in the bottom, in the top, and the difference
colnames(tftarget_1)[1] = "Name"
cancer.TF = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.TF = melt(cancer.TF)
cancer.TF$combo = paste(cancer.TF$Name,cancer.TF$target_gene,cancer.TF$variable,sep="_")
cancer.TF = cancer.TF[order(cancer.TF$combo),]
cancer.TF$combo = NULL

colnames(tftarget_1)[1] = "tfs"
colnames(tftarget_1)[2] = "Name"
cancer.target = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.target = melt(cancer.target)
cancer.target$combo = paste(cancer.target$tfs,cancer.target$Name,cancer.target$variable,sep="_")
cancer.target = cancer.target[order(cancer.target$combo),]
cancer.target$combo = NULL

cancer.join = cbind(cancer.TF,cancer.target$value)
colnames(cancer.join)[4] = "TF.Exp"
colnames(cancer.join)[5] = "Target.Exp"
cancer.join = cancer.join %>%as.data.table()

cancer.cor.bottom = cancer.join[,.(cor=cor.test(TF.Exp,Target.Exp,method="spearman")$estimate),by=c("Name","target_gene")]
cancer.pval.bottom = cancer.join[,.(pval=cor.test(TF.Exp,Target.Exp,method="spearman")$p.value),by=c("Name","target_gene")]
cancer.cor.bottom = merge(cancer.cor.bottom,cancer.pval.bottom,by=c("Name","target_gene"))

cancer.cor = cancer.cor.top %>% select(TF="Name",TargetGene= "target_gene",Top_cor = "cor",Top_pval="pval") %>% left_join(select(cancer.cor.bottom,TF="Name",TargetGene= "target_gene",Bottom_cor="cor",Bottom_pval="pval"),by=c("TF","TargetGene"))
  
cancer.cor = cancer.cor[cancer.cor$TF!=cancer.cor$TargetGene,]
cancer.cor$Delta = abs(cancer.cor$Top_cor-cancer.cor$Bottom_cor)
print(length(cancer.cor[,1]))
print(head(cancer.cor$Top_pval))
print(length(cancer.cor$Top_pval))
cancer.cor$Top_padj = p.adjust(cancer.cor$Top_pval,method="BH")
cancer.cor$Bottom_padj = p.adjust(cancer.cor$Bottom_pval,method="BH")
cancer.cor$lncRNA = lncs

cancer.cor.emp = cancer.cor
rm(cancer.cor)
cancer.cor.emp.sig = cancer.cor.emp[cancer.cor.emp$Delta>=0.45,]

#make sure either the top or bottom correlation is significant 
cancer.cor.emp.sig = cancer.cor.emp.sig[abs(cancer.cor.emp.sig$Top_cor)>0.4 | abs(cancer.cor.emp.sig$Bottom_cor)>0.4,]

cancer.cor.emp.sig$fisher_top = 0.5*log((1+cancer.cor.emp.sig$Top_cor)/(1-cancer.cor.emp.sig$Top_cor))
cancer.cor.emp.sig$fisher_bottom = 0.5*log((1+cancer.cor.emp.sig$Bottom_cor)/(1-cancer.cor.emp.sig$Bottom_cor))

total_len = (length(cancer.dat[,grep("TARGET",colnames(cancer.dat))]))-3
print(total_len)
#total_len = 38

cancer.cor.emp.sig$rewire = abs((cancer.cor.emp.sig$fisher_top-cancer.cor.emp.sig$fisher_bottom)/(sqrt((1*2)/total_len)))

cancer.cor.emp.sig$rewire_pvalue = (2*(1-pnorm(abs(cancer.cor.emp.sig$rewire))))

cancer.cor.emp.sig$combo = paste(cancer.cor.emp.sig$TF,cancer.cor.emp.sig$TargetGene,sep="_")

rm(cancer.cor.bottom,cancer.cor.top,cancer.dat,cancer.join,cancer.pval.bottom,cancer.pval.top,cancer.target,cancer.TF)

print(head(cancer.cor.emp.sig))
print("Triplets Pass Inital Threshold")
print(length(cancer.cor.emp.sig$combo))
#save(cancer.cor.emp.sig,file=paste("./",cancer,"/",lncs,"_",cancer,"_Triplet_Permutated_NumFP_Spear.rda",sep=""))

#--------------------Permutation---------------------------

cancer.dat = cancer_linreg[,colnames(cancer_linreg)%in%c(bottomlnc,toplnc)]
#order them based on the low to high expression
cancer.dat = cancer.dat[,c(bottomlnc,toplnc)]
cancer.dat$Name = rownames(cancer.dat)

#only keep the TFs that dont target the lncRNA in question 
tflncs = tftarget[tftarget$target_gene==lncs,]
tftarget_1 = tftarget[!(tftarget$tfs%in%tflncs$tfs),]
tftarget_1 = tftarget_1[,c(1,2)]

#determine how correlated the tftargets are in the bottom, in the top, and the difference
colnames(tftarget_1)[1] = "Name"
cancer.TF = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.TF = melt(cancer.TF)
cancer.TF$combo = paste(cancer.TF$Name,cancer.TF$target_gene,cancer.TF$variable,sep="_")
cancer.TF$combo2 = paste(cancer.TF$Name,cancer.TF$target_gene,sep="_")
cancer.TF = cancer.TF[cancer.TF$combo2%in%cancer.cor.emp.sig$combo,]
cancer.TF = cancer.TF[order(cancer.TF$combo),]
cancer.TF$combo = NULL
cancer.TF$combo2 = NULL

colnames(tftarget_1)[1] = "tfs"
colnames(tftarget_1)[2] = "Name"
cancer.target = merge(cancer.dat,tftarget_1,by=c("Name"))
cancer.target = melt(cancer.target)
cancer.target$combo = paste(cancer.target$tfs,cancer.target$Name,cancer.target$variable,sep="_")
cancer.target$combo2 = paste(cancer.target$tfs,cancer.target$Name,sep="_")
cancer.target = cancer.target[cancer.target$combo2%in%cancer.cor.emp.sig$combo,]
cancer.target = cancer.target[order(cancer.target$combo),]
cancer.target$combo = NULL
cancer.target$combo2 = NULL

cancer.join = cbind(cancer.TF,cancer.target$value)
colnames(cancer.join)[4] = "TF.Exp"
colnames(cancer.join)[5] = "Target.Exp"
cancer.join = cancer.join %>%as.data.table()
cancer.join$Label = "None"
cancer.join[cancer.join$variable%in%bottomlnc,]$Label = "Bottom"
cancer.join[cancer.join$variable%in%toplnc,]$Label = "Top"
#cancer.join$combo = paste(cancer.join$Name,cancer.join$target_gene,sep="_")
#cancer.join = cancer.join[cancer.join$combo%in%cancer.cor.emp.sig$combo,]


print(head(cancer.join))
print(length(cancer.join$Label))
#perform the shuffling then split to perform the correlation 

all_pvalues = c()
all_pvalues2 = c()
print(Sys.time())

for (ii in 1:100){
	#print("Step0")
	#print(Sys.time())
	cancer.perm <- cancer.join %>% group_by(Name,target_gene) %>% mutate(Target.Exp=sample(Target.Exp)) %>% as.data.table()
        cancer.join.top = cancer.perm[cancer.perm$Label=="Top",]
        cancer.join.bottom = cancer.perm[cancer.perm$Label=="Bottom",]
	cancer.cor.top = cancer.join.top[,.(cor=cor.test(TF.Exp,Target.Exp,method="spearman")$estimate),by=c("Name","target_gene")]
	cancer.cor.bottom = cancer.join.bottom[,.(cor=cor.test(TF.Exp,Target.Exp,method="spearman")$estimate),by=c("Name","target_gene")]
	#print("Step1")
	#print(Sys.time())
	cancer.cor = cancer.cor.top %>% select(TF="Name",TargetGene= "target_gene",Top_cor = "cor") %>% left_join(select(cancer.cor.bottom,TF="Name",TargetGene= "target_gene",Bottom_cor="cor"),by=c("TF","TargetGene"))
	#print("Step2")
	#print(Sys.time())
	cancer.cor$Delta = abs(cancer.cor$Top_cor-cancer.cor$Bottom_cor)
	#filter out triplets that dont pass the threshold
	#cancer.cor.sig = cancer.cor[cancer.cor$Delta>=0.45,]
	#make sure either the top or bottom correlation is significant 
	#cancer.cor.sig = cancer.cor.sig[abs(cancer.cor.sig$Top_cor)>0.4 | abs(cancer.cor.sig$Bottom_cor)>0.4,]	
	cancer.cor.sig =cancer.cor
	cancer.cor.sig$fisher_top = 0.5*log((1+cancer.cor.sig$Top_cor)/(1-cancer.cor.sig$Top_cor))
	cancer.cor.sig$fisher_bottom = 0.5*log((1+cancer.cor.sig$Bottom_cor)/(1-cancer.cor.sig$Bottom_cor))
	total_len = (length(cancer.dat[,grep("TARGET",colnames(cancer.dat))])/2)-3
	#print(total_len)
	#total_len = 38
	print(ii)
	#cancer.cor.sig$combo = paste(cancer.cor.sig$TF,cancer.cor.sig$TargetGene,sep="_")
	#cancer.cor.sig = cancer.cor.sig[cancer.cor.sig$combo%in%cancer.cor.emp.sig$combo,]
	#print(length(cancer.cor.sig$combo))
	cancer.cor.sig$rewire = abs((cancer.cor.sig$fisher_top-cancer.cor.sig$fisher_bottom)/(sqrt((1*2)/total_len)))
	cancer.cor.sig = cancer.cor.sig[order(cancer.cor.sig$rewire,decreasing=T),]
	cancer.cor.sig$rank = 1:length(cancer.cor.sig$rewire)
	cancer.cor.sig$rewire_pvalue = (2*(1-pnorm(abs(cancer.cor.sig$rewire))))
	#all_pvalues = c(all_pvalues,min(cancer.cor.sig$rewire_pvalue))
	#n = 11
	all_pvalues = c(all_pvalues,min(cancer.cor.sig$rewire_pvalue)) #group this by target because the permuation should be lncRNA to target gene  

	#all_pvalues = c(all_pvalues,cancer.cor.sig[cancer.cor.sig$rank==n,]$rewire_pvalue)
}

print(total_len)

#all_pvalues = data.frame(all_pvalues)
#colnames(all_pvalues) = paste("V",c(1:100),sep="")
#colnames(all_pvalues) = paste("V",c(1:2),sep="")
#rownames(all_pvalues) = cancer.cor.sig$combo

#all_pvalues = all_pvalues[rownames(all_pvalues)%in%cancer.cor.emp.sig$combo,]
#all_pvalues = all_pvalues[as.character(cancer.cor.emp.sig$combo),]

print(all_pvalues)

#print(all_pvalues2)

#all_pvalues =all_pvalues2

cancer.cor.emp.sig$padj_perm = unlist(lapply(cancer.cor.emp.sig$rewire_pvalue,function(x) (length(all_pvalues[all_pvalues<=as.numeric(x)])+1)/(length(all_pvalues)+1)))


#cancer.cor.emp.sig$padj_perm = c()

#for (ii in 100){
#	adj_pvalue = length(all_pvalues[ii,][all_pvalues[ii,]>cancer.cor.emp.sig$rewire_pvalue[ii]])/(length(all_pvalues[ii]))
#	cancer.cor.emp.sig$padj_perm = c(cancer.cor.emp.sig$padj_perm,adj_pvalue)
#}

save(cancer.cor.emp.sig,file=paste("./",cancer,"/",lncs,"_",cancer,"_Triplet_Permutated_Spear_0420.rda",sep=""))


