args <- commandArgs()
print(args)

cancer = args[6]

myfiles = list.files(paste("./",cancer,"/",sep=""))
#myfiles = myfiles[3:455]

print(myfiles)

all.cor = c()

total = 0

for (files in myfiles){
	print(paste("./",cancer,"/",files,sep=""))
	load(paste("./",cancer,"/",files,sep=""))
	total = total+length(cancer.cor.emp.sig[,1])
	cancer.cor.emp.sig = cancer.cor.emp.sig[cancer.cor.emp.sig$padj_perm<=0.1,]
	all.cor = rbind(all.cor,cancer.cor.emp.sig)

}

print("All Triplets")
print(total)

#all_NBL.cor = do.call(rbind, lapply(myfiles, function(x) load(x)))
save(all.cor,file=paste(cancer,"_lncMod_Triplets_Combined_Permuted_0420.rda",sep=""))

