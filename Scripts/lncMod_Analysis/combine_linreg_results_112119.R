args <- commandArgs()
cancer = args[6]

myfiles = list.files(paste("./",cancer,"_Permute",sep=""))

#myfiles= myfiles[2:4]
#print(myfiles)
#linreg_results = do.call(rbind, lapply(myfiles, function(x) load(x)))
load(paste("/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_112119/",cancer,"_TFTarget_Pairs_0420.rda",sep=""))

tfgene = table(as.character(all_pairs$genes))

linreg_results = c()
total_length = 0
num_files = 0
total_length_match = 0
for (files in myfiles){
	num_files = num_files+1
	load(paste("./",cancer,"_Permute/",files,sep=""))
	#print(temp)
	total_length = total_length+length(temp[,1])
	if(length(tfgene[temp[,2]]) == length(temp[,1])){
		
		total_length_match = total_length_match+1

	}
	temp$BH = p.adjust(temp$pvalue,method="BH")
	temp = temp[as.numeric(temp$BH)<0.00001,]
	linreg_results = rbind(linreg_results,temp)
}
print("Files")
print(num_files)

print("Counted")
print(total_length_match)

print("All Possible Interactions")
print(total_length)

head(linreg_results)
print(dim(linreg_results))
#save(linreg_results,file=paste(cancer,"_TFTarget_LinReg_Combined_0420.rda",sep=""))
