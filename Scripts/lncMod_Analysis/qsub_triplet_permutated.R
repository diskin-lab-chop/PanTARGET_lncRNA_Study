args <- commandArgs()
print(args)

cancer = args[6]
num = args[7]

print(paste("../linReg_Run_Nov2019/",cancer,"_DE_lncRNA_Expression_0420.rda",sep=""))
load(paste("../linReg_Run_Nov2019/",cancer,"_DE_lncRNA_Expression_0420.rda",sep=""))
head(rownames(DE_lncs))

for (lncs in rownames(DE_lncs)){
        temp = paste("echo 'Rscript triplet_perlncRNA_permute.R ",lncs," ",cancer," ",num,"' | qsub -cwd -l mem_free=12G -l m_mem_free=12G -l h_vmem=12G -e /dev/null -o /dev/null",sep="")
	print(temp)
	system(temp)	      
}                                                                          
