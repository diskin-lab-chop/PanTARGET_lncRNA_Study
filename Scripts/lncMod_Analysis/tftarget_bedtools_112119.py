import os

for files in os.listdir("/home/modia/Apexa/TFTarget_JASPAR/TF_Bed"):
	print("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_112119_PCG_lnc_promoter.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > PanTarget_Promoter_"+files+"' | qsub -cwd -l mem_free=50G -l m_mem_free=50G -l h_vmem=50G")	
	os.system("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_112119_PCG_lnc_promoter.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > PanTarget_Promoter_"+files+"' | qsub -cwd -l mem_free=50G -l m_mem_free=50G -l h_vmem=50G")
