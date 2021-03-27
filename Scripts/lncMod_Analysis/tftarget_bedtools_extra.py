import os

for files in os.listdir("/home/modia/Apexa/TFTarget_JASPAR/TF_Bed"):
	print("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_112119_extra_novel_promoter_0420.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > ./PanTarget_Promoter_Motifs_Extra/PanTarget_Promoter_Extra"+files+"' | qsub -cwd -l mem_free=30G -l m_mem_free=30G -l h_vmem=30G")	
	os.system("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_112119_extra_novel_promoter_0420.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > ./PanTarget_Promoter_Motifs_Extra/PanTarget_Promoter_Extra"+files+"' | qsub -cwd -l mem_free=30G -l m_mem_free=30G -l h_vmem=30G")
