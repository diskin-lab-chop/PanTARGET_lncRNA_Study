README for PanTarget lncMod analysis

-------------Cluster analysis/scripts--------------

Directory: 1) /home/modia/Apexa/TFTarget_JASPAR 2) /home/modia/Apexa/lncRNA_Analysis/PanTarget_StringTie_PostProcessing_Part2/lncMod_PanTarget


Things to consider: 

Determine all TF binding in a particular cancer
	Options: 
	A)
		1. ChIP-seq for all available TF's in a representative cell line or patient samples
		2. Motif analysis + Open chromatin + TF expression - restricted to promoter regions of genes - This is described in lab notebook on 3/4/19 and will be detailed with any modifications below.
		3. Motif analysis + TF expression 
		4. Use expression and the ARACNe-AP method to infer regulator - target gene relationships (goes beyond just binding and motif and does expression inference - method is black boxish)
	B)
		After associating pair, you need to perform linear regression to determine if the target gene expression is regulated/correlates with expression of the TF

The option that was chosen: A2 - subbing open chromatin with promoter region to enable pantarget to be considered instead of just the cancers that have available DNAse/ATAC-seq data. 

JASPAR data downloaded from: http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/ 
These are all human. 

How the motifs are scored and determined to be present: http://genome.ucsc.edu/cgi-bin/hgc?hgsid=755963619_AW58u9aEPu0jFlJAXTbHVSazpAJv&c=chr2&l=16078590&r=16085037&o=16081692&t=16081705&g=hub_186875_JasparTFBS&i=ASCL1 - According to this - only motifs with FIMO p-value less than 0.05 were kept. 


Steps in lncMod as performed on Respublica cluster and modified from initial run as described in lab notebook 3/4/19:

Step 1:

Overlap TF motif to target gene promoter region. 

This involves motifs for TFs globally across the genome: /home/modia/Apexa/TFTarget_JASPAR/TF_Bed

import os

for files in os.listdir("/home/modia/Apexa/TFTarget_JASPAR/TF_Bed"):
        print("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_5119_PCG_lnc_promoter.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > PanTarget_Promoter_"+files+"' | qsub -cwd -l mem_free=50G -l m_mem_free=50G -l h_vmem=50G")
        os.system("echo 'bedtools intersect -a PanTarget_TranscriptAnnotation_5119_PCG_lnc_promoter.bed -b ../TF_Bed/"+files+" -wo | cut -f 4,10 | uniq > PanTarget_Promoter_"+files+"' | qsub -cwd -l mem_free=50G -l m_mem_free=50G -l h_vmem=50G")
   
This overlap found here:
/home/modia/Apexa/TFTarget_JASPAR/PanTarget_Promoter_Motifs_51419

python tftarget_bedtools.py

Note: Strand not considered for overlap, given that a TF regulating a gene can regulate from either strand, not just strand of the gene in question. 

Combine promoter motif overlap results: 
cat PanTarget_Promoter_*tsv.bed > PanTarget_Promoter_All_Predicted_TF_51419.bed

Get the pairs all into one file and subset based on the genes that are actually expressed in that cancer subtype. 

Merge all the results in the promoter_overlap_files folder:

	echo 'Rscript reduce_pairs_pantarget_112119.R AML 1 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o AML_Reduce.log -e AML_Reduce.error
	echo 'Rscript reduce_pairs_pantarget_112119.R TALL 2 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o TALL_Reduce.log -e TALL_Reduce.error
	echo 'Rscript reduce_pairs_pantarget_112119.R BALL 3 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o BALL_Reduce.log -e BALL_Reduce.error
	echo 'Rscript reduce_pairs_pantarget_112119.R RT 4 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o RT_Reduce.log -e RT_Reduce.error
	echo 'Rscript reduce_pairs_pantarget_112119.R WT 5 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o WT_Reduce.log -e WT_Reduce.error
	echo 'Rscript reduce_pairs_pantarget_112119.R NBL 6 PanTarget_TotalGeneExpression_Cutoff_112119.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o NBL_Reduce.log -e NBL_Reduce.error

	Updated: echo 'Rscript reduce_pairs_pantarget_112119.R NBL 6 PanTarget_TotalGeneExpression_Cutoff_030520.rda' |  qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -o NBL_Reduce.log -e NBL_Reduce.error

Step 2:

/home/modia/Apexa/lncRNA_Analysis/PanTarget_StringTie_PostProcessing_Part2/linReg_Run_Nov2019 - Run linear regression for all genes against all relevant TFs: 

echo 'Rscript PanTarget_linreg_qsub2.R AML 1' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G
echo 'Rscript PanTarget_linreg_qsub2.R TALL 2' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G
echo 'Rscript PanTarget_linreg_qsub2.R BALL 3' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G
echo 'Rscript PanTarget_linreg_qsub2.R RT 4' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G
echo 'Rscript PanTarget_linreg_qsub2.R WT 5' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G
echo 'Rscript PanTarget_linreg_qsub2.R NBL 6' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G

Requires the following script: linreg_for_tftarget.R and the expression and tf target pair data: Expression data frame per cancer, genes anno, and TF Target Pairs matrix.

All of these qsub this per gene: 
Example:
echo 'Rscript linreg_for_tftarget.R ENSG00000000003.10 NBL 6 ' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G

Then run: combine_linreg_results_112119.R - to combine all the results, here we removed an padj < 0.00001 (idea being that if each Target gene has 200 TFs to consider and 7000-10000 target genes are considered then this gives us about 2 false positives) as an initial filter and to reduce the size of the combined output file: All_NBL_LinReg_Results.rda. Additionally, we notice that some target genes can have up to 70-100 TFs that appear to be highly correlated in expression. (This be because of indirect regulation or pathway similarity as opposed to direct binding, additionally certain families of transcription factors tend to have similar motifs.)

echo 'Rscript combine_linreg_results_112119.R AML' | qsub -cwd - only takes 1-5 minutes to run - run per cancer.


Step 3: 

Run the lncMod aspect - sorting the patients based on lncRNA expression and then perform correlation in the top/bottom 25% - cutoffs will be determined 

Folder: /home/modia/Apexa/lncRNA_Analysis/PanTarget_StringTie_PostProcessing_Part2/lncMod_PanTarget/lncRNA_Triplets_Run_Nov2019

echo 'Rscript qsub_triplet.R RT 4' | qsub -cwd -e RT_Qsub_Triplet.error -o RT_Qsub_Triplet.log
echo 'Rscript qsub_triplet.R AML 1' | qsub -cwd -e AML_Qsub_Triplet.error -o AML_Qsub_Triplet.log
echo 'Rscript qsub_triplet.R TALL 2' | qsub -cwd -e TALL_Qsub_Triplet.error -o TALL_Qsub_Triplet.log
echo 'Rscript qsub_triplet.R BALL 3' | qsub -cwd -e BALL_Qsub_Triplet.error -o BALL_Qsub_Triplet.log
echo 'Rscript qsub_triplet.R WT 5' | qsub -cwd -e WT_Qsub_Triplet.error -o WT_Qsub_Triplet.log
echo 'Rscript qsub_triplet.R NBL 6' | qsub -cwd -e NBL_Qsub_Triplet.error -o NBL_Qsub_Triplet.log

Requires: triplet_perlncRNA.R -which requires the output files from Step2: All_NBL_LinReg_Results.rda - the significant value- FDR threshold set at 1e-5 for 1million interactions this would then be like 10 false positives; 2nd file is the sample expression: NBL_lncMod_Expression.rda

Example of what gets run from the qsub_triplet script: echo 'Rscript triplet_perlncRNA_permute.R ENSG00000130600.11 NBL 6' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G - This is per lncRNA examined. 

This article explains how to do permutation based p-values while controlling for a certain level of FDR or number of false positives.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2592303/

echo 'Rscript qsub_triplet_permutated.R AML 1' | qsub -cwd -e AML_Qsub_Triplet_Permutated_91019.error -o AML_Qsub_Triplet_Permutated_112119.log
echo 'Rscript qsub_triplet_permutated.R TALL 2' | qsub -cwd -e TALL_Qsub_Triplet_Permutated_91019.error -o TALL_Qsub_Triplet_Permutated_112119.log
echo 'Rscript qsub_triplet_permutated.R BALL 3' | qsub -cwd -e AML_Qsub_Triplet_Permutated_91019.error -o BALL_Qsub_Triplet_Permutated_112119.log
echo 'Rscript qsub_triplet_permutated.R RT 4' | qsub -cwd -e RT_Qsub_Triplet_Permutated_91019.error -o RT_Qsub_Triplet_Permutated_112119.log
echo 'Rscript qsub_triplet_permutated.R WT 5' | qsub -cwd -e WT_Qsub_Triplet_Permutated_91019.error -o WT_Qsub_Triplet_Permutated_112119.log
echo 'Rscript qsub_triplet_permutated.R NBL 6' | qsub -cwd -e NBL_Qsub_Triplet_Permutated_91019.error -o NBL_Qsub_Triplet_Permutated_112119.log

Combines the above results - filters out non significant triplets. padj < 0.1.
echo 'Rscript combine_triplet_permute.R AML' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e AML_Combined_Permute.error -o AML_Combined_Permute.log
echo 'Rscript combine_triplet_permute.R BALL' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e BALL_Combined_Permute.error -o BALL_Combined_Permute.log
echo 'Rscript combine_triplet_permute.R TALL' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e TALL_Combined_Permute.error -o TALL_Combined_Permute.log
echo 'Rscript combine_triplet_permute.R NBL' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e NBL_Combined_Permute.error -o NBL_Combined_Permute.log
echo 'Rscript combine_triplet_permute.R RT' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e RT_Combined_Permute.error -o RT_Combined_Permute.log
echo 'Rscript combine_triplet_permute.R WT' | qsub -cwd -l mem_free=15G -l m_mem_free=15G -l h_vmem=15G -e WT_Combined_Permute.error -o WT_Combined_Permute.log


