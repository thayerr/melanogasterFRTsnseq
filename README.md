# melanogasterFRTsnseq
# Analysis of single nuclei RNA sequencing of the Drosophila melanogaster female somatic reproductive tract (uterus, glands, sperm storage organs), in their unmated state

Please see the corresponding publication: Thayer et al 2024 "Regional specialization, polyploidy, and seminal fluid transcripts in the Drosophila female reproductive tract" 

Major objectives of this analytical pipeline include: a) to identify cell types of the five female somatic reproductive organs and their transcriptomes. b) to identify marker genes, subsequently used to annotate the cell types. c) use the cell-level transcriptomes to describe biological attributes of the cell types.

Files included and their relationships:  
1- FRT_unmated_sn_reduction_GitHubv.Rmd This is code used to inspect data quality and filter and clean the data, including by running SoupX to remove background reads and DoubletFinder to identify and remove likely doublets. Also includes initial clustering.  
2- R code used for analysis, including defining clusters and exploring cluster attributes   
3- average_per_cluster.py A custom python script used to created the Dataset S3 and Dataset S4 supplemental data files accompanying the publication. These files report the average gene expression by cell type, in terms of logCPM and percent expressing cells.  
4- fully analyzed, interactable Seurat object. This can be used to easily explore the dataset, generate visualizations, check expression levels for genes of interest, and so on.   
