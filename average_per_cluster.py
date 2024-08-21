#!/usr/bin/env python

# usage: average_per_cluster.py

# Written by Rachel Thayer August 2023, updated August 2024

# Purpose: parse the seurat output data files to give two summary tables of gene expression by cell type. Specifically, one table gives the mean logCPM per gene per cluster. The second table gives the percent of expressing cells per gene per cluster. 

# Inputs: This requires counts files made by r/Seurat. FRT_unmated_sn_analysis_final.R

# Runtime on the order of 28 hours

# Usage notes! 
# you must update the hard-coded file paths to sync with your own directory layout
# This program writes the output file's header to a different file; you need to concatenate the header and content files together afterward in terminal. 
		# cat header_per_cluster_renamed.txt percent_expressing_per_cluster_renamed.txt > percent_expressing_per_cluster_renamed_colnames.txt


import statistics

# load data from Seurat: read counts; droplet barcode IDs and the cell types to which each barcode has been assigned
in_counts = '/data/rthayer/FRT_unmated_nuclei/seurat_data.csv' #this is from the "data" slot in seurat, which is counts-like data, back-transformed, and also log-normalized. Values are always > 0. It's the same data normalization used within seurat by the 'find markers' command
in_cell_types = '/data/rthayer/FRT_unmated_nuclei/seurat_metadata_renamed.csv'


outfile_cpm = '/data/rthayer/FRT_unmated_nuclei/logCPM_per_cluster_renamed.txt'
outfile_perct = '/data/rthayer/FRT_unmated_nuclei/percent_expressing_per_cluster_renamed.txt'
outfile_header = '/data/rthayer/FRT_unmated_nuclei/header_per_cluster_renamed.txt'

# define names of clusters (i.e. putative cell types). These must match the file format in in_cell_types.
clusters = ["00-OV-1","01-OV-2","02-OV-3","03-SR-prox1","04-SR-prox2","05-SR-prox3","06-SR-dist1","07-SR-dist2","08-SSC","10-ST-ep","11-ST-PV-duct","12-PSC","14-UT-to-ST","15-UT-PE","16-UT-to-SR","17-UT-ant1","18-UT-ant2","19-UT-mid1","20-UT-mid2","21-UT-mid3","22-UT-mid4","23-UT-SVI","24-UT-post","25-FB","28-MUS1","29-MUS2","30-MUS3","31-MUS4","32-MUS5","33-hemocyte","34-glia","35-trachea","36-neuron","37-oenocyte","38-artefact","39-artefact-2","40-artefact-3","41-SSC-like","42-PSC-like"]


#For each cluster: get its barcodes, then get gene counts from its cells
# this piece writes out a file with rows per cluster and columns per gene; header column of gene IDs is regrettably the last line of the file

gene_list = ['ClusterID'] #initialize a holder for the parsed data
for numbe in clusters:               
    current_cluster_nucs =[]
    outfile_cpm_lines = [numbe]
    outfile_perct_lines =[numbe]
    with open(in_cell_types, 'r') as fsr:
        next(fsr)
        for line in fsr:
            seur_cls = line.split(',')[26].strip() # column index 26 has the cell type name assigned
            if seur_cls == numbe:
                barcode = line.split(',')[0][0:20] #the second indexing [0:20]removes a newline character
                current_cluster_nucs.append(barcode) # build a list of all barcodes that belong to the current cell type
    fsr.close()
    keep_columns=[]
    #get gene counts for the current cluster in the loop
    with open(in_counts, 'r') as f2:
        header = f2.readline()
        for n,barcode in enumerate(header.split(',')):
            if barcode in current_cluster_nucs:
                keep_columns.append(n)
        for line in f2: #each line is count data for 1 gene
            line_holder=[]
            for n in keep_columns: #get the column numbers that correspond to barcodes that belong to the current cell type. This allows us to keep only read counts from those columns in the next lines
                line_holder.append(line.split(',')[n+1]) #data lines are shifted one column longer owing to the gene name at the beginning of each
            gene=line.split(',')[0]
            if numbe == '37-oenocyte': #save the order of gene names for the outfile, but only need to do it once, rather than for every cycle of the loop (i.e. re-doing it for every cell type.) Doesn't matter which cluster to write on; I've just picked oenocytes randomly.
                gene_list.append(gene)
            newvector=','.join(str(x) for x in line_holder)
            workable = [float(i) for i in line_holder[1:]] # convert counts values to float so that they can be averaged
            avg_counts = round(statistics.mean(workable), 4) # average expression (logCPM) for the gene for the current cell type, calculated to 4 decimal places
            outfile_cpm_lines.append(avg_counts)
            per_greater_one = round(sum(1 for i in workable if i >1)/ len(workable)*100, 4) # calculate for what percent of the cells in the current cell type, the current gene is expressed at a level > logCPM of 1. 
            #per_non_zero = sum(1 for i in workable if i >0)/ len(workable) #Alternatively, calculate percent of cells that express the gene above logCPM of 0. This seldom differs much from % greater expression than 1
            outfile_perct_lines.append(per_greater_one)
    f2.close()
    fh3=open(outfile_cpm, "a")
    writeme='\t'.join([str(i) for i in outfile_cpm_lines])
    fh3.write(writeme)
    fh3.write('\n')
    fh3.close()
    fh4=open(outfile_perct, "a")
    writemeB='\t'.join([str(i) for i in outfile_perct_lines])
    fh4.write(writemeB)
    fh4.write('\n')
    fh4.close()
    
header_out = '\t'.join(gene_list)
header_out = header_out + '\n'
fh5=open(outfile_header, "w")
fh5.write(header_out)
fh5.close()




