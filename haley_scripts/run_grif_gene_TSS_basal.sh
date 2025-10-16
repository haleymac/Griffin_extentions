# Run griff_gene_TSS.py in your current working directory - copy this and run it wherever you want your griffin output 
# Will create griffin results in this working directory 
# Arguments:
#       --genes A list of genes you are interested in and want to see coverage over. Script will run griffin on each gene individually, and also the genes as a merged set 
#       --geneset_name Name to call the combined geneset ex. 'Housekeeping_genes'
# Note that your samples.yaml file in the nucleosome profiling snakemake directory will already need to be set up 

# Removed VGLL1 because griffin chokes on it

python /projects/pangen/analysis/hmac/applications/Griffin/haley_scripts/griff_gene_TSS.py \
  --genes SCEL FAM83A KRT15 DHRS9 GPR87 ANXA8L1 S100A2 LY6D \
         TNS4 KRT6A KRT6C SLC2A1 KRT17 SERPINB4 SERPINB3 \
         SPRR1B SPRR3 FGFBP1 LEMD1 CST6 KRT7 CTSV \
         AREG UCA1 \
  --geneset_name basal

