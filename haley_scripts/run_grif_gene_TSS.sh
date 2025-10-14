# Run griff_gene_TSS.py in your current working directory - copy this and run it wherever you want your griffin output 
# Will create griffin results in this working directory 
# Arguments:
#       --genes A list of genes you are interested in and want to see coverage over. Script will run griffin on each gene individually, and also the genes as a merged set 
#       --geneset_name Name to call the combined geneset ex. 'Housekeeping_genes'

python /projects/pangen/analysis/hmac/applications/Griffin/haley_scripts/griff_gene_TSS.py \
  --genes MT-CO1 MT-ND4 MT-ATP6 MT-CO2 MT-RNR2 \
  --geneset_name PAN21_top5
