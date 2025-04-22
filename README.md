

# Griffin Extension
Griffin is a  flexible framework for nucleosome profiling of cell-free DNA - see their github: [https://github.com/adoebley/Griffin]

It works pretty nicely but is fussy to set up so I wrote some code to make running it easier via command line argument.


## Description
In this repo I have @adoebley's original repository contents, as well as some additional scripts in haley_scripts that make setting up a griffin run easier. 



### Scripts in haley_scripts:

#### grif_gene_TSS.py

This script is designed to run griffin's nucleosome profiling on the TSSs of a given set of genes. All it needs to run is a list of genes.
It will: 
1. Create site_files.tsv's that griffin needs to run in the working directory
    - these will be the +- 1000bp from the TSS of your genes 
    - there will be one for each individual gene, and one for the complete gene set 
2. Create your results/ directory in the working directory 
3. Create a results.yaml in Griffin/snakemakes/griffin_nucleosome_profiling/config containing the path to your new results dir - './results/'
    - I modified the original config to not include the results output dir, and modified the snakefile to take this new yaml instead. Doing this allows us to re-generate it each time we run this workflow on a new set of genes. 
4. Create a new sites.yaml in Griffin/snakemakes/griffin_nucleosome_profiling/config with the TSS sites for your genes
5. Calls griffin_nucleosome_profiling.snakefile to run on a slurm cluster 
    - results from running this will be in your working directory results/ folder...

To run the script (assuming you have griffin set up as directed on their Github): 

1. make a working directory you want your results in: 
```bash
mkdir griffin_run
cd griffin_run
```

2. Run script
```bash 
python /path/to/Griffin/haley_scripts/griff_gene_TSS.py \
  --genes ACTB GAPDH RPLP0 \
  --geneset_name housekeeping
```
**CLI Arguments:**
       --genes A list of genes you are interested in and want to see coverage over. Script will run griffin on each gene individually, and also the genes as a merged set 
       --geneset_name Name to call the combined geneset ex. 'Housekeeping_genes'

