# ------------------------------------------------------------------------------
#                ______,---'__,---'
#          _,-'---_---__,---'
#    /_    (,  ---____',
#    /  ',,   `, ,-'
#   ;/)   ,',,_/,'
#   | /\   ,.'//\
#   `-` \ ,,'    `.
#        `',   ,-- `.
#        '/ / |      `,         _
#        //'',.\_    .\\      ,{==>-
#     __//   __;_`-  \ `;.__,;'
#   ((,--,) (((,------;  `--' jv
# ------------------------------------------------------------------------------
# FILENAME : griff_gene_TSS.py
#
# AUTHOR : haleymac
#
# DATE : April 22, 2025
#
# DESCRIPTION :
#
#   -  Run griffin on a given set of genes - looking at coverage over their TSSs
#   -  Uses samples in /projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config/samples.GC.yaml
#       - You will need to manually update if you don't want to just see PAN21 and PAN27
#   - Takes CLI arguments:
#       --genes A list of genes you are interested in and want to see coverage over. Script will run griffin on each gene individually, and also the genes as a merged set 
#       --geneset_name Name to call the combined geneset ex. 'Housekeeping_genes'
#   - Outputs sites_files, results in working directory script is called in 
# ------------------------------------------------------------------------------



import pandas as pd
import os
import yaml
from pathlib import Path
import subprocess
from yaml.representer import SafeRepresenter
import argparse




# -----------------------------
# CLI ARGUMENT PARSING
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Prepare and run Griffin nucleosome profiling on a given set of gene TSSs.")
    parser.add_argument("--genes", nargs="+", required=True, help="List of genes to look at TSS coverage over")
    parser.add_argument("--geneset_name", required=True, help="Name of the geneset (used to name merged site file)")
    return parser.parse_args()



# -----------------------------
# STEP 1: MAKE GENE INPUT TSVS
# -----------------------------
def make_gene_input_tsvs(gene_list, output_dir, outfile):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs('logs/cluster', exist_ok=True)
    gtf = pd.read_csv('/projects/pangen/analysis/hmac/reference/hg38_ref/genome.ucsc.regions.txt', sep='\t')
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    df = gtf[gtf['name2'].isin(gene_list) & gtf['chrom'].isin(chroms)]
    df = df[['chrom', 'txStart', 'name2', 'strand']]
    df.rename(columns={'chrom': 'Chrom', 'txStart': 'Start', 'name2': 'Gene', 'strand': 'Strand'}, inplace=True)
    df['End'] = df['Start'] + 1
    df['position'] = df['Start']
    df['uniprotId'] = 'unknown'
    column_order = ["Chrom", "Start", "End", "uniprotId", "Gene", "position", 'Strand']
    df = df[column_order].drop_duplicates()
    gene_counts = df["Gene"].value_counts()
    gene_counter = {}

    def rename_gene(gene):
        if gene_counts[gene] > 1:
            gene_counter[gene] = gene_counter.get(gene, 0) + 1
            return f"{gene}_{gene_counter[gene]}"
        else:
            return gene

    df["Gene"] = df["Gene"].apply(rename_gene)

    for _, row in df.iterrows():
        row_df = pd.DataFrame([row])
        row_df.to_csv(os.path.join(output_dir, f"{row['Gene']}.tsv"), sep='\t', index=False)

    df.to_csv(f"{output_dir}/{outfile}", sep="\t", index=False)
    return df

# -----------------------------
# STEP 2: WRITE sites.yaml in /projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling
# -----------------------------
def write_sites_yaml():
    tsvdir = Path(".") / "site_files"
    output_file = Path("/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config/sites.yaml")

    with output_file.open("w") as file:
        file.write("site_lists:\n")

        def write_paths(indir):
            for site in os.listdir(indir):
                name = site.split('.')[0]
                full_path = (indir / f"{name}.tsv").resolve()
                file.write(f"    {name}: {full_path}\n")

        write_paths(tsvdir)

# -----------------------------
# STEP 3: WRITE results.yaml
# -----------------------------

# Trick PyYAML into writing a plain scalar string
class PlainString(str):
    pass

def plain_string_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='')

yaml.add_representer(PlainString, plain_string_representer)

def write_griffin_config(output_dir):
    output_dir = Path(output_dir).resolve()
    cwd = Path.cwd()
    results_dir = cwd / "results"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    config = {
        'results_dir': str(results_dir),
    }

    with open(output_dir / 'results.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


# -----------------------------
# STEP 4: RUN Snakemake
# -----------------------------
def run_snakemake():
    cmd = [
        "snakemake", "-s", "/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/griffin_nucleosome_profiling.snakefile",
        "--latency-wait", "60",
        "--keep-going",
        "--rerun-incomplete",
        "--cluster-config", "/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config/cluster_slurm.yaml",
        "--cluster", "sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}",
        "-j", "40"
    ]
    subprocess.run(cmd, check=True)

# -----------------------------
# MAIN SCRIPT LOGIC
# -----------------------------
if __name__ == "__main__":
    args = parse_args()
    gene_list = args.genes
    geneset_name = args.geneset_name
    
    tsvdir = 'site_files'
    outfile = f'{geneset_name}.tsv'
    make_gene_input_tsvs(gene_list, tsvdir, outfile)
    write_sites_yaml()
    write_griffin_config('/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config')
    run_snakemake()
