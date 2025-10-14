
import os
import sys
import argparse
import pandas as pd
import pysam
import pybedtools
import pyBigWig
import numpy as np
import time
import yaml 
from multiprocessing import Pool



sample_name = "PAN21_WBC"
bam_path = "/projects/pangen/analysis/jtopham/fragx/samples/PAN21/PanGen-PAN21-WBC.bam"
GC_bias_path = "/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_GC_and_mappability_correction/results/GC_bias/PAN21_WBC.GC_bias.txt"
mappability_bias_path = "none"
tmp_dir = "tmp"
ref_seq_path = "/projects/pangen/analysis/jtopham/fragx/fragx_hg38_ref.fa"
mappability_bw = "/projects/pangen/analysis/hmac/applications/Griffin/Ref/k100.Umap.MultiTrackMappability.bw"
chrom_sizes_path = "/projects/pangen/analysis/hmac/applications/Griffin/Ref/hg38.standard.chrom.sizes"
sites_yaml = "/projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config/sites.yaml"
griffin_scripts_dir = "/projects/pangen/analysis/hmac/applications/Griffin/scripts"
chrom_column = "Chrom"
position_column = "position"
strand_column = "Strand"
chroms =  ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
norm_window =  [-5000, 5000]
sz_range = [100, 200]
map_q = 20
number_of_sites = "none"
sort_by = "none"
ascending = "none"
CPU = 8


mappability_correction = "False"


##sample specific params for testing
# sample_name = 'MBC_1041_1_ULP'
# bam_path = '../../../../griffin_revisions_1/MBC_copy_bams/bam_file_copies/MBC_1041_1_ULP_recalibrated.bam'
# GC_bias_path = '../../../../griffin_revisions_1/GC_correction/MBC_ULP_GC_and_mappability_correction/results/GC_bias/MBC_1041_1_ULP.GC_bias.txt'

# mappability_bias_path = '../../../../griffin_revisions_1/GC_correction/MBC_ULP_GC_and_mappability_correction/results/mappability_bias/MBC_1041_1_ULP.mappability_bias.txt'
# mappability_correction = 'True'

# # mappability_bias_path = 'none'
# # mappability_correction = 'False'

# tmp_dir = 'tmp'

# ref_seq_path = '/fh/fast/ha_g/grp/reference/GRCh38/GRCh38.fa'
# mappability_bw='../../../../griffin_revisions_1/genome/k100.Umap.MultiTrackMappability.hg38.bw'
# chrom_sizes_path = '/fh/fast/ha_g/grp/reference/GRCh38/hg38.standard.chrom.sizes'

# # #additional params for testing
# sites_yaml = '/fh/fast/ha_g/user/adoebley/projects/griffin_revisions_1/MBC/CNA_correction_100kb_ATAC_np/config/sites.yaml'
# griffin_scripts_dir = '../'

# chrom_column = 'Chrom'
# position_column = 'position'
# strand_column = 'Strand'
# chroms = ['chr'+str(m) for m in np.arange(1,23)]

# norm_window = [-5000, 5000] 
# # norm_window = [-5000, 5000] 
# sz_range = [100, 200]
# map_q = 20

# number_of_sites = 500
# sort_by = 'Chrom'
# #sort_by = 'peak.count'
# ascending = 'False'

# # number_of_sites = 'none'
# # sort_by = 'none'
# # ascending = 'none'

# CPU = 6



#define global parameters and open global files
########################################
#GET GC BIAS
########################################
#open the GC_bias file 
GC_bias = pd.read_csv(GC_bias_path, sep='\t')
#get rid of extremely low GC bias values
#these fragments will now be excluded 
#these fragments are extremely rare so it is difficult to get a good estimate of GC bias
GC_bias['smoothed_GC_bias'] = np.where(GC_bias['smoothed_GC_bias']<0.05,np.nan,GC_bias['smoothed_GC_bias'])

GC_bias = GC_bias[['length','num_GC','smoothed_GC_bias']]
GC_bias = GC_bias.set_index(['num_GC','length']).unstack()

#convert to a dictionary
GC_bias = GC_bias.to_dict()

#get rid of values where the num_GC is greater than the length (included due to the way I made the dict)
GC_bias2 = {}
for key in GC_bias.keys():
    length = key[1]
    GC_bias2[length] = {}
    for num_GC in range(0,length+1):
        bias = GC_bias[key][num_GC]
        GC_bias2[length][num_GC]=bias
GC_bias = GC_bias2 
del(GC_bias2)

#snakemake should create these folders, but if not using the snakemake, this is needed
tmp_sample_dir = tmp_dir+'/'+sample_name
if not os.path.exists(tmp_sample_dir): 
    os.mkdir(tmp_sample_dir)

tmp_pybedtools = tmp_sample_dir+'/tmp_pybedtools'
if not os.path.exists(tmp_pybedtools): 
    os.mkdir(tmp_pybedtools)
pybedtools.set_tempdir(tmp_pybedtools)

tmp_bigWig = tmp_sample_dir+'/tmp_bigWig'
if not os.path.exists(tmp_bigWig): 
    os.mkdir(tmp_bigWig)
    
    
    
    
def import_and_filter_sites(site_name,site_file,strand_column,chrom_column,position_column,chroms,ascending,sort_by,number_of_sites):
    import pandas as pd
    current_sites = pd.read_csv(site_file,sep='\t')
    if strand_column not in current_sites.columns:
        current_sites[strand_column]=0 
            
    #throw out sites that aren't on the selected chroms
    current_sites = current_sites[current_sites[chrom_column].isin(chroms)]
    
    #select the sites to use if specified
    if sort_by.lower()=='none': #if using all sites 
        print(site_name,'processing all '+str(len(current_sites))+' sites')
    
    else: #othewise sort by the specified column
        current_sites=current_sites.sort_values(by=sort_by,ascending=ascending).reset_index(drop=True)#sort and reset index
        current_sites=current_sites.iloc[0:int(number_of_sites)]
        print(site_name+'\tprocessing',len(current_sites),'sites\t('+
              str(sort_by),'range after sorting: ',min(current_sites[sort_by]),'to',
              str(max(current_sites[sort_by]))+')')

    current_sites = current_sites[[chrom_column,position_column,strand_column]]
    current_sites['site_name']=site_name
    return(current_sites)


#import the site_lists
with open(sites_yaml,'r') as f:
    sites = yaml.safe_load(f)
sites = sites['site_lists']


print(sites)

all_sites = pd.DataFrame()
for site_name in sites.keys():
    site_file = sites[site_name]
    current_sites = import_and_filter_sites(site_name,site_file,strand_column,chrom_column,position_column,chroms,ascending,sort_by,number_of_sites)
    all_sites = all_sites.append(current_sites, ignore_index=True).copy()
    sys.stdout.flush()




def define_fetch_interval(name_to_print,sites,chrom_column,position_column,chroms,chrom_sizes_path,upstream_bp,downstream_bp):
    import pandas as pd
    import numpy as np
    #separate fw and reverse sites
    fw_markers = ['+',1,'1']
    rv_markers = ['-',-1,'-1']
    fw_sites = sites[sites['Strand'].isin(fw_markers)].copy()
    rv_sites = sites[sites['Strand'].isin(rv_markers)].copy()

    undirected_sites = sites[~(sites['Strand'].isin(fw_markers+rv_markers))].copy()

    if len(rv_sites)+len(fw_sites)+len(undirected_sites)==len(sites):
        print(name_to_print+' (fw/rv/undirected/total): '+
              str(len(fw_sites))+'/'+
              str(len(rv_sites))+'/'+
              str(len(undirected_sites))+'/'+
              str(len(sites)))
    else: #I don't think this should ever happen...
        print('total fw sites:\t\t'+str(len(fw_sites)))
        print('total rv sites:\t\t'+str(len(rv_sites)))
        print('total undirected sites:'+'\t'+str(len(undirected_sites)))
        print('total sites:\t\t'+str(len(sites)))
        sys.exit('Problem with strand column')

    #set up to fetch a window extending across the desired window
    fw_sites['fetch_start'] = fw_sites[position_column]+upstream_bp
    fw_sites['fetch_end'] = fw_sites[position_column]+downstream_bp

    undirected_sites['fetch_start'] = undirected_sites[position_column]+upstream_bp
    undirected_sites['fetch_end'] = undirected_sites[position_column]+downstream_bp

    #for reverse sites, flip the window
    rv_sites['fetch_start'] = rv_sites[position_column]-downstream_bp
    rv_sites['fetch_end'] = rv_sites[position_column]-upstream_bp
    
    #merge fw and reverse back together and sort them back into the original order
    sites = fw_sites.append(rv_sites).append(undirected_sites).sort_index()
    sites = sites.sort_values(by = [chrom_column,position_column]).reset_index(drop=True)

    chrom_sizes = pd.read_csv(chrom_sizes_path, sep='\t', header=None)
    chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
    chrom_sizes = chrom_sizes.set_index(0)

    adjusted_ends_df = pd.DataFrame()

    for chrom in chroms:
        length = chrom_sizes.loc[chrom][1]
        current = sites[sites[chrom_column]==chrom].copy()
        current['fetch_start'] = np.where(current['fetch_start']<0,0,current['fetch_start'])
        current['fetch_end'] = np.where(current['fetch_end']>length,length,current['fetch_end'])    
        adjusted_ends_df = adjusted_ends_df.append(current)
    adjusted_ends_df = adjusted_ends_df.sort_values(by = [chrom_column,position_column]).reset_index(drop=True)
    adjusted_ends_df = adjusted_ends_df.copy()

        # Flip fetch_start and fetch_end if fetch_start > fetch_end
    adjusted_ends_df['fetch_start'], adjusted_ends_df['fetch_end'] = np.where(
        adjusted_ends_df['fetch_start'] <= adjusted_ends_df['fetch_end'],
        (adjusted_ends_df['fetch_start'], adjusted_ends_df['fetch_end']),
        (adjusted_ends_df['fetch_end'], adjusted_ends_df['fetch_start'])
    )
    return(adjusted_ends_df)

#number of bp to fetch upstream and downstream of the site
upstream_bp = norm_window[0]-sz_range[0] #this should be negative
downstream_bp = norm_window[1]+sz_range[0] #this should be positive
all_sites = define_fetch_interval('Total sites',all_sites,chrom_column,position_column,
                                                    chroms,chrom_sizes_path,upstream_bp,downstream_bp)



start_time = time.time()
all_sites_bed = pybedtools.BedTool.from_dataframe(all_sites[[chrom_column,'fetch_start','fetch_end']])
all_sites_bed = all_sites_bed.sort()

print(all_sites_bed)

    
