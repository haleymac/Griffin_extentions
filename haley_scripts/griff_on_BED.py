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
# FILENAME : griff_on_BED.py
#
# AUTHOR : haleymac
#
# DATE : April 28, 2025
#
# DESCRIPTION :
#
#   -  Run griffin on sites defined in a BED file- looking at coverage over the intervals in the BED as one 'site'
#   -  Uses samples in /projects/pangen/analysis/hmac/applications/Griffin/snakemakes/griffin_nucleosome_profiling/config/samples.GC.yaml
#       - You will need to manually update if you don't want to just see PAN21 and PAN27
#   - Takes CLI arguments:
#   - Outputs sites_files, results in working directory script is called in 
# ------------------------------------------------------------------------------



import pandas as pd

def add_strand_from_start_end(bed3_path, output_path=None):
    """
    Reads a BED3 file, adds a 'Strand' column:
      - '+' if start < end
      - '-' if start > end
    Optionally writes to a new file.
    """
    # Load BED3
    bed = pd.read_csv(bed3_path, sep='\t', header=None, names=['Chrom', 'Start', 'End'])

    # Add Strand based on Start and End
    bed['Strand'] = bed.apply(lambda row: '+' if row['Start'] < row['End'] else '-', axis=1)

    # Optionally write to file
    if output_path:
        bed.to_csv(output_path, sep='\t', index=False, header=True)

    return bed





