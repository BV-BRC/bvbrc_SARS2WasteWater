import click
import glob
import pandas as pd
import os


def get_last_segment(path):
    return path.rsplit('/', 1)[-1]

def parse_assembly_stats_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Create a dictionary to store parsed values
    parsed_values = {}
    # Iterate through each line and extract key-value pairs
    for line in lines:
        key, value = line.strip().split('\t')
        parsed_values[key] = value
    return parsed_values

def complile_stats(path):
    # Step 1: Get the sample IDs
    # set up a structure like the snakemake wildcards
    # make a list of the input files 
    # single end unzipped
    uz_se_file_list = glob.glob("staging/se_reads/*.fastq")
    # single end zipped
    z_se_file_list = glob.glob("staging/se_reads/*.fastq.gz")
    # paired end unzipped
    uz_pe_file_list = glob.glob("staging/pe_reads/*_R1.fastq")
    # paired end zipped
    z_pe_file_list = glob.glob("staging/pe_reads/*_R1.fastq.gz")
    # all together now
    file_list = uz_se_file_list + z_se_file_list + uz_pe_file_list + z_pe_file_list


    # Extract the wildcard part from each file path
    wildcards = [os.path.basename(file_path).split('.')[0] for file_path in file_list]

    modified_wildcards = []
    # For each item in the input file list, remove read 1 so the sample ID is clean
    for s in wildcards:
        # Check if the string contains "_R1"
        if '_R1' in s:
            # If it does, remove "_R1"
            modified_string = s.replace('_R1', '')
        else:
            modified_string = s
        
        # Append the modified string (clean sample ID) to the new list
        modified_wildcards.append(modified_string)

    # make a dataframe with one row for each clean sample id
    df_samples = pd.DataFrame(modified_wildcards, columns=['sampleID'])

    # make a column for each file that we are tracking
    df_samples['assembly - sorted bam']=0
    df_samples['assembly - complete'] = 0
    df_samples['freyja - variants_file'] = 0
    df_samples['freyja - depth'] = 0
    df_samples['freyja - sample_demixing_results'] = 0
    df_samples['freyja - sample_variant_plot'] = 0
    df_samples['freyja - sample_lineage_plot'] = 0

    # # Step 3: Check the existence of output files and update the dataframe
    for sample_id, idx in zip(df_samples['sampleID'], df_samples.index):
        if os.path.exists(f'output/{sample_id}/assembly/{sample_id}.sorted.bam'):
            df_samples['assembly - sorted bam'][idx] = 1
        if os.path.getsize(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv') != 0:
            df_samples['assembly - complete'][idx] = 1
        ## parse the assembly stats file ##
            assembly_stats = parse_assembly_stats_file(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv')
        # Every line in the assembly stats becomes a new column with the parsed value for that sample in the row
            for key, value in assembly_stats.items():
                df_samples.at[idx, key] = value
        ## sample level ##
        # Freyja variants
        if os.path.exists(f'output/{sample_id}/{sample_id}_freyja_variants.tsv'):
            df_samples['freyja - variants_file'][idx] = 1
        # Freyja depth
        if os.path.exists(f'output/{sample_id}/{sample_id}_freyja.depths'):
            df_samples['freyja - depth'][idx] = 1
        # Freyja demix file created
        if os.path.exists(f'output/{sample_id}/{sample_id}_demixing_result.csv'):
            df_samples['freyja - sample_demixing_results'][idx] = 1
        # sample variant level plot 
        if os.path.exists(f'output/{sample_id}/{sample_id}_variant_plot.png'):
            df_samples['freyja - sample_variant_plot'][idx] = 1
        if os.path.exists(f'output/{sample_id}/{sample_id}_lineage_plot.png'):
            df_samples['freyja - sample_lineage_plot'][idx] = 1
        # clean up primer path only show the bed file and remove the path
        df_samples['primers'][idx] = get_last_segment(df_samples['primers'][idx])
        # TO DO:parse samtools flagstats file
        # parse_samtools_flagstats(f'output/{sample_id}/assembly/{sample_id}_flagstat.txt')
    # run level ##
    df_samples['multisample - demix_results'] = 0
    df_samples['multisample - summarized_variants'] = 0
    df_samples['multisample - lineage_plot'] = 0

    if  os.path.exists("output/multisample_demix_results.tsv"):
        df_samples['multisample - demix_results']= 1
    if  os.path.exists("output/multisample_summarized_variants_plot.png"):
        df_samples['multisample - summarized_variants']= 1
    if  os.path.exists("output/multisample_lineage_plot.png"):
        df_samples['multisample - lineage_plot']= 1

    # the reference col shows an internal path - dropping because it does not help user
    df_samples = df_samples.drop(['reference'], axis='columns')

    # reorder columns

    df_samples = df_samples[["sampleID",
    "assembly - sorted bam",
    "assembly - complete",
    "freyja - variants_file",
    "freyja - depth",
    "freyja - sample_demixing_results",
    "freyja - sample_variant_plot",
    "freyja - sample_lineage_plot",
    "multisample - demix_results",
    "multisample - summarized_variants",
    "multisample - lineage_plot",
    "depth_mean",
    "depth_median",
    "depth_stdv",
    "depth_min",
    "depth_max",
    "n_count",
    "n_blocks",
    "fasta_length",
    "primers",
    "primer_count",
    "mapped_reads",
    "unmapped_reads",
    "primer_trim_count",
    "primer_trim_pct",
    "variant_count"]]


    # write out to CSV - without the index column
    df_samples.to_csv(path, index = False)


@click.command()
@click.argument("path")
def cli(path):
    complile_stats(path)

if __name__ == '__main__':
    cli()
