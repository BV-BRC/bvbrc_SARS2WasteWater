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
    df_samples['assembly - complete'] = "Incomplete"
    df_samples['Freyja - Analysis'] = "Incomplete"
    df_samples['Freyja - Visualization'] = "Incomplete"

    # # Step 3: Check the existence of output files and update the dataframe
    for sample_id, idx in zip(df_samples['sampleID'], df_samples.index):
        if os.path.getsize(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv') != "Incomplete":
            df_samples['assembly - complete'][idx] = "Complete"
        ## parse the assembly stats file ##
            assembly_stats = parse_assembly_stats_file(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv')
        # Every line in the assembly stats becomes a new column with the parsed value for that sample in the row
            for key, value in assembly_stats.items():
                df_samples.at[idx, key] = value
        if os.path.exists(f'output/{sample_id}/{sample_id}_freyja_result.tsv'):
            df_samples['Freyja - Analysis'][idx] = "Complete"
        if os.path.exists(f'output/{sample_id}/{sample_id}_lineage_plot.png'):
            df_samples['Freyja - Visualization'][idx] = "Complete"
        # clean up primer path only show the bed file and remove the path
        df_samples['primers'][idx] = get_last_segment(df_samples['primers'][idx])
        # TO DO:parse samtools flagstats file
        # parse_samtools_flagstats(f'output/{sample_id}/assembly/{sample_id}_flagstat.txt')

    # the reference col shows an internal path - dropping because it does not help user
    df_samples = df_samples.drop(['reference'], axis='columns')

    # reorder columns

    df_samples = df_samples[[
    "sampleID",
    "assembly - complete",
    "Freyja - Analysis",
    "Freyja - Visualization",
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
