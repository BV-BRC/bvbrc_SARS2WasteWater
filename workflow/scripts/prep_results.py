import click
import json
import pandas as pd
import os


def get_last_segment(path):
    return path.split('/')[-1]

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
    # Step 1: get the sample ids from the config file 
    with open('config.json', 'r') as file:
        data = json.load(file)

        # Extracting sample IDs from the "paired_end_libs" and "single_end_libs" lists
        paired_end_sample_ids = [lib['sample_id'] for lib in data['params']['paired_end_libs']]
        single_end_sample_ids = [lib['sample_id'] for lib in data['params']['single_end_libs']]
        srr_sample_ids = [lib['sample_id'] for lib in data['params']['srr_libs']]

        # Combine all sample IDs into one list
        all_sample_ids = paired_end_sample_ids + single_end_sample_ids + srr_sample_ids

        unique_sample_ids = list(set(all_sample_ids))


    # make a dataframe with one row for each clean sample id
    df_samples = pd.DataFrame(unique_sample_ids, columns=['sampleID'])
    df_samples['Assembly'] = "Incomplete"
    df_samples['Freyja - Analysis'] = "Incomplete"
    df_samples['Freyja - Visualization'] = "Incomplete"

    # # Step 2: Check the existence of output files and update the dataframe
    for sample_id, idx in zip(df_samples['sampleID'], df_samples.index):
        if os.path.getsize(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv') != "Incomplete":
            df_samples['Assembly'][idx] = "Complete"
        ## parse the assembly stats file ##
            assembly_stats = parse_assembly_stats_file(f'output/{sample_id}/assembly/{sample_id}.statistics.tsv')
        # Every line in the assembly stats becomes a new column with the parsed value for that sample in the row
            for key, value in assembly_stats.items():
                df_samples.at[idx, key] = value
        else:
            df_samples['Assembly'] = "Incomplete"
        if os.path.exists(f'output/{sample_id}/freyja/{sample_id}_freyja_result.tsv'):
            df_samples['Freyja - Analysis'][idx] = "Complete"
        if os.path.exists(f'plots/lineages_plot.html'):
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
    "Assembly",
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

    # rename columns
    df_samples.rename(columns = {"sampleID":"Sample ID", "depth_mean":"Depth Mean", \
                                "depth_median":"Depth Median", "depth_stdv": "Depth Standard Deviation", \
                                "depth_min":"Depth Minimum", "depth_max":"Depth Maximum", \
                                "n_count":"Total N Count", "n_blocks":"N Blocks", \
                                "fasta_length":"Fasta Length", "primers":"Primers", \
                                "primer_count":"Primer Count","mapped_reads":"Mapped Reads",
                                "unmapped_reads":"Unmapped Reads","primer_trim_count":"Primer Trim Count", \
                                "primer_trim_pct":"Percentage of Primers Trimmed","variant_count":"Variant Count" }, inplace = True)
    # write out to CSV - without the index column
    df_samples.to_csv(path, index = False, sep="\t")


@click.command()
@click.argument("path")
def cli(path):
    complile_stats(path)

if __name__ == '__main__':
    cli()