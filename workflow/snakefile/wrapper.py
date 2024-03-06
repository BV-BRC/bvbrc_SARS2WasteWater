#!/opt/patric-common/runtime/bin/python3
import csv
import json
import os
import pandas as pd
import re
import sys
import shutil
import subprocess

#
# python wrapper for the SARS2Waterwater analysis pipeline
# is designed to be executed by App-SARS2Wastewater perl script
#

def check_for_error_msgs(raw_msg, file, out_msg):
    with open(file) as f:
        if raw_msg in line:
            write_to_job_failed_file(out_msg)
        else:
            pass

def check_input_fastqs(input_dir, filename):
    input_path = f"{input_dir}/{filename}"
    if os.path.isfile(input_path) == True:
        return input_path
    else:
        msg = f"Error FASTQ input path: {input_path} does not exisit'. \n check on {input_path} \n"
        sys.stderr.write(msg)
        sys.exit(1)
        pass


def check_sample_metadata_csv(input_info, config_file, input_dir, staging_metadata_file):
    if staging_metadata_file != "0":
        if os.path.isfile(staging_metadata_file) == True:
            metadata_df = pd.read_csv(staging_metadata_file)
            #edited_staging_metadata_file = ""
            # check if the date column and sample column are the same length
            if len(metadata_df["Read_file"]) != len(metadata_df["sample_collection_datetime"]):
                msg = f"Error {staging_metadata_file} is missing a filename or date'. \n check on {staging_metadata_file} \n"
                sys.stderr.write(msg)
            else:
                ### rewrite the file with <sample_id>.tsv as required by Freyja ###
                # create a new column 'Sample' with sample IDs - using 'Sample' for Freyja's plotting function
                metadata_df["Sample"] = metadata_df["Read_file"].apply(lambda x: input_info.get(x))
                # add the file ending to the sample_id for Freyja 
                metadata_df["Sample"] = metadata_df["Sample"].astype(str) + ".tsv"
                # reorder the columns
                metadata_df = metadata_df[["Sample", "Read_file", "sample_collection_datetime"]]
                # drop the additional 'Read_file' column no longer needed 
                metadata_df.drop(columns=["Read_file"], inplace=True)
                edited_sample_metadata_csv_path = os.path.join(input_dir, "formatted_sample-metadata.csv")
                metadata_df.to_csv(edited_sample_metadata_csv_path, index=False)
        else:
            msg = f"Error sample metadata CSV {staging_metadata_file} does not exisit'. \n check on {staging_metadata_file} \n"
            sys.stderr.write(msg)
            sys.exit(1)
    else:
        # assign edited metadata to empty string
        edited_sample_metadata_csv_path = ""
        msg = f"Metadata CSV not provided. \n"
        sys.stderr.write(msg)
    ## add edited_sample_metadata_csv value with path to config JSON ##
    with open(config_file, 'r') as file:
        config = json.load(file)
    # add updated metadata_csv_path
    config["edited_sample_metadata_csv"] = edited_sample_metadata_csv_path
    # write the updated data back to the file
    with open('config.json', 'w') as file:
        json.dump(config, file, indent=4)
    return 


def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict


def post_processing_check(all_sample_ids, output_dir):
    dict_samples = {}
    complete = []
    incomplete = []
    #missing_sample_plot
    # missing_demix
    for sample_name in all_sample_ids: 
        # missing demix results
        # check for missing demix results
        freyja_path = f"{output_dir}/{sample_name}/{sample_name}_demixing_result.csv"
        if os.path.isfile(freyja_path) == True:
            msg = f"{sample_name}"
            complete.append(msg)
        else:
            dict_samples[sample_name] = False
            msg = f"{sample_name}"
            incomplete.append(msg)
        sample_plot = 
    ## write the complete and incomplete samples to the error report for user to see ##
    if len(incomplete) != 0:
        msg = f"Freyja results produced for the following samples: {complete}. \n \
                Freyja results missing: Review the following samples: **{incomplete}** \n \n \
                If you have any questions about this job please go to the jobs page \n \
                select this job. Then on the green action bar, select 'report issue'  \n \
                A team member will assist you with the analysis. \n"
        write_to_job_failed_file(msg)
        # check for errors then add them to Job Failed text file 
        for sample in incomplete:
            ### Check demix output ###
            raw_msg = "demix: Solver error encountered, most likely due to insufficient sequencing depth."
            out_msg = f"The Freyja demix command gave this error for {sample}: \n 'demix: Solver error encountered, \
                most likely due to insufficient sequencing depth.'"
            demix_error_file = f"tmp/{sample}_demix_error.txt"
            check_for_error_msgs(raw_msg, demix_error_file, out_msg)
            ### Check sample plots ##
            # variant plot
            raw_msg = "ERROR: No samples matching coverage requirements, so no plot will be generated."
            out_msg = "The Freyja plot command gave this error for {sample}: \n 'ERROR: No samples matching \
                coverage requirements, so no plot will be generated.' \
                The default coverage requirement is 60% \n"
            sample_plot_error_variant_file = f"tmp/{sample}_sample_variant_plot_error.txt"
            check_for_error_msgs(raw_msg, sample_plot_error_variant_file, out_msg)
            # linage plot
            raw_msg = "ERROR: No samples matching coverage requirements, so no plot will be generated."
            out_msg = "The Freyja plot command with flag '--lineage' gave this error for {sample}: \n 'ERROR: No samples matching \
                coverage requirements, so no plot will be generated.' \
                The default coverage requirement is 60% \n"
            sample_plot_error_lineage_file = f"tmp/{sample}_sample_lineage_plot_error.txt"
            # lineage plot coverage depth error
            # lineage plot too many lineages error
            raw_msg = "Exception: Too many lineages to show. Use --config to group."
            out_msg = "The Freyja plot command with flag '--lineage' gave this error for {sample}: \n \
                'Exception: Too many lineages to show.'"
            check_for_error_msgs(raw_msg, sample_plot_error_lineage_file, out_msg)

        sys.exit(1)
    else:
        msg = f"Freyja results produced for the following samples: {complete}. \n"
        sys.stderr.write(msg)
        return True


def preprocessing_check(output_dir, input_dict):
    ## Using the input dictionary instead of path from staging directory due to file name confusion
    # Parse "sample_id" from "paired_end_libs"
    paired_sample_ids = [item["sample_id"] for item in input_dict.get("paired_end_libs", [])]

    # Parse "sample_id" from "single_end_libs"
    single_sample_ids = [item["sample_id"] for item in input_dict.get("single_end_libs", [])]

    # Parse "sample_id" from "srr_libs"
    srr_sample_ids = [item["sample_id"] for item in input_dict.get("srr_libs", [])]
    # Merge all sample IDs into one list
    all_sample_ids = paired_sample_ids + single_sample_ids + srr_sample_ids

    # Edit the sample ids to match the sample ids defined in set-up-sample-dictionary(input_dir input_dict output_dir cores)
    clean_sample_ids = []
    for sample_id in all_sample_ids:
        # Define a regular expression pattern to match all special characters except underscore
        pattern = r"[^a-zA-Z0-9_]"
        # Use the re.sub() function to replace all matches of the pattern with an empty string
        sample_id = re.sub(pattern, "", sample_id)
        clean_sample_ids.append(sample_id)
    all_sample_ids = clean_sample_ids
        
    # file check for preprocessing/sars-one-codex
    dict_samples = {}
    complete = []
    incomplete = []
    wrong_sequence_type = []


    # check if iVar bam exisits for each file 
    for sample_name in all_sample_ids:
        iVar_bam = f"{output_dir}/{sample_name}/assembly/{sample_name}.ivar.bam"
        stats_path = f"{output_dir}/{sample_name}/assembly/{sample_name}.statistics.tsv"

        ### check for incorrect sequencing type ###
        # if the wrong sequencing type is given the stats file will be empty
        if os.path.isfile(stats_path) == True:
            if os.path.getsize(stats_path) == 0:
                wrong_sequence_type.append(sample_name)
            else:
                msg = f"assembly complete for {sample_name}"
                sys.stderr.write(msg)
        else:
            wrong_sequence_type.append(sample_name)
        ### check for assembly results ###
        if os.path.isfile(iVar_bam) == True:
            complete.append(sample_name)
        else:
            dict_samples[sample_name] = False
            incomplete.append(sample_name)

    ### Write out messages and warnings to JobFailed.txt ###
    # if both wrong sequencing type and empty ivar bam file
    if len(incomplete) != 0 and len(wrong_sequence_type) !=0:
            msg = f"WARNING: Errors in FASTQ proccessing \n \
            \n CHECK SEQUENCING TYPE: This service only accepts **AMPLICON** based sequencing. \n \
            Please check the sequencing type for the following samples: **{wrong_sequence_type}** \n \n"
            write_to_job_failed_file(msg)
            msg = f" WARNING: CHECK SAMPLE(S) \n \
                Primer trimming and variant calling is complete for the following samples: {complete}. \n \
                CHECK SAMPLE: Primer trimming and variant calling is INCOMPLETE for the following samples: **{incomplete}** \n \
                Please review the sample(s) by uploading the FASTQ file(s) to the FastqUtils service.\n \
                Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing \n"
            write_to_job_failed_file(msg)
    ## if an iVAR trimmed bam file does not exist but statstics.tsv exisits 
    elif len(incomplete) != 0:
        msg = f"WARNING: CHECK SAMPLE(S) \n \
                Primer trimming and variant calling is complete for the following samples: {complete}. \n \
                CHECK SAMPLE: Primer trimming and variant calling is INCOMPLETE for the following samples: **{incomplete}** \n \
                Please review the sample(s) by uploading the FASTQ file(s) to the FastqUtils service.\n \
                Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing \n"
        write_to_job_failed_file(msg)
        sys.exit(1)
    ## if the statistics file is empty it is probably the wrong sequencing type write inormative message and stop the analysis 
    elif len(wrong_sequence_type) !=0:
        msg = f"WARNING: Errors in FASTQ proccessing \n \
                CHECK SEQUENCING TYPE: \n \
                \n This service only accepts amplicon based sequencing. \n \
                Please check the sequencing type for the following samples: **{wrong_sequence_type}** \n \
                Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing"
        write_to_job_failed_file(msg)
        sys.exit(1)
    else:
        msg = "preprocessing complete \n"
        sys.stderr.write(msg)
        return True, all_sample_ids


def run_snakefiles(input_dict, input_dir, output_dir,  config):
    input_dir = config["input_data_dir"]
    # SNAKEMAKE_PATH = config["snakemake"]
    # SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"
    SNAKEMAKE_PATH = "/opt/patric-common/runtime/artic-ncov2019/bin/snakemake"
    SNAKEFILE_DIR = "/home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_SARS2WasteWater/workflows/snakefile/"
    common_params = [
        SNAKEMAKE_PATH,
        "--cores", str(config["cores"]),
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going"
        ]

    if config["cores"] == 1:
        common_params.append("--debug")
        msg = "starting paired end reads\n"
        sys.stderr.write(msg)
    # process any paired reads
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_freyja_snakefile")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
        msg = "starting single end reads\n"
        sys.stderr.write(msg)
    # process any single reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_freyja_snakefile")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
    # once the snake files complete, check for the iVar bam files
    iVar_check = preprocessing_check(output_dir, input_dict)
    ## check that the iVar files exist ###
    file_check = iVar_check[0]
    all_sample_ids = iVar_check[1]

    if  file_check == True:
        # if all of the iVar bam files exist then the final step will trigger.
        # paired and single end reads will be processed together by Freyja
        print('starting freyja command')
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "freyja_snakefile")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
    # final file check
    freyja_check = post_processing_check(all_sample_ids, output_dir)
    if freyja_check == True:
        msg = "Wastewater Analysis is complete \n"
        sys.stderr.write(msg)
    return


def set_up_sample_dictionary(input_dir, input_dict, output_dir, cores):
    # set up the sample dictionary
    input_info = {}
    #### paired reads ####
    to_copy = []
    if len(input_dict["paired_end_libs"]) != 0:
        paired_sample_dict = {}
        ws_paired_reads = []
        ws_paired_reads = input_dict["paired_end_libs"]

        pe_path= f"{input_dir}/pe_reads"
        os.makedirs(pe_path, exist_ok = True)

        for i in range(len(ws_paired_reads)):
            # get filepaths and file names for both reads
            read1_filename = ws_paired_reads[i]["read1"].split("/")[-1]
            read1_filepath = check_input_fastqs(input_dir, read1_filename)

            # read 2
            read2_filename = ws_paired_reads[i]["read2"].split("/")[-1]
            read2_filepath = check_input_fastqs(input_dir, read2_filename)
            # pull out sample id
            sample_id = ws_paired_reads[i]["sample_id"].split("/")[-1]
            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
            # collect pairs of clean sample ids and file names for sample-meta-data-file
            # input_info[sample_id] = read1_filename
            input_info[read1_filename] = sample_id
            # set paths for zipped and unzipped files 
            if read1_filename.endswith(".gz"):
                pe_r1_samplename = f"{sample_id}_R1.fastq.gz"
                pe_r2_samplename = f"{sample_id}_R2.fastq.gz"
            else:
                pe_r1_samplename = f"{sample_id}_R1.fastq"
                pe_r2_samplename = f"{sample_id}_R2.fastq"
            paired_sample_dict[read1_filename] = pe_r1_samplename
            paired_sample_dict[read2_filename] = pe_r2_samplename

            to_copy.append([read1_filepath, f"{input_dir}/pe_reads/{pe_r1_samplename}"])
            to_copy.append([read2_filepath, f"{input_dir}/pe_reads/{pe_r2_samplename}"])

    #### single reads ####
    if len(input_dict["single_end_libs"]) != 0:
        single_end_sample_dict = {}
        ws_single_end_reads = []
        ws_single_end_reads = input_dict["single_end_libs"]
        ### relative path running from service-script ###
        se_path= f"{input_dir}/se_reads"
        os.makedirs(se_path, exist_ok = True)

        for i in range(len(ws_single_end_reads)):
            se_filename = ws_single_end_reads[i]["read"].split("/")[-1]
            se_filepath = check_input_fastqs(input_dir, se_filename)
            sample_id = ws_single_end_reads[i]["sample_id"]
            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
            # collect pairs of clean sample ids and file names for sample-meta-data-file
            input_info[se_filename] = sample_id
            if se_filename.endswith(".gz"):
                se_samplename = f"{sample_id}.fastq.gz"
            else:
                se_samplename = f"{sample_id}.fastq"
            single_end_sample_dict[se_filename] = se_samplename
            to_copy.append([se_filepath, f"{input_dir}/se_reads/{se_samplename}"])

    for (src, dest) in to_copy:
        shutil.copy(src, dest)

    # export the sample dictionary to .CSV
    with open(f"{output_dir}/sample_key.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["User sample name", "Analysis sample name"])
        if len(input_dict["paired_end_libs"]) != 0:
            for key, value in paired_sample_dict.items():
                writer.writerow([key, value])
        if len(input_dict["single_end_libs"]) != 0:
            for key, value in single_end_sample_dict.items():
                writer.writerow([key, value])
    return input_info


def write_to_warning_file(message):
    with open("warning.txt", "a") as warning_file:
        warning_file.write(message + "\n")

def write_to_job_failed_file(message):
    with open("JobFailed.txt", "a") as job_failed_file:
        job_failed_file.write(message + "\n")


# run the script from service-script/app_App-SARS2Wastewater.pl
# It takes a single argument which is the pathname of a config.json file
# This contains the app parameters in the params slot.
def main(argv):
    config_file = argv[0]
    print("Wrapper command recieved \n")
    try:
        fh = open(config_file)
        config = json.load(fh)
        fh.close()
    except IOError:
        msg = f"Could not read params file {config_file}"
        sys.stderr.write(msg)
        sys.exit(1)
    except json.JasonDecodeError as e:
        msg = f"JSON parse error in pyfile {config_file}: {e}"
        sys.stderr.write(msg)
        sys.exit(1)

    input_dict = config["params"]
    input_dir = config["input_data_dir"]
    output_dir = config["output_data_dir"]
    staging_metadata_file = config["staging_sample_metadata_path"]

    # rename files according to sample ids, copy files to pe_reads and se_reads in the staging directory
    input_info = set_up_sample_dictionary(input_dir, input_dict, output_dir, min(8, int(config["cores"])))
    # if a user provided metadata csv exists, make sure it matches the reads given
    check_sample_metadata_csv(input_info, config_file, input_dir, staging_metadata_file)
    # run the snakefiles
    run_snakefiles(input_dict, input_dir, output_dir, config)
if __name__ == "__main__":
    main(sys.argv[1:])