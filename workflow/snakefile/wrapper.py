#!/opt/patric-common/runtime/bin/python3
import csv
import json 
import os
import re
import sys
import shutil
import subprocess

# def check_input_fastqs(input_dir, filename, files_to_zip):
def check_input_fastqs(input_dir, filename):
    input_path = f"{input_dir}/{filename}"
    if os.path.isfile(input_path) == True:
        return input_path
    else:
        msg = f"Error {input_path} does not exisit'. \n check on {input_path} \n"
        sys.stderr.write(msg)
        sys.exit(1)
        pass

def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict


def post_processing_check(all_sample_ids, output_dir):
    dict_samples = {}
    complete = []
    incomplete = []
    for sample_name in all_sample_ids:
        # check for missing Freyja results
        freyja_path = f"{output_dir}/{sample_name}/{sample_name}.freyja"
        if os.path.isfile(freyja_path) == True:
            msg = f"{sample_name}"
            complete.append(msg)
        else:
            dict_samples[sample_name] = False
            msg = f"{sample_name}"
            incomplete.append(msg)
    ## write the complete and incomplete samples to the error report for user to see ##
    if len(incomplete) != 0:
        msg = f"Freyja results produced for the following samples: {complete}. \n \
                Freyja results missing: Review the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
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

    # check if iVar bam exisits for each file 
    for sample_name in all_sample_ids:
        iVar_bam = f"{output_dir}/{sample_name}/assembly/{sample_name}.ivar.bam"


        if os.path.isfile(iVar_bam) == True:
            complete.append(sample_name)
        else:
            dict_samples[sample_name] = False
            incomplete.append(sample_name)
    if len(incomplete) != 0:
        msg = f"Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing \n \
                Primer trimming and variant calling is complete for the following samples: {complete}. \n \
                Primer trimming and variant calling is INCOMPLETE for the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
        sys.exit(1)
        ## return False
    else:
        msg = "preprocessing complete \n"
        sys.stderr.write(msg)
        return True, all_sample_ids


def run_snakefiles(input_dict, input_dir, output_dir,  config):

    input_dir = config["input_data_dir"]
    SNAKEMAKE_PATH = config["snakemake"]
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"

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
    
    # process any paired reads
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_freyja_snakefile")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)

    # process any single reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_freyja_snakefile")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
    # once the snake files complete, check for the iVar bam files
    iVar_check = preprocessing_check(output_dir, input_dict)
    ### check that the iVar files exist ###
    file_check = iVar_check[0]
    all_sample_ids = iVar_check[1]

    if  file_check == True:
        # if all of the iVar bam files exist then the final step will trigger.
        # paired and single end reads will be processed together by Freyja
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
    #### paired reads ####
    to_copy = []
    if len(input_dict["paired_end_libs"]) != 0:
        paired_sample_dict = {}
        ws_paired_reads = []
        ws_paired_reads = input_dict["paired_end_libs"]

        pe_path= f"{input_dir}/pe_reads"
        os.makedirs(pe_path, exist_ok = True)

        for i in range(len(ws_paired_reads)):
            read1_filename = ws_paired_reads[i]["read1"].split("/")[-1]
            read1_filepath = check_input_fastqs(input_dir, read1_filename)

            read2_filename = ws_paired_reads[i]["read2"].split("/")[-1]
            read2_filepath = check_input_fastqs(input_dir, read2_filename)
            sample_id = ws_paired_reads[i]["sample_id"].split("/")[-1]
            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
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
            #sample_ids.append(sample_id)

            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
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
    return


# run the script from service-script/app_taxonomic_classification perl script
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
        print(f"Could not read params file {config_file}")
        sys.exit(1)
    except json.JasonDecodeError as e:
        print(f"JSON parse error in pyfile {config_file}: {e}")
        sys.exit(1)

    input_dict = config["params"]
    input_dir = config["input_data_dir"]
    output_dir = config["output_data_dir"]

    set_up_sample_dictionary(input_dir, input_dict, output_dir, min(8, int(config["cores"])))
    # run the snakefiles
    run_snakefiles(input_dict, input_dir, output_dir, config)

if __name__ == "__main__":
    main(sys.argv[1:])