import base64
import json
import os
import pandas as pd


def generate_html(output_html_path, workflow_dir):
    barcode_version = ""
    warning_text = ""
    bvbrc_logo_path = os.path.join(workflow_dir, "bv-brc-header-logo-bg.png")
    base64_string = image_to_base64(bvbrc_logo_path)
    bvbrc_logo_base64 = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
    ## read in the job stats
    job_stats_df = pd.read_csv("output/job_stats.tsv", sep="\t", header=0)
    # Define the columns for analysis profress dataframe
    progress_cols = ['Sample ID', 'Assembly', 'Freyja - Analysis', 'Freyja - Visualization']
    # Define the columns for the asssembly stats
    sats_cols = ['Sample ID', 'Depth Mean', 'Depth Median', 'Depth Standard Deviation',
                'Depth Minimum', 'Depth Maximum', 'Total N Count', 'N Blocks',
                'Fasta Length', 'Primers', 'Primer Count', 'Mapped Reads',
                'Unmapped Reads', 'Primer Trim Count', 'Percentage of Primers Trimmed',
                'Variant Count']
    # Create the progress DataFrame
    assembly_progress_df = job_stats_df[progress_cols]
    # Create the assembly DataFrame
    stats_df = job_stats_df[sats_cols]
    # Extract table data from dataframe
    assembly_table_headers, assembly_table_rows = generate_progress_table_html(assembly_progress_df)
    stats_table_headers, stats_table_rows = generate_table_html(stats_df)
    # Find standard plots 
    lineage_plot = "output/plots/lineage_plot.svg"
    variant_plot = "output/plots/variants_plot.svg"
    if os.path.exists(lineage_plot):
        standard_plot_svg_header = "<h3>Lineage and Variant Abundance by Sample</h4>"
        standard_plot_description = "<p>Below are stacked bar graphs showing the relative abundance of \
        variants (left) and lineages (right) across the samples.<p>"
        variant_svg = read_svg(variant_plot)
        lineage_svg = read_svg(lineage_plot)
    else:
        standard_plot_svg_header = ""
        standard_plot_description = ""
        variant_svg = ""
        lineage_svg = ""
    # Time series plots - Day
    lineage_day = "output/plots/lineages_by_day_plot.svg"
    variant_day = "output/plots/variants_by_day_plot.svg"
    # Time series plots - Week
    lineage_week = "output/plots/lineages_by_week_plot.svg"
    variant_week = "output/plots/variants_by_week_plot.svg"
    # Time series plots - Month
    lineage_month = "output/plots/lineages_by_month_plot.svg"
    variant_month = "output/plots/variants_by_month_plot.svg"
    if os.path.exists(lineage_day):
        time_series_header = "<h3>Time Series Plots</h3>"
        day_plot_svg_header = "<h3>Lineage and Variant Abundance by Date</h4>"
        day_plot_description = "<p>Below are the smoothed stacked bar graphs showing the relative abundance \
        of variants (left) and lineages (right) by date, generated by aggregating the results from one or more \
        samples collected on the same date. This plot may help reveal short-term variations and spikes in data \
        that might be linked to specific events or daily human activities.</p>"
        day_variant_svg = read_svg(variant_day)
        day_lineage_svg = read_svg(lineage_day)
    else:
        time_series_header = ""
        day_plot_svg_header = ""
        day_plot_description = ""
        day_variant_svg = ""
        day_lineage_svg = ""
    if os.path.exists(lineage_week):
        week_plot_svg_header = "<h3>Lineage and Variant Abundance by Week</h4>"
        week_plot_description = "<p>The stacked bar graphs below show the relative abundance of variants (left) and \
        lineages (right) by epiweek, generated by aggregating the results from one or more samples collected in the \
        same epiweek. An 'epiweek', short for epidemiological week, is a standard method of grouping days into weeks \
        for the purposes of public health and epidemiological tracking. An epiweek serves to create a consistent and \
        comparable method of collecting and analyzing data across different time periods and regions. An epiweek \
        typically begins on a Sunday and ends on a Saturday, consisting of seven days in total.<p>"
        week_variant_svg = read_svg(variant_week)
        week_lineage_svg = read_svg(lineage_week)
    else:
        week_plot_svg_header = ""
        week_plot_description = ""
        week_variant_svg = ""
        week_lineage_svg = ""
    if os.path.exists(lineage_month):
        month_plot_svg_header = "<h3>Lineage and Variant Abundance by Month</h4>"
        month_plot_description = "<p>The stacked bar graphs below show the relative abundance of variants (left) and \
        lineages (right), aggregated by collection month.<p>"
        month_variant_svg = read_svg(variant_month)
        month_lineage_svg = read_svg(lineage_month)
    else:
        month_plot_svg_header = ""
        month_plot_description = ""
        month_variant_svg = ""
        month_lineage_svg = ""
    #month_header = "<h3>Lineage and Variants Summarized by Month</h4>"

    ### write the report ###
    if os.path.exists("output/warning.txt"):
        warning_header="<h3>Analysis Warnings</h4>"
        # Read the warning text
        with open("output/warning.txt") as f:
            warning_text = f.read()
    else:
        warning_header="<h3>Analysis Warnings</h4>"
        warning_text="No warnings were generated during your analysis."

    # assembly plots
    image_paths = find_assembly_plots("assembly_images")
    # barcode version
    barcode_version=read_barcode_version("barcode_version.txt")
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
<<<<<<< HEAD
        <style>
            body {{ font-family: Roboto, sans-serif; }}
            header {{
                display: flex;
                justify-content: space-between;
                align-items: center;
                padding: 10px 20px;
            }}
            header > a img {{
                max-width: 225px;  /* Maximum width */
                max-height: 225px;  /* Maximum height */
                width: auto;
                height: auto;
            }}
            .title {{
                font-size: 36px;  /* Adjust the size of the title text */
                font-family: 'Roboto', sans-serif;
                font-weight: bold;
            }}
=======
        <title>SARS-CoV-2 Wastewater Report</title>
        <style>
            body {{ font-family: Roboto, sans-serif; }}
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
            .warning {{ color: black; }}
            table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
            th, td {{ padding: 5px; text-align: left; }}
            img {{ width: 100%; max-width: 600px; height: auto; }}
            .svg-container {{
                    display: flex;
                    justify-content: space-around;
                    align-items: center;
                    padding: 10px 0;
                }}
                .svg-plot {{
                    flex: 1 1 auto;
                    margin: 10px;
                    overflow: hidden;
                }}
                .svg-plot svg {{
                    display: block;
                    max-width: 100%;
                    height: auto;
                    transform: scale(0.99);
                }}
            .image-row {{
                display: flex;
                flex-wrap: wrap;
                justify-content: flex-start;
            }}
            .image-container {{
                width: 33%; /* Each image container takes up one-third of the row */
                padding: 5px; /* Padding around the images */
                box-sizing: border-box;
            }}
            img {{
                width: 100%; /* Make the image expand to fill the container */
                max-width: 600px; /* Maximum width of the image */
                height: auto; /* Maintain aspect ratio */
            }}
        </style>
    </head>
    <body>
<<<<<<< HEAD
        <header>
        <div class="title">SARS-CoV-2 Wastewater Report</div>
        <a href="https://www.bv-brc.org/">
            {bvbrc_logo_base64}
            </a>
        </header>
        <p>This report presents the findings from a comprehensive analysis of wastewater samples using the \
        SARS-CoV-2 Wastewater Analysis Service at the BV-BRC[1], aimed at detecting and quantifying lineages \
        and variants of concern (VOC) of the SARS-CoV-2 virus. The service analyzes raw reads by aligning them \
        to the reference genome (Wuhan-Hu-1) and then performs variant analysis using Freyja. This report \
        provides the description of the overall analysis workflow, sample processing status, key variant calling \
        and alignment statistics, and sequencing depth coverage plots. It also provides lineage and VOC abundance \
        plots by sample, date, week, and month for tracking the prevalence and distribution of different variants \
        over time to aid public health response.</p>
        <h3>Analysis Workflow</h4>
        <p>The service accepts one or more FASTQ files, either uploaded directly by the user or retrieved via SRA Run \
        Accessions.  These files should contain short reads generated by amplicon-based sequencing of wastewater \
        samples and require metadata about the primers used and sample identifiers. Collection dates are optional. \
        The raw reads are aligned to the reference genome (Wuhan-Hu-1, NC_045512) with Minimap2 with MiniMap2[2]. Then \
        SAMtools[3] converts the aligned FASTQs into BAM files. SAMtools also sorts the aligned BAM files by the leftmost \
        coordinates. Then the primers and low-quality sequences are trimmed by iVAR[4].  FastQC[5] offers a range of \
        quality assessments for the raw FASTQ files, as well as the aligned and sorted BAM files.  Note: that \
        wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants<p>
        <br>
        <p>The relative lineage abundance analysis is performed using Freyja[6], which uses lineage-determining mutational \
        “barcodes” derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, \
        non-negative) de-mixing problem. Using SNV frequency and sequencing depth at each position in the genome, Freyja \
        returns an estimate of the true lineage abundances in each sample. The results are summarized by sample or \
        aggregated by date, week, and month and presented as lineage and variant abundance plots. <p>
        <!-- Analysis progress results CSV Table -->
        <h3>Sample Processing</h4>
        <p>The table below tracks the progress of each sample through three main stages of analysis workflow: aligning the \
        reads and variant calling, Freyja - Analysis, and Freyja - Visualization. If a sample is labeled as incomplete for \
        any of the stages, please refer to the assembly and alignment statistics table in the following section or the \
        warnings at the end of the report.</p>
=======
        <h1>SARS-Cov-2 Wastewater Analysis Report</h1>
        <p>This report presents the findings from a comprehensive analysis of wastewater samples using the \
        SARS-CoV-2 Wastewater Analysis Service at the BV-BRC, aimed at detecting and quantifying lineages and \
        variants of concern (VOC) of the SARS-CoV-2 virus. The service assembles raw reads by aligning them to \
        the reference genome (Wuhan-Hu-1) using One Codex pipeline and then performs variant analysis using Freyja.\
        The report provides the description of the overall analysis workflow, sample processing status, key assembly \
        and alignment statistics, and sequencing depth and coverage plots. It also provides lineage and VOC abundance \
        plots by sample, date, week, and month for tracking the prevalence and distribution of different variants over \
        time to aid public health response.</p>
        <h3>Analysis Workflow</h4>
        <p>The service takes as input one or more FASTQ files containing short reads generated by amplicon-based \
        sequencing of wastewater samples or via SRA Run Accessions, along with metadata about the primers used, sample \
        identifiers collection dates.  The raw reads are first analyzed for quality using FastQC, trimmed for primers \
        and low quality sequences using iVar, and then aligned to the reference genome (Wuhan-Hu-1) using the OneCodex \
        pipeline for generating consensus sequence and variant calling. The relative linage and variant abundance analysis \
        is performed using Freyja, which uses lineage-determining mutational “barcodes” derived from the UShER global \
        phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. Using SNV \
        frequency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage \
        abundances in each sample. The results are summarized by sample or aggregated by date, week, and month and presented \
        as lineage and variant abundance plots.<p>
        <!-- Analysis progress results CSV Table -->
        <h3>Sample Processing</h4>
        <p>The table below tracks the progress of each sample through three main stages of analysis workflow: Assembly, \
        Freyja - Analysis, and Freyja - Visualization. If a sample is labeled as incomplete for any of the stages, please \
        refer to the assembly and alignment statistics table in the following section or the warnings at the end of the \
        report.</p>
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
        <table>
            <thead>
                <tr>
                    {assembly_table_headers}
                </tr>
            </thead>
            <tbody>
                {assembly_table_rows}
            </tbody>
        </table>
        <br>
        <!-- Assembly results CSV Table -->
        <h3>Primer Trimming and Alignment Statistics</h4>
        <p>The following table provides and overview of the key statistics from the primer trimming, removal of low-quality \
        sequence, alignment to the Wuhan-Hu-1 reference genome, and variant calling. Please refer to the user guide for a \
        detailed description of the values in the table.</p>
        <table>
            <thead>
                <tr>
                    {stats_table_headers}
                </tr>
            </thead>
            <tbody>
                {stats_table_rows}
            </tbody>
        </table>
        <p>The columns are as follows<p>
        <ul>
            <li><strong>Depth Mean, Median, and Standard Deviation</strong>: These statistics describe the average \
<<<<<<< HEAD
                coverage (<em>depth</em>) of reads across the virus genome, and the middle value of coverage \
                (<em>median</em>) and how much this coverage varies from the average (<em>standard deviation</em>).</li>
            <br>
            <li><strong>Depth Minimum and Maximum</strong>:These indicate the lowest and highest coverage observed \
                in the sequencing of the viral genome.</li>
            <br>
            <li><strong>Total N Count and N Blocks</strong>: 'N' refers to ambiguous bases in the DNA sequence that \
                could not be called. This metric shows the total count ambiguous bases and the number of blocks that \
=======
                coverage (<em>depth</em>) of sequencing across the virus genome, and the middle value of coverage \
                (<em>median</em>) and how much this coverage varies from the average (<em>standard deviation</em>).</li>
            <br>
            <li><strong>Depth Minimum and Maximum</strong>:These indicate the lowest and highest coverage observed \
                in the sequencing of the virus genome.</li>
            <br>
            <li><strong>Total N Count and N Blocks</strong>: 'N' refers to bases in the DNA sequence that could not \
                be identified. This metric shows the total count of such unidentified bases and the number of blocks \
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
                they are grouped into.</li>
           <br>
            <li><strong>Fasta Length</strong>:The length of the consensus sequence in a FASTA format file. Note: that \
                wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants.</li>
            <br>
<<<<<<< HEAD
            <li><strong>Primers and Primer Count</strong>: Primers are short sequences used to initiate DNA amplification. \
                The primer count indicates how many primers were used in the amplification reaction.</li>
            <br>
            <li><strong>Mapped and Unmapped Reads</strong>: Mapped reads are sequences that could be aligned to the \
                reference genome (Wuhan-Hu-1), whereas unmapped reads contain sequences that could not be aligned, \
                indicating potential contamination.</li>
=======
            <li><strong>Primers and Primer Count</strong>: Primers are short sequences used to initiate DNA synthesis. \
                The primer count indicates how many times primers were used in the sequencing process.</li>
            <br>
            <li><strong>Mapped and Unmapped Reads</strong>: Mapped reads are sequences that could be aligned to the \
                reference genome (Wuhan-Hu-1), whereas unmapped reads could not be aligned, indicating potential \
                variations or errors.</li>
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
            <br>
            <li><strong>Primer Trim Count and Percentage of Primers Trimmed</strong>:  This shows how many primers \
                were removed during data cleanup to reduce errors and the percentage of total primers this count represents.</li>
            <br>
<<<<<<< HEAD
            <li><strong>Variant Count</strong>: The total number of different viral SARS-CoV-2 variants identified in the sample.</li>
        </ul>
        <br>
        <h3>Coverage Plots</h3>
        <p>Below are the coverage plots for each sample. The X-axis represents the positions along the reference genome and the Y-axis \
        shows the sequencing depth or coverage at each position. The plots show areas of high and low coverage relative to the \
        Wuhan-hu-1 reference genome.  Note: that wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants.<p>
=======
            <li><strong>Variant Count</strong>: The total number of different viral variants identified in the sample.</li>
        </ul>
        <br>
        <h3>Coverage Plots</h3>
        <p>Below are coverage plots are generated to show areas of high and low coverage relative to the Wuhan-hu-1 reference genome.\
        Note that wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants.<p>
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
        <div class="image-row">
            """
    # Loop over each image path, convert to Base64, and add an <img> tag to the HTML
    for image_path in image_paths:
        base64_string = image_to_base64(image_path)
        img_tag = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
        html_template += img_tag

    # Close the HTML tags
    html_template += """
    </div>
        <!-- SVG Images Placement -->
        {standard_plot_svg_header}
        {standard_plot_description}
        <!--  Standard plots -->
        <div class="svg-container">
            <div class="svg-plot">{variant_svg}</div>
            <div class="svg-plot">{lineage_svg}</div>
        </div>
        <br>
        {time_series_header}
        <!-- Day Plots -->
        {day_plot_svg_header}
        {day_plot_description}
            <div class="svg-container">
            <div class="svg-plot">{day_variant_svg}</div>
            <div class="svg-plot">{day_lineage_svg}</div>
        </div>
        <!-- week Plots -->
        {week_plot_svg_header}
        {week_plot_description}
            <div class="svg-container">
            <div class="svg-plot">{week_variant_svg}</div>
            <div class="svg-plot">{week_lineage_svg}</div>
        </div>
        <!-- month Plots -->
        {month_plot_svg_header}
        {month_plot_description}
            <div class="svg-container">
            <div class="svg-plot">{month_variant_svg}</div>
            <div class="svg-plot">{month_lineage_svg}</div>
        </div>
        <div>
        <h3>Barcode Version</h3>
        <p>Freyja uses lineage-determining mutational "barcodes" derived from the UShER global phylogenetic \
        tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. The barcodes \
        are generated by the authors of Freyja. The BV-BRC will keep these barcodes up to date. The barcode \
        version can impact your results. These samples were analyzed with barcode version {barcode_version}.<p>
        </div>
        <!-- Warning Text -->
        <div class="warning">
            {warning_header}
            {warning_text}
        </div>
    <h3>References</h3>
    <ol>
<<<<<<< HEAD
    <li>1.	Wattam AR, Davis JJ, Assaf R, Boisvert S, Brettin T, Bun C, Conrad N, Dietrich EM, Disz T, Gabbard JL, et al. 2017. Improvements to PATRIC, the all-bacterial Bioinformatics Database and Analysis Resource Center. Nucleic Acids Res 45:D535-D542.</a></li>
    <ol>
=======
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
    <li>Etherington, G.J., R.H. Ramirez-Gonzalez, and D. MacLean, bio-samtools 2: a package for analysis and visualization of sequence and alignment data with SAMtools in Ruby. Bioinformatics, 2015. 31(15): p. 2565-2567.</a></li>

    <li>Grubaugh, N.D., Gangavarapu, K., Quick, J. et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019). https://doi.org/10.1186/s13059-018-1618-7</a></li>

    <li>Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018. 34(18): p. 3094-3100.</a></li>

    <li>Karthikeyan, S., Levy, J.I., De Hoff, P. et al. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature 609, 101–108 (2022). https://doi.org/10.1038/s41586-022-05049-6</a></li>
    <ol>
    </body>
    </html>
    """

    # Format the template with dynamic data
    formatted_html = html_template.format(
<<<<<<< HEAD
        bvbrc_logo_base64=bvbrc_logo_base64,
=======
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883
        barcode_version=barcode_version,
        standard_plot_svg_header=standard_plot_svg_header,
        standard_plot_description=standard_plot_description,
        lineage_svg=lineage_svg,
        variant_svg=variant_svg,
        time_series_header=time_series_header,
        day_plot_svg_header=day_plot_svg_header,
        day_plot_description=day_plot_description,
        day_variant_svg=day_variant_svg,
        day_lineage_svg=day_lineage_svg,
        week_plot_svg_header=week_plot_svg_header,
        week_plot_description=week_plot_description,
        week_variant_svg=week_variant_svg,
        week_lineage_svg=week_lineage_svg,
        month_plot_svg_header=month_plot_svg_header,
        month_plot_description=month_plot_description,
        month_variant_svg=month_variant_svg,
        month_lineage_svg=month_lineage_svg,
        ##### ### ######
        assembly_table_headers=assembly_table_headers,
        assembly_table_rows=assembly_table_rows,
        stats_table_headers=stats_table_headers,
        stats_table_rows=stats_table_rows,
        warning_header=warning_header,
        warning_text=warning_text,
    )
    # Write to file
    with open(output_html_path, 'w') as file:
        file.write(formatted_html)
    print(f"Generated HTML report at {output_html_path} including SVG images.")


def generate_progress_table_html(df):
    # Generate table headers
    headers = ''.join(f'<th>{header}</th>' for header in df.columns)
    
    # Generate table rows
    rows = ''.join(
        f'<tr>{"".join(generate_progress_cell(row[column]) for column in df.columns)}</tr>'
        for _, row in df.iterrows()
    )
    return headers, rows

def generate_progress_cell(value):
    # Check if the value is "Complete" or "Incomplete"
    if value == "Complete":
        color = "lightgreen; opacity: 0.7"
    elif value == "Incomplete":
        color = "lightcoral; opacity: 0.7"
    else:
        color = "transparent"  # No color for other values

    # return f'<td style="background-color: {color};">{value}</td>'
    return f'<td style="background-color: {color}; color: black;">{value}</td>'

# Function to generate the HTML table rows
def generate_table_html(df):
    # Generate table headers
    headers = ''.join(f'<th>{header}</th>' for header in df.columns)

    # Generate table rows
    rows = ''.join(
        f'<tr>{" ".join(f"<td>{row[column]}</td>" for column in df.columns)}</tr>'
        for _, row in df.iterrows()
    )
    
    return headers, rows

### this is how you can find the png plots 
def find_assembly_plots(directory):
    assembly_plots = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.lower().endswith('detail.png'):  # Ensures case insensitivity
                full_path = os.path.join(dirpath, filename)
                assembly_plots.append(full_path)
    return assembly_plots

def image_to_base64(file_path):
    """Read binary file and return base64-encoded string"""
    with open(file_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string

def read_barcode_version(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        return ""

def read_svg(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        return ""

# step 0 get paths set up
with open("config.json", 'r') as file:
        config = json.load(file)
output_dir = config["output_data_dir"]
<<<<<<< HEAD
workflow_dir = config["workflow_dir"]
html_report_path = os.path.join(output_dir, "SARS2Wastewater_report.html")


generate_html(html_report_path, workflow_dir)
=======
html_report_path = os.path.join(output_dir, "SARS2Wastewater_report.html")


generate_html(html_report_path)
>>>>>>> 029e4959d9aab47399306fa6dabeb42c7fef9883



