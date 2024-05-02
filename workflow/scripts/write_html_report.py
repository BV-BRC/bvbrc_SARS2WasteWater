import base64
import json
import os
import pandas as pd
import re

# find assembly plots ending in detail.png
def find_assembly_plots(directory):
    assembly_plots = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.lower().endswith('detail.png'):  # Ensures case insensitivity
                full_path = os.path.join(dirpath, filename)
                assembly_plots.append(full_path)
    return assembly_plots

def generate_progress_cell(value):
    # Check if the value is "Complete" or "Incomplete"
    if value == "Complete":
        color = "lightgreen; opacity: 0.7"
    elif value == "Incomplete":
        color = "lightcoral; opacity: 0.7"
    else:
        color = "transparent"  # No color for other values
    return f'<td style="background-color: {color}; color: black;">{value}</td>'

def generate_progress_table_html(df):
    # Generate table headers
    headers = ''.join(f'<th>{header}</th>' for header in df.columns)
    
    # Generate table rows
    rows = ''.join(
        f'<tr>{"".join(generate_progress_cell(row[column]) for column in df.columns)}</tr>'
        for _, row in df.iterrows()
    )
    return headers, rows

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

def read_plotly_html(plot_path):
    # Read the content from 'Variant_Plot_Interactive.html'
    with open(plot_path, 'r') as file:
        plotly_html_content = file.read()
    # Extract everything within the <body> tags
    extracted_content = re.findall(r'<body>(.*?)</body>', plotly_html_content, re.DOTALL)

    # Assuming extracted_content contains our needed Plotly graph initialization scripts
    plotly_graph_content = extracted_content[0] if extracted_content else ''
    return plotly_graph_content

def write_html_report(workflow_dir):
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
    # set up the logo
    bvbrc_logo_path = os.path.join(workflow_dir, "bv-brc-header-logo-bg.png")
    base64_string = image_to_base64(bvbrc_logo_path)
    bvbrc_logo_base64 = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
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
    #set up the barcodes
    barcode_version=read_barcode_version("barcode_version.txt")

    # set up the freyja plots 
    lineage_plot = "plots/lineages_plot.html"
    variant_plot = "plots/variants_plot.html"

    standard_plot_header, standard_plot_description, standard_lineage_plot, standard_variant_plot = (
        ("<h3>Lineage and Variant Abundance by Sample</h3>",
        "<p>Below are stacked bar graphs showing the relative abundance of variants (left) and lineages (right) across the samples.<p>",
        read_plotly_html(lineage_plot),
        read_plotly_html(variant_plot))
        if os.path.exists(lineage_plot) else
        ("", "", "", "")
    )
    # Time series plots - Day
    lineage_day = "plots/lineages_by_day_plot.html"
    variant_day = "plots/variants_by_day_plot.html"
    time_series_header, day_plot_header, day_plot_description, day_lineage_plot, day_variant_plot = (
        ("<h3>Time Series Plots</h3>",
        "<h3>Lineage and Variant Abundance by Date</h3>",
        "<p>Below are the smoothed stacked bar graphs showing the relative abundance \
            of variants (left) and lineages (right) by date, generated by aggregating the results from one or more \
            samples collected on the same date. This plot may help reveal short-term variations and spikes in data \
            that might be linked to specific events or daily human activities.</p>",
        read_plotly_html(lineage_day),
        read_plotly_html(variant_day))
        if os.path.exists(lineage_day) else
        ("", "", "", "", "")
    )
    # Time series plots - Week
    lineage_week = "plots/lineages_by_week_plot.html"
    variant_week = "plots/variants_by_week_plot.html"
    week_plot_header, week_plot_description, week_lineage_plot, week_variant_plot = (
        ("<h3>Lineage and Variant Abundance by Week</h3>",
        "<p>The stacked bar graphs below show the relative abundance of variants (left) and \
            lineages (right) by epiweek, generated by aggregating the results from one or more samples collected in the \
            same epiweek. An 'epiweek', short for epidemiological week, is a standard method of grouping days into weeks \
            for the purposes of public health and epidemiological tracking. An epiweek serves to create a consistent and \
            comparable method of collecting and analyzing data across different time periods and regions. An epiweek \
            typically begins on a Sunday and ends on a Saturday, consisting of seven days in total.<p>",
        read_plotly_html(lineage_week),
        read_plotly_html(variant_week))
        if os.path.exists(lineage_week) else
        ("", "", "", "")
    )
    # Time series plots - month
    lineage_month = "plots/lineages_by_month_plot.html"
    variant_month = "plots/variants_by_month_plot.html"
    month_plot_header, month_plot_description, month_lineage_plot, month_variant_plot = (
        ("<h3>Lineage and Variant Abundance by Month</h3>",
        "<p>The stacked bar graphs below show the relative abundance of variants (left) and \
            lineages (right), aggregated by collection month.</p>",
        read_plotly_html(lineage_month),
        read_plotly_html(variant_month))
        if os.path.exists(lineage_month) else
        ("", "", "", "")
    )
    ### write the report HTML ###
    html_template = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <style>
                body {{ font-family: Roboto, sans-serif; color: black; }}
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
                        color: black;
                    }}
                    .warning {{ color: black; }}
                    table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
                    th, td {{ padding: 5px; text-align: left; }}
                    img {{ width: 100%; max-width: 600px; height: auto; }}
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
                    .img {{
                        width: 100%; /* Make the image expand to fill the container */
                        max-width: 600px; /* Maximum width of the image */
            </style>
        </head>
        <body>
            <header>
                <div class="title">SARS-CoV-2 Wastewater Report</div>
                <a href="https://www.bv-brc.org/">
                    {bvbrc_logo_base64}
                    </a>
            </header>
            <p>This report presents the findings from a comprehensive analysis of wastewater samples using the \
                SARS-CoV-2 Wastewater Analysis Service at the BV-BRC [1], aimed at detecting and quantifying lineages \
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
                The raw reads are aligned to the reference genome (Wuhan-Hu-1, NC_045512) with Minimap2 with MiniMap2 [2]. Then \
                SAMtools [3] converts the aligned FASTQs into BAM files. SAMtools also sorts the aligned BAM files by the leftmost \
                coordinates. Then the primers and low-quality sequences are trimmed by iVAR [4].  FastQC [5] offers a range of \
                quality assessments for the raw FASTQ files, as well as the aligned and sorted BAM files.  Note: that \
                wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants<p>

                <p>The relative lineage abundance analysis is performed using Freyja [6], which uses lineage-determining mutational \
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
                        coverage (<em>depth</em>) of reads across the virus genome, and the middle value of coverage \
                        (<em>median</em>) and how much this coverage varies from the average (<em>standard deviation</em>).</li>
                    <br>
                    <li><strong>Depth Minimum and Maximum</strong>:These indicate the lowest and highest coverage observed \
                        in the sequencing of the viral genome.</li>
                    <br>
                    <li><strong>Total N Count and N Blocks</strong>: 'N' refers to ambiguous bases in the DNA sequence that \
                        could not be called. This metric shows the total count ambiguous bases and the number of blocks that \
                        they are grouped into.</li>
                <br>
                    <li><strong>Fasta Length</strong>:The length of the consensus sequence in a FASTA format file. Note: that \
                        wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants.</li>
                    <br>
                    <li><strong>Primers and Primer Count</strong>: Primers are short sequences used to initiate DNA amplification. \
                        The primer count indicates how many primers were used in the amplification reaction.</li>
                    <br>
                    <li><strong>Mapped and Unmapped Reads</strong>: Mapped reads are sequences that could be aligned to the \
                        reference genome (Wuhan-Hu-1), whereas unmapped reads contain sequences that could not be aligned, \
                        indicating potential contamination.</li>
                    <br>
                    <li><strong>Primer Trim Count and Percentage of Primers Trimmed</strong>:  This shows how many primers \
                        were removed during data cleanup to reduce errors and the percentage of total primers this count represents.</li>
                    <br>
                    <li><strong>Variant Count</strong>: The total number of different viral SARS-CoV-2 variants identified in the sample.</li>
                </ul>
                <br>
                <h3>Coverage Plots</h3>
                <p>Below are the coverage plots for each sample. The X-axis represents the positions along the reference genome and the Y-axis \
                shows the sequencing depth or coverage at each position. The plots show areas of high and low coverage relative to the \
                Wuhan-hu-1 reference genome.  Note: that wastewater consensus sequences generated from this workflow are likely to contain a mixture of variants.<p>
                <div class="image-row">
    """.format(
        bvbrc_logo_base64=bvbrc_logo_base64,
        assembly_table_headers=assembly_table_headers,
        assembly_table_rows=assembly_table_rows,
        stats_table_headers=stats_table_headers,
        stats_table_rows=stats_table_rows,
        )
    # Loop over each image path, convert to Base64, and add an <img> tag to the HTML
    for image_path in image_paths:
        base64_string = image_to_base64(image_path)
        img_tag = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
        html_template += img_tag
    html_template += """
            </div>
        </body>
        </html>"""
    second_half = """
    <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <style>
                .plot-container {{
                    display: flex;
                    justify-content: space-around;
                    }}
                .plot {{
                        width: 50%; /* Each plot takes up half of the container width */
                    }}
                    .caption {{
                        text-align: center;
                        font-size: 14px;
                        margin-top: 20px;
                        color: black;
                        font-family: Roboto;
                    }}
                body {{ font-family: Roboto, sans-serif; color: black; }}
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
                        color: black;
                    }}
                    .warning {{ color: black; }}
                    table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
                    th, td {{ padding: 5px; text-align: left; }}
                    img {{ width: 100%; max-width: 600px; height: auto; }}
            </style>
        </head>
        <body>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            {standard_plot_header}
            {standard_plot_description}
            <div class="plot-container">
                <div class="plot" id="plot1">
                    {standard_variant_plot} <!-- Direct embedding of the plot content -->
                </div>
                <div class="plot" id="plot2">
                    {standard_lineage_plot} <!-- Direct embedding of the plot content -->
                </div>
            </div>
            {day_plot_header}
            {day_plot_description}
            <div class="plot-container">
                <div class="plot" id="plot1">
                    {day_variant_plot} <!-- Direct embedding of the plot content -->
                </div>
                <div class="plot" id="plot2">
                    {day_lineage_plot} <!-- Direct embedding of the plot content -->
                </div>
            </div>
            {week_plot_header}
            {week_plot_description}
            <div class="plot-container">
                <div class="plot" id="plot1">
                    {week_variant_plot} <!-- Direct embedding of the plot content -->
                </div>
                <div class="plot" id="plot2">
                    {week_lineage_plot} <!-- Direct embedding of the plot content -->
                </div>
            </div>
            {month_plot_header}
            {month_plot_description}
            <div class="plot-container">
                <div class="plot" id="plot1">
                    {month_variant_plot} <!-- Direct embedding of the plot content -->
                </div>
                <div class="plot" id="plot2">
                    {month_lineage_plot} <!-- Direct embedding of the plot content -->
                </div>
            </div>
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
            <li>1.	Wattam AR, Davis JJ, Assaf R, Boisvert S, Brettin T, Bun C, Conrad N, Dietrich EM, Disz T, Gabbard JL, et al. 2017. Improvements to PATRIC, the all-bacterial Bioinformatics Database and Analysis Resource Center. Nucleic Acids Res 45:D535-D542.</a></li>
            <ol>
            <li>Etherington, G.J., R.H. Ramirez-Gonzalez, and D. MacLean, bio-samtools 2: a package for analysis and visualization of sequence and alignment data with SAMtools in Ruby. Bioinformatics, 2015. 31(15): p. 2565-2567.</a></li>
        
            <li>Grubaugh, N.D., Gangavarapu, K., Quick, J. et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019). https://doi.org/10.1186/s13059-018-1618-7</a></li>
        
            <li>Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018. 34(18): p. 3094-3100.</a></li>
        
            <li>Karthikeyan, S., Levy, J.I., De Hoff, P. et al. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature 609, 101–108 (2022). https://doi.org/10.1038/s41586-022-05049-6</a></li>
            <ol>
        </body>
        </html>
    """.format(
        barcode_version=barcode_version,
        standard_plot_header=standard_plot_header,
        standard_plot_description=standard_plot_description,
        standard_variant_plot=standard_variant_plot,
        standard_lineage_plot = standard_lineage_plot,
        time_series_header=time_series_header,
        day_plot_header=day_plot_header, 
        day_plot_description=day_plot_description,
        day_variant_plot=day_variant_plot,
        day_lineage_plot=day_lineage_plot,
        week_plot_header=week_plot_header,
        week_plot_description=week_plot_description,
        week_lineage_plot=week_lineage_plot,
        week_variant_plot=week_variant_plot,
        month_plot_header=month_plot_header,
        month_plot_description=month_plot_description,
        month_lineage_plot=month_lineage_plot,
        month_variant_plot=month_variant_plot,
        warning_header=warning_header,
        warning_text=warning_text,
        )
    return html_template, second_half
# step 0 get paths set up
with open("config.json", 'r') as file:
        config = json.load(file)
output_dir = config["output_data_dir"]
workflow_dir = config["workflow_dir"]

html_template, second_half = write_html_report(workflow_dir)

html_report_path = os.path.join(output_dir, "SARS2Wastewater_report.html")
with open(html_report_path, 'w') as file:
    file.write(html_template)
    file.write(second_half)
print(f"Generated HTML report at {html_report_path}.")


