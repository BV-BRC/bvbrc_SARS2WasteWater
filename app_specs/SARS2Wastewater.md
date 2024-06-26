
# Application specification: SARS2Wastewater

This is the application specification for service with identifier SARS2Wastewater.

The backend script implementing the application is [App-SARS2Wastewater.pl](../service-scripts/App-SARS2Wastewater.pl).

The raw JSON file for this specification is [SARS2Wastewater.json](SARS2Wastewater.json).

This service performs the following task:   Assemble SARS2 reads into a consensus sequence

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_libs |  | group  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| primers | Primer set to use for assembly | enum  | :heavy_check_mark: | ARTIC |
| minimum_base_quality_score | The minimum base quality score | int  |  | 20 |
| minimum_genome_coverage | The minimum genome coverage | int  |  | 60 |
| agg_minimum_lineage_abundance | Minimum lineage abundance for the plot command | float  |  | 0.01 |
| minimum_coverage_depth | The minimum coverage depth minimum.  | int  |  | 0 |
| confirmedonly | Excludes unconfirmed lineages from the analysis. | bool  |  | 0 |
| minimum_lineage_abundance | Minimum lineage abundance | float  |  | 0.001 |
| coverage_estimate | coverage estimate value | int  |  | 10 |
| timeseries_plot_interval | Timeseries plot interval | string  |  | 0 |
| primer_version | Version number for primer | string  |  |  |
| barcode_csv | Custom barcode path | string  |  |  |
| sample_metadata_csv | Sample metadata csv | string  |  | 0 |
| keep_intermediates | Keep all intermediate output from the pipeline | bool  |  | 1 |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| debug_level | Debug level | int  |  | 0 |

