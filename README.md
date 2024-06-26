# SARS-Cov-2 Wastewater Analysis Service¶


## Overview
The SARS-CoV-2 Wastewater Analysis service is a comprehensive analysis of wastewater aimed at detecting and quantifying lineages and variants of concern (VOC) of the SARS-CoV-2 virus.

The service analyzes raw short amplicon reads by aligning them to the reference genome (Wuhan-Hu-1) and then performs variant analysis using Freyja. Freyja is a tool to identify and recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses lineage-determining mutational “barcodes” with information from the UShER global phylogenetic tree. We manage updating the barcodes to provide you up to date variant and lineage assignments. The results of this analysis workflow include sample processing status, key variant calling and alignment statistics, and sequencing depth coverage plots. It also provides lineage and VOC abundance plots by sample, date, week, and month for tracking the prevalence and distribution of different variants over time to aid public health response.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There is one application service specification defined here:

1.  [SARS-Cov-2 Wastewater Analysis Service](https://github.com/nicolegobo/bvbrc_SARS2Wastewater/blob/master/app_specs/SARS2Wastewater.json): Service that that provides the backend for the BV-BRC web inerface; it takes reads as input.

The code in this module provides the BV-BRC application service wrapper scripts for the wastewater service:

| Script name | Purpose |
| ----------- | ------- |
| [App-SARSS2Wastewater.pl](service-scripts/App-SARS2Wastewater.pl) | App script for the [SARS-Cov-2 Wastewater Analysis Service](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_wastewater_analysis_service.html) |

## See also

* [SARS-Cov-2 Wastewater Analysis Service](https://www.bv-brc.org/app/SARS2Wastewater)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_wastewater_analysis_service.html)
* [SARS-Cov-2 Wastewater Analysis Tutorial](https://www.bv-brc.org/docs/tutorial/sars_cov_2_wastewater/sars_cov_2_wastewater.html)

## References

Wattam AR, Davis JJ, Assaf R, Boisvert S, Brettin T, Bun C, Conrad N, Dietrich EM, Disz T, Gabbard JL, et al. 2017. Improvements to PATRIC, the all-bacterial Bioinformatics Database and Analysis Resource Center. Nucleic Acids Res 45:D535-D542.

Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018. 34(18): p. 3094-3100.

Etherington, G.J., R.H. Ramirez-Gonzalez, and D. MacLean, bio-samtools 2: a package for analysis and visualization of sequence and alignment data with SAMtools in Ruby. Bioinformatics, 2015. 31(15): p. 2565-2567.

Grubaugh, N.D., Gangavarapu, K., Quick, J. et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019). https://doi.org/10.1186/s13059-018-1618-7

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Karthikeyan, S., Levy, J.I., De Hoff, P. et al. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature 609, 101–108 (2022). https://doi.org/10.1038/s41586-022-05049-6

