# SARS-Cov-2 Wastewater Analysis Service

## Overview

The SARS-CoV-2 Wastewater Analysis service is a comprehensive analysis of wastewater aimed at detecting and quantifying lineages and variants of concern (VOC) of the SARS-CoV-2 virus. 

The service accepts raw short amplicon reads from wastewater samples. The service analyzes reads by aligning them to the reference genome (Wuhan-Hu-1) and then analyzes the variants in the sample using [Freyja](https://andersen-lab.com/secrets/code/). Below is an overview of the service.

![SARSWastewater overview](../images/sars_wastewater_service/image_1_workflow_image.png "SARSWastewater overview")



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [SARS2Wastewater](app_specs/SARS2Wastewater.md)


## See also

* [SARS-Cov-2 Wastewater Analysis Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_wastewater_analysis_service.html)
  * [SARSWastewater Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/SARSWastewater.html)
  * [SARSWastewater Service Tutorial](https://www.bv-brc.org/docs/../../tutorial/sarswastewater/sarswastewater.html.html)

![submission page](../images/sars_wastewater_service/image_2_submission_page.png "submission page")



## References

1.	Wattam AR, Davis JJ, Assaf R, Boisvert S, Brettin T, Bun C, Conrad N, Dietrich EM, Disz T, Gabbard JL, et al. 2017. Improvements to PATRIC, the all-bacterial Bioinformatics Database and Analysis Resource Center. Nucleic Acids Res 45:D535-D542.
2.	Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018. 34(18): p. 3094-3100.
3.	Etherington, G.J., R.H. Ramirez-Gonzalez, and D. MacLean, bio-samtools 2: a package for analysis and visualization of sequence and alignment data with SAMtools in Ruby. Bioinformatics, 2015. 31(15): p. 2565-2567.
4.	Grubaugh, N.D., Gangavarapu, K., Quick, J. et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol 20, 8 (2019). https://doi.org/10.1186/s13059-018-1618-7
5.	Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
6.	Karthikeyan, S., Levy, J.I., De Hoff, P. et al. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature 609, 101–108 (2022). https://doi.org/10.1038/s41586-022-05049-6


