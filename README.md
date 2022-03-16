# Mutational Signatures Data

This repository contains Python scripts to ingest and process data for [mutational signature analysis](https://en.wikipedia.org/wiki/Mutational_signatures), from public datasets in the [International Consortium for Cancer and Genomics (ICGC) Data Portal](https://dcc.icgc.org/), with an emphasis on adherence to [FAIR (Findable, Accessible, Interoperable, and Reusable) principles](https://www.go-fair.org/fair-principles/).

## Consuming ICGC Data
The jupyter notebook `notebooks/Consuming ICGC Data.ipynb` shows how to consume simple somatic mutation data from the ICGC Data Portal via its [APIs](https://docs.icgc.org/portal/api-endpoints/).

## Processing ICGC simple somatic mutation data
The jupyter notebook `notebooks/Processing ICGC Data.ipynb` shows how to process the simple somatic mutation (SSM) data that we downloaded from the ICGC Data Portal. After processing the SSM dataset, we generate a mutational spectra matrix, which tallies the counts for each context (as determined by the base on the 5' and the 3' ends of the mutated allele) of single base substitution. To obtain the context of the substitution, we require the GRCh37 reference genome. We obtain this genome via the [UCGC Genome Browser API](https://hgdownload.soe.ucsc.edu/downloads.html). The mutational spectra matrix for the BRCA datasets is stored in `data/mutational_spectra_wgs/`. This dataset can be readily used for mutational signature analysis.

## Simulating MSK-IMPACT 410 gene panel data from whole genome MAF files
The jupyter notebook `notebooks/Simulating MSK-IMPACT 410.ipynb` shows how to simulate a mutational spectra from MSK-IMPACT 410 gene panel from the MAF files generated from the ICGC SSM whole genome analysis dataset.


## Consuming MSK-IMPACT data and preparing its mutational spectra matrix 
The jupyter notebook `notebooks/Consuming and Preparing MSK-IMPACT Data.ipynb` shows how to use the [cBioPortal web API](https://www.cbioportal.org/webAPI) to download MSK-IMPACT data from the 2017 paper: [Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients](https://pubmed.ncbi.nlm.nih.gov/28481359/). The downloaded data is then converted into a mutational spectra matrix and stored at `data/msk_impact_nature_2017/`.