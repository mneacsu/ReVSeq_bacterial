# ReVSeq - Bacterial Reads Analysis

This repository contains code for the analysis of bacterial reads detected in ReVSeq samples.

## Extract bacterial info from a set of samples 

The `process_plate.py` script extracts the following information for each sample: the number & percentage of bacterial reads, the percentage of bacterial reads mapping to each of the 6 large phyla, potentially pathogenic species, and optionally, virulence factors present in the sample (using [ARIBA](https://github.com/sanger-pathogens/ariba/wiki)).

Example usage:

```
python3 -W ignore process_plate.py kraken2_reports_dir raw_reads_dir metadata_file output_dir --search_VF -t 16
```
## Visualize species composition

[KronaTools](https://github.com/marbl/Krona/wiki) can be used to visualize the species composition based on kraken2 reports or any tables containing taxons and their quantifications.

Parameters:  
`-m` is the index of the quantifications column (= 3 in kraken2 reports)
`-t` is the index of the taxon id column (= 5 in kraken2 reports)
`-c` allows combining multiple reports into a single Krona chart. Unless specified, one chart is created for each input report.

Example usage:

```
ktImportTaxonomy -m 3 -t 5 -c -o all_samples.html kraken2_results/*/*/*_genus_filtered.tsv
```

## Miscellaneous

`get_detected_genera.py` - list all genera found in kraken2 reports
`taxonomy.R` - filter list of taxonomic names (e.g. keep only bacterial genera)
`genera_statistics.py` - rank bacterial genera by read abundance across all samples
`statistics_by_sample.py` - get percentage of bacterial reads mapping to the top 20 most abundant genera for each sample
`format_metadata.py` - format metadata table (reduce to one row per sample, get PCR/sequencing type, etc.)