# FRBAL: Fast Radio Burst Analysis Library

This code searches for Fast Radio Bursts in data collected from the Green Bank 20m radio telescope. There are two main steps to the pipeline: the first is a conversion of the raw observation files to .fil filterbank format (which makes them a lot more compact). The second is the search for strong dispersed pulses using the Fast Discrete Dispersion Measure Transform (FDMT; Zackay & Ofek 2017), with some crude shell scripts to run it in parallel across the academic computing cluster at Green Bank.

## File Conversion Pipeline
This pipeline is run on the 20m-data machine.
NOTE: there have been problems with this pipeline causing load issues when a high-resolution observation is being collected. I added ionice, but it is still probably wise to avoid running it during observations.

Search for long observations:
`find_long_observations.sh`

If you want to peek while it's running:
`tail -f ~/long_observations_sorted.txt`

Final results, sorted by size
`less long_observations_sorted.txt`

Select some observations
`cp long_observations_sorted.txt long_observations_selected.txt`

Remove the size field
`vim long_observations_selected.txt`

Convert them to fil files
`cat long_observations_selected.txt | xargs -L 1 fits2fil.sh`

## FRB Search Pipeline

Either use scratch/fil/process_results.sh to process all .fil files found in /scratch/fil (not very smart, just checks to see if there is a corresponding directory in /scratch/results). Or run `dist_fdmt_search.sh filename.fil` on a particular filename.
