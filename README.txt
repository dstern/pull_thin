This is a little utility for post-processing .tsv files produced by the Multiplexed Shotgun Genotyping (msg) pipeline

DEPENDENCIES

Python 2.7

USAGE

python pull_thin_tsv.0.5.1.py <config_file: default = pt.cfg> 

If no config_file provided, program looks for pt.cfg in same folder

INFORMATION
For info on msg, see 
github.com/JaneliaSciComp/msg

This utility will perform the following tasks
1 - Pull out specific individuals from across multiple tsv files
2 - Reduce each file to thinned data (defined by parameter difffac), keeping the flanking values with values < difffac
3 - Can replace the NAs in the files with priors, which is useful for some downstream applications, such as R-QTL


