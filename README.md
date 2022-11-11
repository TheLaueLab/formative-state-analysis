# formative-state-analysis

1. identify contact hubs from single cell Hi-C contact map (Hub-identification.r)

This is the R functions used to identify contact hubs from single haploid cells 

Require package: dplyr (https://cran.r-project.org/web/packages/dplyr/index.html)

Input file: .ncc files generated from Nuc_process (https://github.com/tjs23/nuc_processing)
            a table specifying the start and the end of each chromosome of the given organism, here for example, single_cell_structure_model_chr_start_and_end.csv (this is for mouse mm10)

Additional notes: 
This is used mainly for single haploid cells, which only have only allele for any given loci, but can be generalised to the diploid case. 
If wishing to use this for diploid cells, please make sure you disambiguate the diploid contact map first before putting it into this pipeline to identify hubs. 
It's recommended to only use unambiguate contacts. 
When running the code, please run the R function definitions from line 11 to line 380 first, and then go to line 380 to set input path and run the main functions 
from line 384 to line 403
    
