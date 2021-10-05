# Methods:
The methods for this project can be divided into two sections, creating a TE annotation, and running TE Density and analyzing the results.

## TE Annotation:
TE annotations were created with [EDTA](https://github.com/oushujun/EDTA)
Analysis scripts for running EDTA can be found at [Fragaria TE Annotations](https://github.com/sjteresi/Fragaria_TE_Annotations).
The `Makefile` contains a record of what commands and scripts were performed; scripts can be located within the `src/` directory.
EDTA v1.9.7 was used to annotate and characterize the TEs in the genome. 
EDTA was run independently for each genome.
A FASTA file and a CDS FASTA file were used as inputs to EDTA for each genome.
The CDS FASTA was generated using [GFFRead](https://github.com/gpertea/gffread) v0.12.6.
EDTA was then run with default options in all cases except for the usage of the `--cds`, `--sensitive 1`, and `--anno` options.

## TE Density Analysis:
TE density was assessed using [TE Density](https://github.com/sjteresi/TE_Density), (Teresi, manuscript in review).
Analysis scripts for running TE Density can be located at [CFRA Syntelog TE Differences](https://github.com/sjteresi/CFRA_Syntelog_TE_Differences). 
The `Makefile` contains a record of what commands and scripts were performed; scripts can be located within the `src/` directory.
Prior to running TE Density, a custom script was used to modify some of the TE groupings in order to reorganize, simplify, and conform to the naming system presented in Wicker et al 2007.
The naming scheme used to reorganize the TE groupings can be found at `src/replace_names_strawberry.py.`
A gene annotation along with the renamed TE annotation were used as primary inputs to TE Density v0.9.
TE Density was run with default options for all genomes, an example script used to run TE Density can be found at `src/TE_Density_562.sb`.

### Biased and Unbiased Gene TE Density Comparisons:
The TE density values of biased and unbiased genes were plotted using [Matplotlib](https://matplotlib.org/) and [Seaborn](https://seaborn.pydata.org/).
Graphs were made for each combination of (*Window* || *TE Type*  || *Upstream/Downstream*) for the relevant list of genes.

### Assessment of Genome-Wide TE Presence Trends in H4:
A dotplot of average TE presence for all genes was generated using `src/H4/new_h4_dotplot.py`.
The average TE density for all genes was calculated for each combination of (*Window* || *TE Type*  || *Upstream/Intragenic/Downstream*) and plotted.
