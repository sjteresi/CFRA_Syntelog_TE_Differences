# scripts for running calculating TE density for various strawberry genomes
# and for performing an analysis of the TE density results
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(realpath $(ROOT_DIR)/data/)
DEV_RESULTS := $(realpath $(ROOT_DIR)/results/)
DEV_562_RESULTS := $(DEV_RESULTS)/graphs/biased_562
DEV_1008_RESULTS := $(DEV_RESULTS)/graphs/biased_1008
DEV_UNBIASED_RESULTS := $(DEV_RESULTS)/graphs/unbiased
DEV_562_GENES := $(DEV_DATA)/maker_annotation.562.gff_sorted.gff
DEV_1008_GENES := $(DEV_DATA)/2339.maker_annotation.gff_sorted.gff
DEV_562_TEs := $(DEV_DATA)/562_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_1008_TEs := $(DEV_DATA)/1008_2339_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_RAW_SYNTELOGS := $(DEV_DATA)/CFRA2339_CFRA562_Syntelogs.csv
DEV_SYNTELOGS := $(DEV_RESULTS)/set_syntelogs.tsv
DEV_562_DENSITY_DATA := $(DEV_RESULTS)/finalized_data/562_TE_Density
DEV_1008_DENSITY_DATA := $(DEV_RESULTS)/finalized_data/1008_2339_TE_Density
DEV_562_GENE_DATA := $(DEV_RESULTS)/filtered_input_data/Cleaned_562_Genes.tsv
DEV_1008_GENE_DATA := $(DEV_RESULTS)/filtered_input_data/Cleaned_2339_Genes.tsv
DEV_562_BIASED_GENES_10X :=$(DEV_DATA)/Fv562_Biased_10x.txt
DEV_1008_BIASED_GENES_10X :=$(DEV_DATA)/Fv2339_Biased_10x.txt
DEV_562_BIASED_GENES_5X :=$(DEV_DATA)/Fv562_5x_Bias.txt
DEV_1008_BIASED_GENES_5X :=$(DEV_DATA)/Fv2339_5x_Bias.txt
DEV_562_BIASED_GENES_2X :=$(DEV_DATA)/Fv562_2x_Bias.txt
DEV_1008_BIASED_GENES_2X :=$(DEV_DATA)/Fv2339_2x_Bias.txt
DEV_PROCESSED_SYNTELOGS := $(DEV_RESULTS)/Syntelog_Matches.tsv
DEV_UNBIASED_GENES := $(DEV_DATA)/Fv2339_n_Fv562_Unbiased.txt

# H4 Related Analyses
DEV_H4_GENES := $(DEV_DATA)/H4/H4_Genes.gtf
DEV_H4_TEs := $(DEV_DATA)/H4/F_vesca_H4_V4.1.fasta.mod.EDTA.TEanno.gff3
DEV_H4_GENE_DATA := $(DEV_RESULTS)/filtered_input_data/Cleaned_H4_Genes.tsv
DEV_H4_DENSITY_DATA := $(DEV_RESULTS)/finalized_data/H4_TE_Density
DEV_H4_RESULTS := $(DEV_RESULTS)/graphs/H4




filter_genes:
	@echo Filtering 562 genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_562_GENES) 562 --o $(ROOT_DIR)/results/filtered_input_data
	@echo 
	@echo Filtering 1008/2339 genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_1008_GENES) 2339 --o $(ROOT_DIR)/results/filtered_input_data

1008_2339_filter_TEs:
	@echo Filtering 1008/2339 TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_1008_TEs) 2339 --o $(ROOT_DIR)/results/filtered_input_data

562_filter_TEs:
	@echo Filtering 562 TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_562_TEs) 562 --o $(ROOT_DIR)/results/filtered_input_data

calculate_TE_Density:
	@echo Running TE Density for 562
	sbatch $(ROOT_DIR)/src/TE_Density_562.sb
	@echo

filter_syntelogs:
	@echo reading raw syntelog file and fixing names
	mkdir -p $(DEV_RESULTS)
	python $(ROOT_DIR)/src/import_syntelogs.py $(DEV_RAW_SYNTELOGS) set_syntelogs.tsv $(DEV_RESULTS)

add_syntelog_HDF5_indices:
	@echo adding indices of TE density HDF5 indices to syntelog table
	mkdir -p $(DEV_RESULTS)
	python $(ROOT_DIR)/src/add_indices_to_syntelog_table.py $(DEV_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) -o $(DEV_RESULTS)

generate_graphs_562:
	@echo Generating TE density syntelog comparison graphs
	mkdir -p $(DEV_562_RESULTS)
	@echo Generating graphs for 562 at 2X
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_562_BIASED_GENES_2X) 562 2 -o $(DEV_562_RESULTS)
	@echo Generating graphs for 562 at 5X
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_562_BIASED_GENES_5X) 562 5 -o $(DEV_562_RESULTS)
	@echo Generating graphs for 562 at 10X 
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_562_BIASED_GENES_10X) 562 10 -o $(DEV_562_RESULTS)

generate_graphs_1008:
	# Functionally equivalent to generate_graphs_562 except for the use of input files
	@echo Generating TE density syntelog comparison graphs
	mkdir -p $(DEV_1008_RESULTS)
	@echo Generating graphs for 1008 at 2X
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_1008_BIASED_GENES_2X) 1008 2 -o $(DEV_1008_RESULTS)
	@echo Generating graphs for 1008 at 5X
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_1008_BIASED_GENES_5X) 1008 5 -o $(DEV_1008_RESULTS)
	@echo Generating graphs for 1008 at 10X 
	python $(ROOT_DIR)/src/compare_biased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_1008_BIASED_GENES_10X) 1008 10 -o $(DEV_1008_RESULTS)

generate_unbiased_graphs:
	@echo Generating TE density syntelog comparison graphs
	@echo 
	mkdir -p $(DEV_UNBIASED_RESULTS)
	python $(ROOT_DIR)/src/compare_unbiased_exp_w_density.py $(DEV_PROCESSED_SYNTELOGS) $(DEV_562_DENSITY_DATA) $(DEV_1008_DENSITY_DATA) $(DEV_562_GENE_DATA) $(DEV_1008_GENE_DATA) $(DEV_UNBIASED_GENES) -o $(DEV_UNBIASED_RESULTS)

# Code relevant to H4
H4_filter_TEs:
	@echo Filtering H4 TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_H4_TEs) H4 --o $(ROOT_DIR)/results/filtered_input_data

H4_filter_genes:
	@echo Filtering H4 genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/H4/import_strawberry_H4_gene_anno.py $(DEV_H4_GENES) H4 --o $(ROOT_DIR)/results/filtered_input_data

H4_dotplot:
	@echo Generating dotplot for H4
	mkdir -p $(DEV_H4_RESULTS)
	python $(ROOT_DIR)/src/H4/new_h4_dotplot.py $(DEV_H4_DENSITY_DATA) $(DEV_H4_GENE_DATA) -o $(DEV_H4_RESULTS)

dev_562_dotplot:
	python $(ROOT_DIR)/src/H4/new_h4_dotplot.py $(DEV_562_DENSITY_DATA) $(DEV_562_GENE_DATA) -o $(DEV_H4_RESULTS)
