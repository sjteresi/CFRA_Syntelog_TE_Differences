# scripts for running calculating TE density for various strawberry genomes
# and for performing an analysis of the TE density results
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(realpath $(ROOT_DIR)/data/)
DEV_562_GENES := $(DEV_DATA)/maker_annotation.562.gff_sorted.gff
DEV_1008_GENES := $(DEV_DATA)/2339.maker_annotation.gff_sorted.gff
DEV_562_TEs := $(DEV_DATA)/562_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_1008_TEs := $(DEV_DATA)/1008.fasta.mod.EDTA.TEanno.gff3
#DEV_TES := $(DEV_DATA)/TEs/Blueberry_EDTA_TEs.gff
#DEV_GENE_EXPRESSION := $(DEV_DATA)/Genes/Blueberry_TPM_All.tsv
#DEV_DENSITY_FILES := $(DEV_DATA)/../../TE_Density_Example_Data/Blueberry/TEs/Density_Data/Draper
#DEV_GENOME := "BlueberryExample"


filter_genes:
	@echo Filtering 562 genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_562_GENES) 562
	@echo 
	@echo Filtering 1008/2339 genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_gene_anno.py $(DEV_1008_GENES) 2339

filter_TEs:
	@echo Filtering 562 TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_562_TEs) 562
	@echo 
	@echo Filtering 1008/2339 TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_strawberry_EDTA.py $(DEV_1008_TEs) 2339

calculate_TE_Density:
	@echo Running TE Density for 562
	sbatch $(ROOT_DIR)/src/TE_Density_562.sb
	@echo
#
#generate_expression_graphs:
#	@echo Generating TE density vs. gene expression graphs
#	@echo 
#	mkdir -p $(ROOT_DIR)/results/graphs
#	python $(ROOT_DIR)/src/compare_expression.py $(DEV_GENE_EXPRESSION) $(DEV_GENES) $(DEV_DENSITY_FILES) -o $(ROOT_DIR)/results/graphs
