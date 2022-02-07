# Vesca_Hybrid_TE_annotations
# __author__ = Beth Alger

#user makes a list of genomes files and a list of cds files - both need to be in the same order 

#user change genome dir and working directory in 1st_run.sb and 2nd_run.sb

#check directories in mk_pangenome_1.sb

#change EDTA directory for all scripts

Order:

#first pass of EDTA. We did NOT include any previously made TE file because it resulted in low DNA TEs after the 2nd run

1st_run.sb

#remove single-copy annotations; aggregate all vesca TE libraries; remove redundants and rename TEs

mk_pangenome.sb

#run EDTA for all genomes again with new pangenome TE library

2nd_run.sb
