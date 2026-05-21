#!/bin/sh

# generate pathway-centered abstraction on Reactome top pathways
python 01_pathway_abstraction_reactome.py

# get uniprot ids per pathway
python 02_uniprot_per_pathways.py

# weight pathway-centered abstraction
python 03_weight_pathway_abstraction.py

# clean intermediate files
rm ../../Results/PathwayAbstraction/*.tsv
rm ../../Results/PathwayAbstraction/*IsAComponentOf.csv
rm ../../Results/PathwayAbstraction/*NextStepPathway.csv
rm ../../Results/PathwayAbstraction/*IsAComponentOf_ResnikER.csv
rm ../../Results/PathwayAbstraction/*NextStepPathway_ERcontent.csv