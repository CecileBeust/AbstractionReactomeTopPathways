#!/bin/bash
set -e  # Arrête le script si une commande échoue

# generate pathway-centered abstraction on Reactome top pathways
echo "[1/3] Running pathway abstraction..."
python3 01_pathway_abstraction_reactome.py

# get uniprot ids per pathway
#echo "[2/3] Getting UniProt IDs per pathway..."
#python3 02_uniprot_per_pathways.py

# weight pathway-centered abstraction
echo "[3/3] Weighting pathway abstraction..."
python3 03_weight_pathway_abstraction.py

# clean intermediate files (ONLY AFTER ALL SCRIPTS HAVE RUN)
echo "Cleaning intermediate files..."
rm -f ../../Results/PathwayAbstraction/01_TopPathways/*.tsv
rm -f ../../Results/PathwayAbstraction/01_TopPathways/*IsAComponentOf.csv
rm -f ../../Results/PathwayAbstraction/01_TopPathways/*NextStepPathway.csv
rm -f ../../Results/PathwayAbstraction/01_TopPathways/*IsAComponentOf_ResnikER.csv
rm -f ../../Results/PathwayAbstraction/01_TopPathways/*NextStepPathway_ERcontent.csv

echo "Done."