## Abstraction of Reactome Top Pathways in BioPAX

This repository contains code to abstract the 29 top pathways of Reactome using different contextual abstraction methods:
- **Pathway-centered abstraction**: abstraction of Reactome pathways extracting the hierarchy and sequence of pathways from the BioPAX graph
    - Scripts: `Scripts/01_PathwayAbstraction` 
    - Abstractions: `Results/PathwayAbstraction/` 
- **SIF (Simple Interaction Format) abstraction**: protein and small molecule-centered abstraction of Reactome pathways 
    - Scripts: `Scripts/02_SIFAbstraction/` 
    - Abstractions: `Results/SIFAbstraction`
- **Revisited pathway co-membership abstraction**: protein-centered abstraction that generates a weighted clique of proteins belonging to the same pathway in Reactome. Co-membership interactions are scored with regard to the biological distance of pathways to which proteins belong.
    - Scripts: `Scripts/03_PathwayComembership` 
    - Results: `Results/PathwayComembership`