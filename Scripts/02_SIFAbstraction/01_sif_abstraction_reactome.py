from SPARQLWrapper import SPARQLWrapper, TURTLE, JSON, CSV
import subprocess
import time
import os
from requests.utils import requote_uri
from urllib.parse import quote
import pandas as pd
from utility import *
import glob

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data', 'BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, 'Data', 'BioPAX/ReactomeTopPathways')
results_dir = os.path.join(current_directory, 'Results/SIFAbstraction') 

os.makedirs(results_dir, exist_ok=True)

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
PREFIX abstraction:<http://abstraction/#>
"""

query_controls_state_change_of = """ 
CONSTRUCT {
    ?enzymeRef abstraction:ControlsStateChangeOf ?protRef
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?catalysis bp3:controlled ?reaction .
    ?catalysis (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .

    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?protLeft .
    ?protLeft rdf:type bp3:Protein .
    ?protLeft bp3:entityReference ?protRef .

    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?protRight .
    ?protRight rdf:type bp3:Protein .
    ?protRight bp3:entityReference ?protRef .
    ?protRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    
    # State change of the product
    ?protRight bp3:feature ?feature .
    ?feature rdf:type/(rdfs:subClassOf*) bp3:EntityFeature .
    FILTER (?enzymeRef != ?protRef)
    FILTER (?protLeft != ?protRight)
}
"""

query_transport_of = """
CONSTRUCT {
    ?enzymeRef abstraction:ControlsTransportOf ?protRef
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?protLeft .
    ?protLeft rdf:type bp3:Protein .
    ?protLeft bp3:entityReference ?protRef .

    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?protRight .
    ?protRight rdf:type bp3:Protein .
    ?protRight bp3:entityReference ?protRef .
    ?protRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    # Define cellular locations for Protein 1 and Protein 2
    ?protLeft bp3:cellularLocation ?cellularLocVocab1 .
    ?protRight bp3:cellularLocation ?cellularLocVocab2 .
    
    ?catalysis bp3:controlled ?reaction .
    ?catalysis (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER (?cellularLocVocab1 != ?cellularLocVocab2)
    FILTER (?enzymeRef != ?protRef)
    FILTER (?protLeft != ?protRight)
}
"""

query_controls_phosphorylation_of = """
CONSTRUCT {
	?enzymeRef abstraction:ControlsPhosphorylationOf ?protRef
	}
WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?catalysis bp3:controlled ?reaction .
    ?catalysis rdf:type/(rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .

    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?protLeft .
    ?protLeft rdf:type bp3:Protein .
    ?protLeft bp3:entityReference ?protRef .
    
    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?protRight .
    ?protRight rdf:type bp3:Protein .
    ?protRight bp3:entityReference ?protRef .
    ?protRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    # Specify that ?protRight is phosphorylated
    ?protRight bp3:feature ?feature .
    ?feature rdf:type bp3:ModificationFeature .
    ?feature bp3:modificationType ?modificationType .
    ?modificationType rdf:type bp3:SequenceModificationVocabulary .
    ?modificationType bp3:term ?modificationTerm .
    FILTER(CONTAINS(LCASE(?modificationTerm), "phospho"))
    FILTER(?enzymeRef != ?protRef)
    FILTER(?protLeft != ?protRight)
}
"""

query_controls_expression_of = """ 
CONSTRUCT {
    ?enzymeRef abstraction:ControlsExpressionOf ?productRef
    }
    WHERE { 
    # Case where the reaction is a Conversion
    {
        VALUES ?NucleicAcid { bp3:Dna bp3:Rna }
        ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
        ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?left .
        ?left rdf:type ?NucleicAcid .
        ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?product .
        ?product rdf:type bp3:Protein .
        ?product bp3:entityReference ?productRef .
        ?productRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        ?catalysis bp3:controlled ?reaction .
        ?catalysis rdf:type/(rdfs:subClassOf*) bp3:Control .
        ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
        ?enzyme rdf:type bp3:Protein .
        ?enzyme bp3:entityReference ?enzymeRef .
        ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    }
    UNION
    # Case where the reaction is a TemplateReaction
    {
        ?reaction rdf:type/(rdfs:subClassOf*) bp3:TemplateReaction .
        ?reaction bp3:product/(bp3:memberPhysicalEntity)* ?product .
        ?product rdf:type bp3:Protein .
        ?product bp3:entityReference ?productRef .
        ?productRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        ?catalysis bp3:controlled ?reaction .
        ?catalysis rdf:type bp3:TemplateReactionRegulation .
        ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
        ?enzyme rdf:type bp3:Protein .
        ?enzyme bp3:entityReference ?enzymeRef .
        ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    }
    FILTER (?enzymeRef != ?productRef)
}
"""

query_catalysis_precedes = """ 
CONSTRUCT {
    ?enzymeRefA abstraction:CatalysisPrecedes ?enzymeRefB
    }
    WHERE {
    ?reaction1 rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction1 bp3:left ?leftMolecule1 .
    ?leftMolecule1 rdf:type bp3:SmallMolecule .
    # get the connecting molecule
    ?reaction1 bp3:right ?connectingMolecule .
    ?connectingMolecule rdf:type bp3:SmallMolecule .
    ?connectingMolecule bp3:entityReference ?connectingMoleculeRef .
    ?connectingMoleculeRef bp3:xref ?connectingMoleculeXref .
    ?connectingMoleculeXref rdf:type bp3:UnificationXref .
    ?connectingMoleculeXref bp3:db "ChEBI" .
    ?connectingMoleculeXref bp3:id ?connectingMoleculeID .
    FILTER (?connectingMoleculeID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))
    ?reaction2 bp3:left ?connectingMolecule .
    ?reaction2 rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction2 bp3:right ?rightMolecule2 .
    ?rightMolecule2 rdf:type bp3:SmallMolecule .

    # Get catalysis of the first reaction
    ?catalysis1 bp3:controlled ?reaction1 .
    ?catalysis1 (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis1 bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzymeProtA .
    ?enzymeProtA rdf:type bp3:Protein .
    ?enzymeProtA bp3:entityReference ?enzymeRefA .
    ?enzymeRefA rdf:type/(rdfs:subClassOf*) bp3:EntityReference .

    # get catalysis of the second reaction
    ?catalysis2 bp3:controlled ?reaction2 .
    ?catalysis2 (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis2 bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzymeProtB .
    ?enzymeProtB rdf:type bp3:Protein .
    ?enzymeProtB bp3:entityReference ?enzymeRefB .
    ?enzymeRefB rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    
    FILTER (?enzymeRefA != ?enzymeRefB)
    FILTER (?leftMolecule1 != ?rightMolecule2)
    FILTER (?leftMolecule1 != ?connectingMolecule)
    FILTER (?rightMolecule2 != ?connectingMolecule)
}
"""

query_in_complex_with = """ 
CONSTRUCT {
    ?prot1Ref abstraction:InComplexWith ?prot2Ref
}
WHERE {
    # Case 1: Complex with several different proteins
    {
        ?complex rdf:type/(rdfs:subClassOf*) bp3:Complex .
        ?complex (bp3:component|bp3:memberPhysicalEntity)* ?prot1 .
        ?complex (bp3:component|bp3:memberPhysicalEntity)* ?prot2 .
        ?prot1 rdf:type bp3:Protein .
        ?prot2 rdf:type bp3:Protein .
        ?prot1 bp3:entityReference ?prot1Ref .
        ?prot2 bp3:entityReference ?prot2Ref .
        ?prot1Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        ?prot2Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        # Avoid self-pairing and duplicate pairs
        FILTER (STR(?prot1Ref) < STR(?prot2Ref))
    }
    UNION
    # Case 2: Same protein with stoichiometry > 1
    {
        ?complex rdf:type/(rdfs:subClassOf*) bp3:Complex .
        # Check for stoichiometry > 1
        ?complex bp3:componentStoichiometry ?stoich .
        ?stoich bp3:physicalEntity ?entity1 .
        ?stoich bp3:stoichiometricCoefficient ?coeff .
        FILTER (?coeff > 1)
        ?complex (bp3:component|bp3:memberPhysicalEntity)* ?prot1 .
        ?prot1 rdf:type bp3:Protein .
        ?prot1 bp3:entityReference ?prot1Ref .
        ?prot1Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        # For self-pairing, use the same entity
        BIND (?prot1Ref AS ?prot2Ref)
    }
    FILTER (?prot1Ref != ?prot2Ref)
}
"""

query_interacts_with = """  
CONSTRUCT {
    ?participant1Ref abstraction:InteractsWith ?participant2Ref
    }
    WHERE {
    VALUES ?participantType { bp3:Protein bp3:SmallMolecule }
    ?MolecularInteraction rdf:type/(rdfs:subClassOf*) bp3:MolecularInteraction .
    ?MolecularInteraction bp3:participant/(bp3:memberPhysicalEntity)* ?participant1 .
    ?participant1 rdf:type ?participantType .
    ?participant1 bp3:entityReference ?participant1Ref .
    ?participant1Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    ?MolecularInteraction bp3:participant/(bp3:memberPhysicalEntity)* ?participant2 .
    ?participant2 rdf:type ?participantType .
    ?participant2 bp3:entityReference ?participant2Ref .
    ?participant2Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER (STR(?participant1Ref) < STR(?participant2Ref))
    FILTER (?participant1Ref != ?participant2Ref)
}
"""

query_neighbor_of = """ 
CONSTRUCT {
    ?participant1Ref abstraction:NeighborOf ?participant2Ref
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .

    ?reaction ((^bp3:controlled/bp3:controller)|bp3:left|bp3:right) ?participant1 .
    ?participant1 (bp3:component|bp3:memberPhysicalEntity)* ?participant1Prot .
    ?participant1Prot rdf:type bp3:Protein .
    ?participant1Prot bp3:entityReference ?participant1Ref .
    ?participant1Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .

    ?reaction ((^bp3:controlled/bp3:controller)|bp3:left|bp3:right) ?participant2 .
    ?participant2 (bp3:component|bp3:memberPhysicalEntity)* ?participant2Prot .
    ?participant2Prot rdf:type bp3:Protein .
    ?participant2Prot bp3:entityReference ?participant2Ref .
    ?participant2Ref rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER (STR(?participant1Ref) < STR(?participant2Ref))
    FILTER (?participant1Ref != ?participant2Ref)
}
"""

query_consumption_controled_by = """  
CONSTRUCT {
    ?reactantRef abstraction:ConsumptionControlledBy ?enzymeRef
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?reactant .
    ?reactant rdf:type bp3:SmallMolecule .
    ?reactant bp3:entityReference ?reactantRef .
    ?reactantRef rdf:type bp3:SmallMoleculeReference .
    ?reactantRef bp3:xref ?reactantRefXref .
    ?reactantRefXref rdf:type bp3:UnificationXref .
    ?reactantRefXref bp3:db "ChEBI" .
    ?reactantRefXref bp3:id ?reactantID .

    FILTER (?reactantID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?product .
    ?product rdf:type bp3:SmallMolecule .
    ?product bp3:entityReference ?productRef .
    ?productRef rdf:type bp3:SmallMoleculeReference .
    ?productRef bp3:xref ?productRefXref .
    ?productRefXref rdf:type bp3:UnificationXref .
    ?productRefXref bp3:db "ChEBI" .
    ?productRefXref bp3:id ?productID .

    FILTER (?productID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))
    
    ?catalysis bp3:controlled ?reaction .
    ?catalysis (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER (?reactantRef != ?enzymeRef)
    FILTER (?reactant != ?product)
}
"""

query_controls_production = """ 
CONSTRUCT {
    ?enzymeRef abstraction:ControlsProductionOf ?productRef
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?reactant .
    ?reactant rdf:type bp3:SmallMolecule .
    ?reactant bp3:entityReference ?reactantRef .
    ?reactantRef rdf:type bp3:SmallMoleculeReference .
    ?reactantRef bp3:xref ?reactantRefXref .
    ?reactantRefXref rdf:type bp3:UnificationXref .
    ?reactantRefXref bp3:db "ChEBI" .
    ?reactantRefXref bp3:id ?reactantID .

    FILTER (?reactantID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))
    
    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?product .
    ?product rdf:type bp3:SmallMolecule .
    ?product bp3:entityReference ?productRef .
    ?productRef rdf:type bp3:SmallMoleculeReference .
    ?productRef bp3:xref ?productRefXref .
    ?productRefXref rdf:type bp3:UnificationXref .
    ?productRefXref bp3:db "ChEBI" .
    ?productRefXref bp3:id ?productID .

    FILTER (?productID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))
    
    ?catalysis bp3:controlled ?reaction .
    ?catalysis (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER(?enzymeRef != ?productRef)
    FILTER(?reactant != ?product)
}
"""

query_controls_transport_of_chemical = """ 
CONSTRUCT {
    ?enzymeRef abstraction:ControlsTransportOfChemical ?SmallMoleculeRef
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?SmallMolecule1 .
    ?SmallMolecule1 rdf:type bp3:SmallMolecule .
    ?SmallMolecule1 bp3:cellularLocation ?cellularLocVocab1 .
    ?SmallMolecule1 bp3:entityReference ?SmallMoleculeRef .
    ?SmallMoleculeRef bp3:xref ?SmallMoleculeRefXref .
    ?SmallMoleculeRefXref rdf:type bp3:UnificationXref .
    ?SmallMoleculeRefXref bp3:db "ChEBI" .
    ?SmallMoleculeRefXref bp3:id ?SmallMoleculeID .

    FILTER (?SmallMoleculeID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?SmallMolecule2 . 
    ?SmallMolecule2 rdf:type bp3:SmallMolecule .
    ?SmallMolecule2 bp3:cellularLocation ?cellularLocVocab2 .
    ?SmallMolecule2 bp3:entityReference ?SmallMoleculeRef .
    ?SmallMoleculeRef bp3:xref ?SmallMoleculeRefXref .
    ?SmallMoleculeRefXref rdf:type bp3:UnificationXref .
    ?SmallMoleculeRefXref bp3:db "ChEBI" .
    ?SmallMoleculeRefXref bp3:id ?SmallMoleculeID .

    ?catalysis bp3:controlled ?reaction .
    ?catalysis (rdf:type/rdfs:subClassOf*) bp3:Control .
    ?catalysis bp3:controller/((bp3:component|bp3:memberPhysicalEntity)*) ?enzyme .
    ?enzyme rdf:type bp3:Protein .
    ?enzyme bp3:entityReference ?enzymeRef .
    ?enzymeRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
    FILTER(?cellularLocVocab1 != ?cellularLocVocab2)
}
"""

query_chemical_affects = """  
CONSTRUCT {
    ?chemicalCatalyzerRef abstraction:ChemicalAffects ?proteinRef .
    }
    WHERE {
    # the small molecule has an effect on the protein state through binding
    {
        ?complex rdf:type bp3:Complex .
        ?complex (bp3:component|bp3:memberPhysicalEntity)* ?proteinComponent .
        ?proteinComponent rdf:type bp3:Protein .
        # modification of the protein
        ?proteinComponent bp3:feature ?feature .
        ?feature rdf:type/(rdfs:subClassOf*) bp3:EntityFeature .
        ?complex (bp3:component|bp3:memberPhysicalEntity)* ?chemicalCatalyzer .
        ?chemicalCatalyzer rdf:type bp3:SmallMolecule .
        ?chemicalCatalyzer bp3:entityReference ?chemicalCatalyzerRef .
        ?chemicalCatalyzerRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        ?chemicalCatalyzerRef bp3:xref ?chemicalCatalyzerRefXref .
        ?chemicalCatalyzerRefXref rdf:type bp3:UnificationXref .
        ?chemicalCatalyzerRefXref bp3:db "ChEBI" .
        ?chemicalCatalyzerRefXref bp3:id ?chemicalCatalyzerID .
        FILTER (?chemicalCatalyzerID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))
        BIND(?proteinComponent AS ?protein1)
    }
    UNION
    # the small molecule has an effect on the protein state by controlling a Conversion reaction
    {
        ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
        ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?protLeft .
        ?protLeft rdf:type bp3:Protein .
        ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?protRight .
        ?protRight rdf:type bp3:Protein .

        ?catalysis bp3:controlled ?reaction .
        ?catalysis rdf:type/(rdfs:subClassOf*) bp3:Control .
        ?catalysis bp3:controller ?chemicalCatalyzer .
        ?chemicalCatalyzer rdf:type bp3:SmallMolecule .
        ?chemicalCatalyzer bp3:entityReference ?chemicalCatalyzerRef .
        ?chemicalCatalyzerRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
        ?chemicalCatalyzerRef bp3:xref ?chemicalCatalyzerRefXref .
        ?chemicalCatalyzerRefXref rdf:type bp3:UnificationXref .
        ?chemicalCatalyzerRefXref bp3:db "ChEBI" .
        ?chemicalCatalyzerRefXref bp3:id ?chemicalCatalyzerID .
        FILTER (?chemicalCatalyzerID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

        # modification of the protein
        ?protRight bp3:feature ?feature .
        ?feature rdf:type/(rdfs:subClassOf*) bp3:EntityFeature .
        BIND(?protLeft AS ?protein2)
    }
    BIND(COALESCE(?protein1, ?protein2) AS ?protein)
    ?protein bp3:entityReference ?proteinRef .
    ?proteinRef rdf:type/(rdfs:subClassOf*) bp3:EntityReference .
}
"""

query_reacts_with = """ 
CONSTRUCT {
    ?smallMolecule1Ref abstraction:ReactsWith ?smallMolecule2Ref .
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:BiochemicalReaction .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?smallMolecule1 .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?smallMolecule2 .
    ?smallMolecule1 rdf:type bp3:SmallMolecule .
    ?smallMolecule1 bp3:entityReference ?smallMolecule1Ref .
    ?smallMolecule1Ref rdf:type bp3:SmallMoleculeReference .
    ?smallMolecule1Ref bp3:xref ?smallMolecule1RefXref .
    ?smallMolecule1RefXref rdf:type bp3:UnificationXref .
    ?smallMolecule1RefXref bp3:db "ChEBI" .
    ?smallMolecule1RefXref bp3:id ?smallMolecule1ID .
    FILTER (?smallMolecule1ID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    ?smallMolecule2 rdf:type bp3:SmallMolecule .
    ?smallMolecule2 bp3:entityReference ?smallMolecule2Ref .
    ?smallMolecule2Ref rdf:type bp3:SmallMoleculeReference .
    ?smallMolecule2Ref bp3:xref ?smallMolecule2RefXref .
    ?smallMolecule2RefXref rdf:type bp3:UnificationXref .
    ?smallMolecule2RefXref bp3:db "ChEBI" .
    ?smallMolecule2RefXref bp3:id ?smallMolecule2ID .
    FILTER (?smallMolecule2ID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    FILTER (?smallMolecule1Ref != ?smallMolecule2Ref)
}
"""

query_used_to_produce = """ 
CONSTRUCT {
    ?smallMolecule1Ref abstraction:UsedToProduce ?smallMolecule2Ref .
    }
    WHERE {
    ?reaction rdf:type/(rdfs:subClassOf*) bp3:Conversion .
    ?reaction bp3:left/(bp3:memberPhysicalEntity)* ?smallMolecule1 .
    ?reaction bp3:right/(bp3:memberPhysicalEntity)* ?smallMolecule2 .
    ?smallMolecule1 rdf:type bp3:SmallMolecule .
    ?smallMolecule1 bp3:entityReference ?smallMolecule1Ref .
    ?smallMolecule1Ref rdf:type bp3:SmallMoleculeReference .
    ?smallMolecule1Ref bp3:xref ?smallMolecule1RefXref .
    ?smallMolecule1RefXref rdf:type bp3:UnificationXref .
    ?smallMolecule1RefXref bp3:db "ChEBI" .
    ?smallMolecule1RefXref bp3:id ?smallMolecule1ID .
    FILTER (?smallMolecule1ID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    ?smallMolecule2 rdf:type bp3:SmallMolecule .
    ?smallMolecule2 bp3:entityReference ?smallMolecule2Ref .
    ?smallMolecule2Ref rdf:type bp3:SmallMoleculeReference .
    ?smallMolecule2Ref bp3:xref ?smallMolecule2RefXref .
    ?smallMolecule2RefXref rdf:type bp3:UnificationXref .
    ?smallMolecule2RefXref bp3:db "ChEBI" .
    ?smallMolecule2RefXref bp3:id ?smallMolecule2ID .
    FILTER (?smallMolecule2ID NOT IN ("CHEBI:15377", "CHEBI:24636", "CHEBI:15378", "CHEBI:15379", "CHEBI:57783",
        "CHEBI:13392", "CHEBI:16474", "CHEBI:58349", "CHEBI:18009", "CHEBI:13390",
        "CHEBI:25523", "CHEBI:25524", "CHEBI:30616", "CHEBI:15422", "CHEBI:16761",
        "CHEBI:456216", "CHEBI:57945", "CHEBI:16908", "CHEBI:57540", "CHEBI:15846",
        "CHEBI:16526", "CHEBI:29888", "CHEBI:35782", "CHEBI:68836", "CHEBI:18361",
        "CHEBI:43474", "CHEBI:35780", "CHEBI:18367", "CHEBI:26078", "CHEBI:77740",
        "CHEBI:28931", "CHEBI:58307", "CHEBI:17877", "CHEBI:17552", "CHEBI:58189",
        "CHEBI:37565", "CHEBI:15996"))

    FILTER (?smallMolecule1Ref != ?smallMolecule2Ref) 
}
"""

query_entity_refs="""
SELECT DISTINCT ?entityRef ?entityRefName ?entityID ?entityRefType
WHERE {
    VALUES ?entityRefType { bp3:ProteinReference bp3:SmallMoleculeReference }
    ?entityRef rdf:type ?entityRefType .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref bp3:id ?entityID .
    ?entityRef bp3:name ?entityRefName .
}
"""

endpoint = "http://localhost:3030/top_pathway"
counter = 1
filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))

for owl_file in sorted(filelist):
    # creating results directory
    top_pathway_results_dir = os.path.join(results_dir, f"{counter:02d}_SIF")
    os.makedirs(top_pathway_results_dir, exist_ok=True)

    # launch Fuseki SPARQL endpoint for top pathway
    ################################################################################
    print(f"Processing file: {owl_file}")
    command = [
        '/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server',
        '--file', owl_file,
        '--file', BioPAX_Ontology_file_path,
        '/top_pathway'
    ]
    print("Fuseki command:", command)
    process = subprocess.Popen(command)
    time.sleep(30)
    ##################################################################################

    # 01-controls-state-change-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_controls_state_change_of)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-01-ControlsStateChangeOf.csv"), format='csv')

    # 02-controls-transport-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_transport_of)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-02-ControlsTransportOf.csv"), format='csv')

    # 03-controls-phosphorylation-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_controls_phosphorylation_of)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-03-ControlsPhosphorylationOf.csv"), format='csv')

    # 04-controls-expression-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_controls_expression_of)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-04-ControlsExpressionOf.csv"), format='csv')

    # 05-catalysis-precedes
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_catalysis_precedes)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-05-CatalysisPrecedes.csv"), format='csv')

    # 06-in-complex-with
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_in_complex_with)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-06-InComplexWith.csv"), format='csv')

    # 07-interacts-with
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_interacts_with)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-07-InteractsWith.csv"), format='csv')

    # 08-neighbor-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_neighbor_of)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-08-NeighborOf.csv"), format='csv')

    # 09-consumption-controlled-by
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_consumption_controled_by)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-09-ConsumptionControlledBy.csv"), format='csv')

    # 10-controls-production-of
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_controls_production)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-10-ControlsProductionOf.csv"), format='csv')

    # 11-controls-transport-of-chemical
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_controls_transport_of_chemical)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-11-ControlsTransportOfChemical.csv"), format='csv')

    # 12-chemical-affects
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_chemical_affects)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-12-ChemicalAffects.csv"), format='csv')

    # 13-reacts-with
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_reacts_with)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-13-ReactsWith.csv"), format='csv')

    # 14-used-to-produce
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_used_to_produce)
    df = save_for_cytoscape(sparql, prefixes, os.path.join(top_pathway_results_dir, f"{counter:02d}-14-UsedToProduce.csv"), format='csv')

    # aggregation of files in a single SIF file
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes+query_entity_refs)
    sparql.setReturnFormat(CSV)
    results = sparql.query().convert()
    with open(os.path.join(results_dir, f"UtilityFiles/{counter:02d}-EntityRefsIDs.csv"), "wb") as f:
        f.write(results)

    results_files = list()
    for file in os.listdir(top_pathway_results_dir):
        results_files.append(file)
    aggregatedNetwork = pd.concat([pd.read_csv(os.path.join(top_pathway_results_dir, file), header=0, sep=",") for file in results_files]) 
    print(aggregatedNetwork.head()) 

    identifier_mapping_file = pd.read_csv(os.path.join(results_dir, f"UtilityFiles/{counter:02d}-EntityRefsIDs.csv"), sep=",", header=0)
    dico_entity_ref = dict()
    for index, row in identifier_mapping_file.iterrows():
        parts = row[0].split("#")
        ref = "reactome:" + parts[-1]
        #ref = row[0].replace("http://www.reactome.org/biopax/96/9612973#", "reactome:")
        if not ref in dico_entity_ref.keys():
            dico_entity_ref[ref] = row[2]

    aggregatedNetworkIDS = pd.DataFrame(columns=["Source", "Interaction", "Target"])
    i = 0
    for index, row in aggregatedNetwork.iterrows():
        #entity1 = row[0]
        #entity2 = row[2]
        parts1 = row[0].split("#")
        entity1 = "reactome:" + parts1[-1]
        parts2 = row[2].split("#")
        entity2 = "reactome:" + parts2[-1]
        interaction = row[1]
        entity1ID = dico_entity_ref[entity1]
        entity2ID = dico_entity_ref[entity2]
        aggregatedNetworkIDS.at[i, "Source"] = entity1ID
        aggregatedNetworkIDS.at[i, "Interaction"] = interaction
        aggregatedNetworkIDS.at[i, "Target"] = entity2ID
        i += 1
    aggregatedNetworkIDS.to_csv(os.path.join(top_pathway_results_dir, f"{counter:02d}-SIF-abstraction.sif"), sep=" ", header=0, index=None)

    process.kill()
    time.sleep(30)
    counter += 1