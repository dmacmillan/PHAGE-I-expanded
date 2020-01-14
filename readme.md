# PHAGE-I-expanded (Proportion of HLA Associated Genomic Escape)

## Purpose

## Input
* A list HLA-associated polymorphisms
* A list of patients, their HLA types, and HIV sequences for a given protein
* An HIV protein (selected from a list)

## Output
Tab-separated values with the following fields/columns:
* `patient_ID` - A unique identifier for a patient
* `HIV_protein` - The HIV protein being analyzed
* `HLA_allele` - The HLA allele that the patient has that is associated with the HIV codon in the following column
* `HIV_CODON` - The HIV codon position associated with the HLA in question (HXB2 coordinates)
* `patient_AA` - The amino acid of the patient at this codon
* `state` - Defines whether the patient's sequence is adapted, non-adapted, possible-adapted, or possible-nonadapted
* `CTL_epitope` - The actual epitope amino acid sequence if the HLA-associated site in question is within 3 amino acids of an optimal CTL epitope (most of these come from LANL)
* `HLA_restriction` - The HLA restriction defined for that epitope, sometimes this matches exactly to the HLA of the patient, and sometimes the epitope was only defined at the two digit resolution level. Other times the epitope is restricted by a closely related HLA allele to the one the patient expresses
* `epitope_coordinates` - The amino acid coordinates of the epitope
* `epitope_source` - Where the epitope was originally defined
* `expanded_HLA_definition` - Whether or not an expanded HLA matching occurred to match the patient HLA (matched a closely related HLA)
* `epitope_position` - Position within the epitope where the HLA-associated site is located. E.g. "9C" means position 9 which is also the C-terminus of the epitope. "C+1" means one amino acid downstream of the C-terminus. "N-1" means one amino acid upsteam of the N-terminus.

## Methods
* Translate patient sequences into amino acids by synonymously resolving mixture nucleotides and store patient data in a hash
* Group input HLA polymorphism list into a hash for quick access
* For each association compare against each patient
* For each matching association, determine if the association is adapted, non-adapted, possible adapted, or possible non-adapted
  * `adapted`: The amino acid of the patient matches exactly the amino acid and position of the association and the association is adapted
  * `non-adapted`: The amino acid of the patient matches exactly the amino acid and position of the association and the association is non-adapted
  * `possible adapted`: Two ways:
    * Both an adapated and non-adapted association exist and but neither match the patient's amino acid at that position
    * Only a non-adapted association exists but it does not match the patient's amino acid at that position
  * `possible non-adapted`: Only adapted associations exist at the given location but none of them match the patient's amino acid
* For each epitope defined in the epitope database for the given protein determine if the patient has an association within +/- 3 nucleotides and if so, record the epitopes and their sources
* If no epitope found, search all supertype expansions for every epitope and try to find if the patient has an association. If so, record the epitope and mark it as "expanded definition"