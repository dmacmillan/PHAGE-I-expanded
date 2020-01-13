# PHAGE-I-expanded (Proportion of HLA Associated Genomic Escape)

## Purpose

## Input

## Output

## Methods
  * Input:
    * A list HLA-associated polymorphisms
    * A list of patients, their HLA types, and HIV sequences for a given protein
    * An HIV protein (selected from a list)
  * Algorithm:
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