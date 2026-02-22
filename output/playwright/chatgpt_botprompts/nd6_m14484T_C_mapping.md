# Protein-level consequence of the mtDNA variant m.14484T>C in human MT-ND6

## Executive summary

Using the revised Cambridge Reference Sequence (**rCRS**) coordinate system implemented as RefSeq **NC_012920.1**, the mitochondrial DNA variant **m.14484T>C** occurs in **MT‑ND6** and is consistently annotated by authoritative curated resources as encoding the missense substitution **p.Met64Val** (one‑letter: **M64V**; three‑letter: **Met→Val**). citeturn29view0turn30view3turn2view0turn8view0  
ClinVar classifies the variant as **Pathogenic** for **Leber hereditary optic neuropathy (LHON)** (listed as “Leber optic atrophy”), with multiple submissions and recent evaluation updates (e.g., Oct 22, 2025 for the LHON condition record). citeturn1view1  

## Reference coordinate system and reference sequences used

The variant nomenclature **m.14484T>C** is interpreted relative to the **rCRS**, which MITOMAP states is the GenBank/RefSeq sequence **NC_012920** (and commonly used as the standard comparison sequence for human mtDNA variant reporting). citeturn29view0turn29view1  
ClinVar’s preferred name for the variant explicitly uses the same reference sequence context: **NC_012920.1(MT‑ND6):m.14484T>C**. citeturn30view3  

## Mapping m.14484T>C to the MT-ND6 codon and translating it

### Gene context needed for codon mapping

MITOMAP’s MT‑ND6 locus description specifies that **mtND6 spans 14149–14673** and is **transcribed off the reverse strand, beginning at 14673**. citeturn27view0  
This strand/origin information is sufficient to map the coordinate **14484** into the MT‑ND6 reading frame.

### Codon position within MT-ND6

If MT‑ND6 translation starts at **14673** and proceeds in codons of three nucleotides along the reverse-strand transcript, then the offset from the start to position 14484 is:

- Offset = 14673 − 14484 = 189 nucleotides  
- 189 / 3 = 63 full codons before this base  
- Therefore, **14484 corresponds to the first base of codon 64**, i.e., it affects **amino acid position 64**.  

This positional conclusion is concordant with the curated protein consequence reported by ClinVar and UniProt (both place the effect at residue 64). citeturn27view0turn31view0turn8view0  

### Nucleotide-to-codon-to-amino-acid step

Authoritative variant records state the reference allele and alternate allele at this coordinate as **T→C** (m.14484T>C). citeturn30view3turn2view0  
Because MT‑ND6 is on the reverse strand (MITOMAP), this corresponds to a **complementary base change in the coding direction**, i.e., **A→G** at the first position of the affected codon. citeturn27view0  

Using the **vertebrate mitochondrial genetic code (NCBI transl_table=2)**, the relevant codon translation (shown in DNA alphabet by NCBI convention) is: citeturn25view0  

- **Reference codon (DNA, coding direction): ATG → Met (M)**  
- **Mutant codon (DNA, coding direction): GTG → Val (V)**  

NCBI’s vertebrate mitochondrial code table explicitly maps **ATG to Met** and **GTG to Val**. citeturn25view0  

So, the explicit codon translation step is:

> **ATG (Met) → GTG (Val)**, producing **Met→Val at residue 64**. citeturn25view0turn31view0turn8view0  

## HGVS protein-level consequence and amino-acid substitution

### HGVS protein notation

ClinVar-supported HGVS consequence for this mtDNA variant includes the RefSeq protein accession and protein-level HGVS:

- **NC_012920.1:m.14484T>C (YP_003024037.1:p.Met64Val)** citeturn31view0  

Therefore, the protein-level consequence is:

- **HGVS (protein): p.Met64Val** citeturn31view0  
- **Amino-acid substitution (three-letter): Met → Val** citeturn31view0  
- **Amino-acid substitution (one-letter): M64V** (also used by MITOMAP and UniProt annotations). citeturn2view0turn8view0  

### Cross-resource concordance

MITOMAP’s LHON mutation table lists **m.14484T>C** in **ND6** with amino-acid change **M64V**. citeturn2view0  
UniProt (archived reviewed entry for MT‑ND6, P03923) records a **VARIANT at position 64: M→V** described as a primary LHON mutation. citeturn8view0  

## Authoritative records, identifiers, and clinical significance notes

Authoritative curated sources used (as requested): entity["organization","MITOMAP","human mtDNA database"]; entity["organization","ClinVar","ncbi variant database"]; entity["organization","UniProt","protein knowledgebase"]; entity["organization","NCBI RefSeq","reference sequence database"]. citeturn29view0turn30view3turn8view0turn10search0  

### Direct authoritative record links (via citations)

- **MITOMAP (LHON primary mutations table)** showing **m.14484T>C; ND6; M64V**. citeturn2view0  
- **MITOMAP (confirmed pathogenic mutations list)** includes **Coding MT‑ND6 LHON m.14484 T>C … M‑V**. citeturn32view0  
- **ClinVar variant record (VCV000009688)** for **NC_012920.1(MT‑ND6):m.14484T>C** (includes curated submissions; one submission explicitly gives **YP_003024037.1:p.Met64Val**). citeturn30view3turn31view0  
- **ClinVar condition record (RCV000010325)** for **Leber optic atrophy (LHON)**: aggregate **Pathogenic**, multiple submitters, with a recent evaluation date shown on the record. citeturn1view1  
- **UniProt MT‑ND6 (P03923)** archived reviewed flatfile entry including **VARIANT 64 M→V (LHON; primary mutation)** and the canonical protein sequence context. citeturn8view0  
- **NCBI RefSeq mitochondrial genome record**: **NC_012920.1** is the human mitochondrial reference sequence used for coordinate numbering in ClinVar/MITOMAP. citeturn10search0turn29view0  

### Clinical significance annotations (high-level)

ClinVar’s LHON condition record (Leber optic atrophy) aggregates to **Pathogenic** with multiple submissions and no conflicts, with record metadata including a recent evaluation date (Oct 22, 2025). citeturn1view1  
MITOMAP lists the variant among LHON primary mutations and also includes it in its confirmed pathogenic mutation listings. citeturn2view0turn32view0  

## Isoform/annotation ambiguity and practical notes for clinical and MD interpretation

There is **no alternative splicing isoform ambiguity** in the usual nuclear sense for MT‑ND6; however, **reference annotations can differ slightly in total protein length** depending on how the terminal region/stop is represented. In the materials retrieved here, UniProt’s reviewed entry presents MT‑ND6 as **174 amino acids**, while MITOMAP’s ND6 locus summary reports **ND6 Total AA = 175**. citeturn8view0turn27view0  
Critically for **m.14484T>C**, curated sources agree on the **residue index 64** and the substitution **Met→Val**, so LHON/MD structure mapping at this site is typically stable across common reference frameworks. citeturn31view0turn2view0turn8view0