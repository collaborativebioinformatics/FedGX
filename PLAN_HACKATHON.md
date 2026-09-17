____________________________

## Initial Plan 

1. Synthetic dataset creation across 3 sites (genotype and phenotype)

- Use LDAK from DougSpeed.com
- Generate genotype across three different populations (ancestries)
- Genotype/phenotype mapping

2. Client pipeline for GWAS software in *Site 1, 2, 3*

- REGENIE
- Instructions for extension to PLINK, GCTA, SAIGE
- Edits to server-side GWAS code to handle different GWAS calls
- Possibly LD structure handling

3. Standardization of summary stats in *Central Analytical Engine*

4. Meta-analyses in *Central Analytical Engine* using GWAMA

- FFX — currently implemented
- RFX — currently buggy
- Optionally develop LD structure weighted meta-analyses

5. Visualisation component (bonus, not core scope)

6. PRS (open question — not yet scoped)
