# Figure 1 Figdraw production prompt

This is a reconstructed final production specification matching the v7 artwork. It is retained for reproducibility and is not represented as a byte-for-byte execution log. The authors supplied the finished artwork as a raster PNG; no editable vector source is claimed.

## Prompt

Create from scratch a publication-quality editable vector schematic for a retinal-layer genetics manuscript. Use a white landscape canvas, restrained Nature-style design, Arial or Helvetica, thin gray arrows, pale fills, and no title, caption, gradients, shadows, cartoons, DNA icons, or decorative elements. Use three equal vertical columns with aligned INPUTS, ANALYSES, and INTERPRETATION cards. Keep all text inside boxes and use no crossing arrows.

Place a top evidence key: T1 Primary statistical evidence; T2 Sensitivity/exploratory; T3 Descriptive context; U External validation required.

COLUMN 1 title: “1 Layer and disease genetics”.
INPUTS: “OCT layer GWAS”; “9 traits; N=38,346–38,485”. “Disease GWAS”; “6 traits; 54,082 cases, 1,306,001 controls”. Do not assign these combined totals to early AMD alone.
ANALYSES: “Genetic architecture”; “LDSC; bidirectional MR”; “54 forward, 9 reverse tests”. “Multivariable and sensitivity analyses”; “Three-exposure MVMR; robustness checks”.
INTERPRETATION: “Directional retinal-layer associations”. T1: forward MR links INL and OPL thickness to early AMD; reverse MR links early AMD liability to PS thickness; no significant genome-wide early AMD LDSC correlation, with some imprecise estimates. T2: GCL–dry AMD is estimator-sensitive; INL and OPL attenuate in MVMR; RPE is rs10737680-dependent and unstable across robust estimators.

COLUMN 2 title: “2 Locus, gene and protein prioritization”.
INPUTS: “OCT + early AMD GWAS”; “eQTLGen whole-blood cis-eQTL”; “31,684 individuals”; “GCL + dry AMD GWAS”; “INTERVAL plasma cis-pQTL”.
ANALYSES: “Colocalization and SMR”; “Single-causal-variant coloc; SMR + HEIDI”. “Exploratory protein analysis”; “GCL–TNXB plasma protein pathway”.
INTERPRETATION: “Shared loci and candidate genes”. T1: “13 tests with PP.H4>0.9”; “29 pairs and 19 genes after global Bonferroni and HEIDI filtering”; “ENSG00000258311: strongest cross-dataset candidate; uncharacterized protein”. “Exploratory protein hypothesis”. T2: “TNXB indirect effect P=0.084807; mediation not confirmed”.

COLUMN 3 title: “3 Cellular and tissue context”.
INPUTS: “11 trait GWAS”; “Six-donor retinal atlas”; “16-gene candidate panel”; “Visium: 5 sections, 4,718 spots, 2 donors”.
ANALYSES: “MAGMA and scDRS”; “66 primary tests”. “Descriptive spatial analysis”; “Raw and partial matched-random nulls”.
INTERPRETATION: T3: nine global FDR-positive scDRS results, mainly non-neuronal glia; PS and RPE rod enrichment; 0/30 matched-random empirical-FDR results and no gene-set-specific spatial enrichment; ENSG00000258311 not detected under the available alias; disease status fully confounded with donor.

Add a bottom dashed U card: “Required external validation”: non-overlapping and cross-ancestry replication; retinal regulatory and protein validation; independent atlases and multi-donor spatial studies; functional perturbation. Do not add results, strengthen causal claims, or depict INL and OPL as independent causal pathways.
