# Project explanation: Cardiometabolic GWAS Evidence Explorer

## Abstract

Genome-wide association studies can identify regions of the genome associated with type 2 diabetes, but an association signal does not by itself identify the variant or gene responsible for the underlying biology. Each region can contain several plausible variants and nearby genes, while fine-mapping, distance-based evidence, gene-ranking models, and molecular-QTL colocalisation answer different questions and may be incomplete. This makes it easy to mistake a high ranking, a large record count, or a missing result for stronger evidence than the source actually provides.

The Cardiometabolic GWAS Evidence Explorer addresses this interpretation problem by organising existing public evidence from Open Targets release 26.09 for the FinnGen R12 type 2 diabetes study `FINNGEN_R12_T2D`. From 368 source credible sets, the pipeline deterministically selected five regions using the smallest reported p-values and a one-megabase spacing rule for lead variants on the same chromosome. It verified nine cached API responses and their hashes, checked pagination and source counts, preserved missing values and source quality-control warnings, and calculated lead-variant distance to the canonical transcription start site of each returned candidate gene. The resulting snapshot contains 65 fine-mapped variant rows, 11 Open Targets Locus-to-Gene candidates, and 44 molecular-QTL colocalisation rows. Thirteen colocalisation rows matched a returned L2G candidate, while 31 were retained as evidence involving other genes at the selected regions.

The explorer presents variant posterior inclusion probability, L2G score, COLOC-PIP H4, and eCAVIAR CLPP as separate source measures. It does not combine them into a new probability, claim that a gene causes diabetes, validate the L2G model, or measure clinical usefulness. Its contribution is a transparent, inspectable workflow for comparing what evidence was returned, seeing where evidence is absent, and tracing every displayed result back to a saved public response. The five-region analysis is a bounded descriptive case study from one Finnish study and should be extended across studies, ancestries, traits, and independent validation data before drawing general conclusions.

## The problem it addresses

A GWAS result usually identifies a **region**, not a confirmed causal gene. Moving from that regional signal to a gene that merits further investigation is difficult for four practical reasons:

1. **Several variants can remain plausible.** Fine-mapping distributes support across variants in a region; it does not always isolate one variant.
2. **Several genes can be plausible.** The closest gene is not automatically the relevant gene, and an L2G score is a model ranking rather than experimental proof.
3. **Molecular evidence is uneven.** Some regions have many QTL and colocalisation records, while others have none returned. Missing evidence must not be interpreted as a score of zero.
4. **The measures are easy to misuse.** PIP, L2G, H4, and CLPP come from different models. Adding or averaging them would create an unsupported score and could double-count related evidence.

The project therefore asks a narrower and defensible question:

> For five selected FinnGen type 2 diabetes regions, what variant, candidate-gene, and molecular-QTL evidence did Open Targets return, and where are the evidence gaps?

## How it works

### 1. Freeze the source snapshot

The collection pipeline queries the official Open Targets GraphQL API and saves the exact request, variables, response body, retrieval time, byte count, and SHA-256 hash. The committed snapshot uses Open Targets data release 26.09, genome build GRCh38, and study `FINNGEN_R12_T2D`.

This matters because Open Targets releases and L2G predictions can change. A saved snapshot makes the result traceable to one specific source state.

### 2. Select a bounded set of regions

The source returned 368 credible sets for the study. The pipeline selected five deterministically by:

- sorting by the smallest reported p-values;
- requiring lead variants on the same chromosome to be at least one megabase apart; and
- using stable identifiers to break ties.

This is a sampling rule for a manageable case study. It is not proof that the five regions are statistically independent or representative of all type 2 diabetes biology.

### 3. Retrieve three evidence layers

For each selected region, the pipeline collects:

- **Fine-mapped variants:** the source credible-set variants and their posterior inclusion probabilities;
- **Candidate genes:** the genes returned by the Open Targets L2G model, including source features and SHAP contributions; and
- **Molecular-QTL colocalisation:** returned gene, tissue, study, method, H3, H4, and CLPP fields.

The included snapshot contains 65 variant rows, 11 candidate-gene rows, and 44 molecular-QTL colocalisation rows.

### 4. Validate before presenting

The pipeline rejects partial or inconsistent evidence instead of silently continuing. It checks:

- cached-response hashes;
- GraphQL errors;
- duplicate identifiers and pages;
- reported versus fetched row counts;
- pagination completeness;
- probability ranges and PIP sums; and
- coordinate compatibility before calculating distances.

Every selected region carries the exact source warning `Study has quality control flag(s)`. The saved response does not explain the warning, so the project displays it without guessing its effect.

### 5. Match molecular evidence to returned candidates

A colocalisation row is matched to an L2G candidate only when the source QTL gene identifier or source target identifier equals the candidate's Ensembl gene identifier. Rows involving other genes remain visible as locus-wide evidence rather than being reassigned to a candidate.

For example, the CCND2 region contains 38 returned molecular-QTL colocalisation rows, but only 11 match genes in the returned L2G candidate list. The remaining 27 rows are retained and labelled as evidence for other genes at the region. This prevents the total row count from being misrepresented as 38 confirmations of CCND2.

### 6. Keep unlike measures separate

The interface places the evidence layers next to one another but does not merge them:

| Measure | What it describes | What it does not establish |
|---|---|---|
| Variant PIP | Source fine-mapping support assigned to a variant within a region | That a gene causes type 2 diabetes |
| L2G score | Open Targets model ranking for a returned gene at a GWAS region | Experimental or clinical validation |
| H4 | COLOC-PIP support for a shared-variant model between two signals | Universal proof of one causal mechanism |
| CLPP | eCAVIAR colocalisation statistic | A value interchangeable with H4 |
| Not returned | No corresponding source record was present in the saved response | A biological value of zero |

## Why this is important

Turning a GWAS association into a biological hypothesis is a key step before functional follow-up or therapeutic target research. Poorly handled evidence can create false confidence: the nearest gene may be assumed to be causal, repeated tissues may be counted as independent confirmation, missing records may be converted to zeros, or correlated measures may be combined twice.

This project makes those decision points visible. A researcher or student can inspect the variant shortlist, see the genes ranked by the existing L2G model, examine the tissues and genes represented in molecular-QTL records, and determine whether a displayed molecular row actually matches a returned candidate. The saved responses, manifest, downloadable tables, and deterministic checks also make each claim auditable.

The project is useful as:

- a teaching tool for understanding the path from a GWAS region to candidate genes;
- an evidence-audit layer before designing a more advanced statistical model;
- a transparent starting point for discussing missing and dependent evidence; and
- a data foundation for future cross-study or Bayesian target-prioritisation work.

## What the project produces

- A self-contained [interactive explorer](https://zinga18018.github.io/cardiometabolic-gwas-evidence-explorer/).
- Three analysis-ready CSV files for variants, L2G candidates, and molecular-QTL colocalisations.
- A raw ZIP containing the manifest and nine byte-identical API response files.
- A five-region [real-data case study](reports/real_data_case_study.md).
- Deterministic builders and tests for the explorer, case study, and exported data.

## Limitations

- The five regions are a deterministic convenience sample from one Finnish study.
- The analysis is descriptive and does not estimate predictive performance.
- Source L2G scores already use correlated features, including distance and colocalisation-related information.
- Maximum H4 or CLPP across tissues is descriptive and is not independent replication.
- The source quality-control warning for every selected region remains unresolved.
- The project does not run a new GWAS, fine-mapping model, Mendelian randomisation analysis, causal-gene model, wet-lab experiment, or clinical validation.
- All derived records remain `public_source_derived_pending_human_review`.

## Primary source documentation

- [Open Targets credible sets and variant PIP](https://platform-docs.opentargets.org/credible-set)
- [Open Targets Locus-to-Gene methodology](https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g)
- [Open Targets colocalisation methods](https://platform-docs.opentargets.org/gentropy/colocalisation)
- [Open Targets GraphQL API](https://platform-docs.opentargets.org/data-access/graphql-api)

## A concise interview explanation

> I built an evidence explorer for five FinnGen type 2 diabetes regions using real public aggregate records from Open Targets. The project starts from a GWAS region and shows three different layers of existing evidence: fine-mapped variants, Open Targets L2G candidate genes, and molecular-QTL colocalisation records. I validated the saved API responses, preserved missing values and source quality-control warnings, and matched molecular rows to candidate genes using Ensembl identifiers. I kept PIP, L2G, H4, and CLPP separate because they answer different questions. The explorer helps a reader understand why a region may point toward a gene, but it does not claim that the gene causes diabetes.
