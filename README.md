# Cardiometabolic GWAS Evidence Explorer

[![tests](https://github.com/Zinga18018/cardiometabolic-gwas-evidence-explorer/actions/workflows/tests.yml/badge.svg)](https://github.com/Zinga18018/cardiometabolic-gwas-evidence-explorer/actions/workflows/tests.yml)

A reproducible research prototype for inspecting **existing** type 2 diabetes GWAS credible sets, Locus-to-Gene (L2G) prioritisation and molecular-QTL colocalisation from the Open Targets Platform.

New to genetics? Start with the [beginner and interview guide](BEGINNER_GUIDE.md). It explains genome, DNA, genes, variants, GWAS, fine-mapping, PIP, QTLs, colocalisation, L2G, H3, H4 and CLPP in plain language.

**Development snapshot:** September 2026<br>
**Initial public release:** 25 September 2026

## Abstract

Genome-wide association studies can identify regions associated with type 2 diabetes, but a regional signal does not by itself identify the responsible variant or gene. This project organises existing public evidence from Open Targets release 26.09 for five selected regions from the FinnGen R12 type 2 diabetes study. It brings together 65 fine-mapped variant rows, 11 Locus-to-Gene candidates, and 44 molecular-QTL colocalisation rows while preserving source identifiers, missing values, pagination checks, response hashes, and quality-control warnings. Variant PIP, L2G score, COLOC-PIP H4, and eCAVIAR CLPP remain separate because they describe different models and cannot be treated as one causal probability. The result is an auditable explorer and downloadable dataset for examining how variant, gene-ranking, and molecular evidence differ across regions. It is a bounded descriptive case study, not causal-gene validation, a clinical target ranking, or an analysis of participant-level data.

Read the detailed [problem, workflow, importance, and limitations](PROJECT_EXPLANATION.md).

## In simple terms

A GWAS can point to a region of the genome associated with a disease, but it usually does not identify the responsible gene. This project examines existing Open Targets records for five selected FinnGen type 2 diabetes regions and places three kinds of clues side by side:

- **Fine-mapping:** which variants remain plausible within the associated region.
- **Locus-to-Gene (L2G):** which genes Open Targets ranks for that region.
- **Molecular-QTL colocalisation:** whether the disease signal and a molecular signal may share an underlying variant in a particular biological context.

The project checks that the saved records are complete and unchanged, matches returned molecular evidence to candidate genes, and keeps missing evidence visible. It does not select a causal gene or turn these different measures into a new combined probability.

## Is this regression or classification?

**This repository is neither a regression model nor a classification model.** It is a deterministic post-GWAS evidence-integration and audit pipeline. It trains no model and predicts no patient outcome.

The upstream FinnGen GWAS uses logistic mixed-model association testing because type 2 diabetes is a binary case/control phenotype. The upstream Open Targets L2G model is trained as a gradient-boosting classifier and its score is used to rank candidate genes. Fine-mapping and colocalisation are probabilistic inference stages. This repository retrieves, validates, matches and visualises those existing outputs. See the [beginner guide](BEGINNER_GUIDE.md#is-this-regression-or-classification) for the full distinction and an interview answer.

Open the [live interactive explorer](https://zinga18018.github.io/cardiometabolic-gwas-evidence-explorer/), use [the offline version](reports/index.html), or read the [real-data case study](reports/real_data_case_study.md). The explorer includes locus selection, gene-score filtering, tissue search and colocalisation filtering. No account, build system, network connection or external JavaScript library is needed to view the saved report.

The exact source-derived tables are available as [65 variant rows](reports/dataset/variants.csv), [11 L2G candidate rows](reports/dataset/l2g_candidates.csv), and [44 molecular-QTL colocalisation rows](reports/dataset/molecular_qtl_colocalisations.csv). The [dataset guide](reports/dataset/README.md) defines every field. The [raw Open Targets snapshot](reports/dataset/open_targets_raw_snapshot.zip) contains the manifest and nine original API response files for this bounded five-region extract.

## Research question

What source evidence was returned for candidate genes at a small set of cardiometabolic association signals, and which apparent gaps are missing data rather than negative biological evidence?

The committed real-data case study asks a more focused question: **How do fine-mapping concentration, L2G rankings, molecular-QTL colocalisation, and evidence availability differ across five selected FinnGen R12 type 2 diabetes loci?**

This project audits public outputs. It does **not** run a new GWAS, fine mapping, Mendelian randomisation or a clinical drug-target validation. PIPs, L2G scores and colocalisation statistics originate from Open Targets. Only the retrieval, validation, derived distances and presentation are implemented here.

## Reproduce the included snapshot

Requirements: Python 3.11 or later. Standard library only; no packages to install.

```bash
python gwas_explorer.py verify
python gwas_explorer.py build
python analysis/build_case_study.py
python analysis/build_case_study.py --check
python analysis/export_real_dataset.py
python analysis/export_real_dataset.py --check
python -m unittest discover -s tests -v
```

Open `reports/index.html` locally. To serve the repository for browser preview:

```bash
python -m http.server 8765
```

Then visit `http://localhost:8765/reports/`.

## Real-data case study

The saved `FINNGEN_R12_T2D` snapshot contains real public aggregate credible-set, L2G, and colocalisation records from Open Targets, not participant-level data or a complete GWAS summary-statistics file. The deterministic case-study build compares five selected loci using **65 variant rows, 11 returned L2G candidates, and 44 molecular-QTL colocalisation rows**.

The analysis preserves PIP, L2G, H4, and CLPP as separate source measures. It matches colocalisation records to candidates using the source Ensembl gene identifiers, retains unmatched records and missing values, and carries every source quality-control warning into the outputs. It does not create a combined score or posterior.

Read the [case-study report](reports/real_data_case_study.md), inspect the [locus summary](reports/locus_summary.csv) and [candidate-level long table](reports/candidate_evidence_long.csv), or view the [evidence-availability figure](reports/evidence_availability.svg). The counts in the figure describe returned records rather than evidence strength. All derived records remain `public_source_derived_pending_human_review`.

## Fetch a new snapshot

```bash
python gwas_explorer.py fetch --data-dir data-refresh-20261001 --report-dir reports-refresh-20261001 --loci 5
```

Choose a **new** snapshot directory for each collection. Existing snapshots cannot be overwritten by the fetch command. The fetch is bounded to 3–8 loci and at most eight pages per endpoint (500 rows per page). A partial evidence fetch is labelled partial; incomplete study discovery fails rather than silently selecting from a truncated list. API release metadata is checked before and after collection.

## Method

1. Query the current official GraphQL schema and source release metadata. Save the exact query, variables, timestamps, raw response bytes and SHA-256 hash for every response.
2. Retrieve all credible-set metadata for `FINNGEN_R12_T2D`. This study was selected as one coherent public type 2 diabetes study, rather than pooling overlapping cohorts.
3. Select the smallest reported p-values deterministically, with at least 1 Mb between lead variants on the same chromosome. This separation is a sampling convention, not an LD-independent-locus guarantee. No biological success criterion is used for selection.
4. Retrieve variant PIPs, source L2G candidates and all returned molecular-QTL colocalisation pages for each selected locus. Reject GraphQL errors, duplicated identifiers, changed pagination counts and invalid probabilities.
5. Derive absolute lead-variant distance to each returned candidate's canonical TSS using GRCh38 coordinates. Use transcript start for positive strand and end for negative strand. Unknown coordinates remain missing.
6. Preserve source scores and null values. Group colocalisation records by the QTL gene identifier, retaining tissue, study and method in the underlying JSON.
7. Generate a compact evidence CSV, complete structured JSON, audit report and self-contained interactive HTML.

The “nearest returned candidate” comparison is limited to genes present in the filtered L2G response. It is **not a genome-wide nearest-gene baseline**. A one-candidate locus trivially agrees and provides no comparative validation. Agreement is descriptive, never accuracy.

## Evidence boundaries

- L2G is an existing machine-learning score, not a new causal posterior computed here. Its constituent features can already include distance and colocalisation; displaying them together does not create independent evidence.
- H4 and eCAVIAR CLPP are kept separate. No universal colocalisation threshold is asserted; the UI slider is exploratory.
- No returned QTL row means `not_returned`, not zero probability and not a negative finding.
- Maximum H4 is descriptive over fetched rows and can come from a biologically irrelevant tissue. Multiple tissues or overlapping QTL studies are not independent replications.
- Every selected locus in the included snapshot carries the exact source `qualityControls` text `Study has quality control flag(s)`. The cached API response provides no additional detail, so this project does not infer the flag's meaning or effect.
- No pLOF burden, Perturb-seq, gene network or therapeutic-actionability model is implemented.
- Public-derived records remain `public_source_derived_pending_human_review`.
- The sample contains loci selected from the smallest reported p-values in one Finnish study. It cannot establish model performance or generalise across ancestries.
- Source study counts and discovery-sample annotations are preserved separately, even when their denominators differ.

## Layout

| Path | Purpose |
|---|---|
| `PROJECT_EXPLANATION.md` | Detailed abstract, problem statement, workflow, importance, example, limitations and interview explanation |
| `BEGINNER_GUIDE.md` | Plain-language genetics glossary, method classification, pipeline walkthrough and defense-ready interview answers |
| `gwas_explorer.py` | Retrieval, integrity checks, validation, evidence compilation and CLI |
| `explorer_template.html` | Offline interactive report template |
| `data/manifest.json` | Exact requests, response hashes, UTC times, selection rule |
| `data/raw/` | Unmodified official API responses, including schema and release metadata |
| `reports/evidence.json` | Complete fetched evidence and derived distances |
| `reports/candidate_evidence.csv` | One row per locus–gene candidate; blanks retain missing values |
| `reports/audit_report.md` | Actual run counts and interpretation boundaries |
| `reports/run_summary.json` | Machine-readable run totals |
| `analysis/build_case_study.py` | Deterministic real-data case-study builder and validation |
| `reports/real_data_case_study.md` | Five-locus descriptive analysis with explicit claim boundaries |
| `reports/locus_summary.csv` | One row per selected locus; source measures remain separate |
| `reports/candidate_evidence_long.csv` | Candidate–colocalisation matches, candidate missingness, and unmatched locus colocs with tissue and method retained |
| `reports/evidence_availability.svg` | Accessible availability chart; counts are not evidence strength |
| `reports/case_study_metadata.json` | Input hash, source release, output hashes and analysis constraints |
| `analysis/export_real_dataset.py` | Deterministic exporter for analysis-ready CSVs and the raw-source archive |
| `reports/dataset/` | Three source-derived CSVs, data dictionary and raw five-region API snapshot ZIP |
| `tests/` | Offline validation tests |
| `.github/workflows/tests.yml` | Verification and unit-test workflow for each push and pull request |
| `.github/workflows/pages.yml` | Deployment of the saved explorer to GitHub Pages |
| `LICENSE` | MIT licence for the original project code |

## Fit with Bayesian target-prioritisation research

This is a data and provenance foundation for studying the variant-to-gene portion of a larger framework: source provenance, uncertain gene mappings, tissue-specific evidence and honest missingness. The next research step would be a clearly specified probabilistic model with independent validation, assumptions about dependent evidence, and evaluation across studies. The current project does not implement that full framework.

## Sources and data licence

- [Open Targets official GraphQL documentation](https://platform-docs.opentargets.org/data-access/graphql-api)
- [Credible-set documentation](https://platform-docs.opentargets.org/credible-set)
- [Locus-to-Gene methodology](https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g)
- [Variant annotation and GRCh38](https://platform-docs.opentargets.org/variant)
- [Open Targets source and licence information](https://platform-docs.opentargets.org/licence)
- [FinnGen T2D study in Open Targets](https://platform.opentargets.org/study/FINNGEN_R12_T2D)
- Buniello, A. et al. (2025). [Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery](https://academic.oup.com/nar/article/53/D1/D1467/7917960). *Nucleic Acids Research*.

Open Targets marks its platform data CC0 and advises checking underlying-source conditions. The source data remain attributed to Open Targets and the originating studies. This repository makes no claim of affiliation with either organisation. No participant-level data are downloaded.

The original project code in this repository is released under the MIT License. Public data and source-derived records remain subject to their original source terms.
