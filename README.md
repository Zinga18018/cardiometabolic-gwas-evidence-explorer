# Cardiometabolic GWAS Evidence Explorer

[![tests](https://github.com/Zinga18018/cardiometabolic-gwas-evidence-explorer/actions/workflows/tests.yml/badge.svg)](https://github.com/Zinga18018/cardiometabolic-gwas-evidence-explorer/actions/workflows/tests.yml)

A reproducible research prototype for inspecting **existing** type 2 diabetes GWAS credible sets, Locus-to-Gene (L2G) prioritisation and molecular-QTL colocalisation from the Open Targets Platform.

**Development snapshot:** September 2026<br>
**Initial public release:** 25 September 2026

Open the [live interactive explorer](https://zinga18018.github.io/cardiometabolic-gwas-evidence-explorer/), use [the offline version](reports/index.html), or read [the audit report](reports/audit_report.md). The explorer includes locus selection, gene-score filtering, tissue search and colocalisation filtering. No account, build system, network connection or external JavaScript library is needed to view the saved report.

## Research question

What evidence supports the candidate genes at a small set of cardiometabolic association signals, and which apparent gaps are missing data rather than negative biological evidence?

This project audits public outputs. It does **not** run a new GWAS, fine mapping, Mendelian randomisation or a clinical drug-target validation. PIPs, L2G scores and colocalisation statistics originate from Open Targets. Only the retrieval, validation, derived distances and presentation are implemented here.

## Reproduce the included snapshot

Requirements: Python 3.11 or later. Standard library only; no packages to install.

```bash
python gwas_explorer.py verify
python gwas_explorer.py build
python -m unittest discover -s tests -v
```

Open `reports/index.html` locally. To serve the repository for browser preview:

```bash
python -m http.server 8765
```

Then visit `http://localhost:8765/reports/`.

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
- The sample contains strong signals from one Finnish study. It cannot establish model performance or generalise across ancestries.
- Source study counts and discovery-sample annotations are preserved separately, even when their denominators differ.

## Layout

| Path | Purpose |
|---|---|
| `gwas_explorer.py` | Retrieval, integrity checks, validation, evidence compilation and CLI |
| `explorer_template.html` | Offline interactive report template |
| `data/manifest.json` | Exact requests, response hashes, UTC times, selection rule |
| `data/raw/` | Unmodified official API responses, including schema and release metadata |
| `reports/evidence.json` | Complete fetched evidence and derived distances |
| `reports/candidate_evidence.csv` | One row per locus–gene candidate; blanks retain missing values |
| `reports/audit_report.md` | Actual run counts and interpretation boundaries |
| `reports/run_summary.json` | Machine-readable run totals |
| `tests/` | Offline validation tests |
| `.github/workflows/tests.yml` | Verification and unit-test workflow for each push and pull request |
| `.github/workflows/pages.yml` | Deployment of the saved explorer to GitHub Pages |
| `LICENSE` | MIT licence for the original project code |

## Fit with Bayesian target-prioritisation research

This is a practical foundation for the variant-to-gene portion of a larger framework: source provenance, uncertain gene mappings, tissue-specific evidence and honest missingness. The next research step would be a clearly specified probabilistic model with independent validation, assumptions about dependent evidence, and evaluation across studies. The current project does not implement that full framework.

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
