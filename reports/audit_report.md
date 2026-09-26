# Cardiometabolic GWAS evidence audit

Snapshot: 2026-09-25T14:04:47.212407Z. Open Targets data release 26.09.

Study: [FINNGEN_R12_T2D](https://platform.opentargets.org/study/FINNGEN_R12_T2D) — Type 2 diabetes, definitions combined.
Selected 5 of 368 credible sets by the smallest reported p-values, requiring lead variants on the same chromosome to be at least 1 Mb apart. This is a convenience sample, not an independent-locus proof or performance benchmark.

| Lead variant | Top returned L2G gene | Variants fetched | PIP sum | L2G candidates | Molecular-QTL colocalisation rows | Source quality controls |
|---|---|---:|---:|---:|---:|---|
| [10_112998590_C_T](https://platform.opentargets.org/credible-set/83b58cde6c7d78f3211dce1b35b83674) | TCF7L2 | 1 | 0.999981 | 1 | 0 | Study has quality control flag(s) |
| [12_4275678_T_G](https://platform.opentargets.org/credible-set/25c86693e33ae778acf27c3255c2925f) | CCND2 | 1 | 1.000000 | 2 | 38 | Study has quality control flag(s) |
| [16_53793567_A_G](https://platform.opentargets.org/credible-set/911b5c2383e73877f5d1a394a5ed728c) | FTO | 56 | 0.950564 | 3 | 2 | Study has quality control flag(s) |
| [20_44354839_C_CA](https://platform.opentargets.org/credible-set/f26ca97da3415b5a097fdf85b538997e) | HNF4A | 1 | 1.000000 | 4 | 0 | Study has quality control flag(s) |
| [6_20703721_A_G](https://platform.opentargets.org/credible-set/8b114fb0fc75e8f34eece8cee51dd41f) | CDKAL1 | 6 | 0.999189 | 1 | 4 | Study has quality control flag(s) |

## What was computed

Parsed and validated source probabilities, checked pagination and cached-response hashes, preserved missing values, and calculated absolute distances from lead variants to canonical transcription start sites for returned candidates. Negative-strand transcripts use their end coordinate. L2G and colocalisation values were supplied by Open Targets; they were not estimated by this project.

## Interpretation boundaries

- A missing colocalisation row means no row was returned by this query and release; it is not evidence against a gene.
- L2G scores are existing model outputs; they are not this project's causal posterior or clinical target probability. Colocalisation H4 and eCAVIAR CLPP are distinct fields, not interchangeable scores.
- The nearest comparison covers only returned L2G candidates, which are filtered upstream; it does not identify the nearest gene across the genome. Agreement is descriptive, not accuracy.
- Maximum H4 is the maximum among fetched molecular-QTL records. It may reflect a tissue unrelated to the disease mechanism. Multiple QTL signals and tissues are not independent replications.
- Source quality-control text is reproduced exactly from the `qualityControls` field. The cached Open Targets API response provides no additional detail about a returned flag, so this project does not infer its meaning or effect.
- No new GWAS, fine-mapping, Mendelian randomisation, drug-target validation, wet-lab experiment or AHBA analysis was performed.
- Study-level nSamples and discoverySamples are retained as separate source annotations and can differ. No denominator was silently changed.
- Every derived record remains public_source_derived_pending_human_review.

## Reproduce

Run `python gwas_explorer.py verify`, `python gwas_explorer.py build`, and `python -m unittest discover -s tests -v` from the repository root. A new download requires a new snapshot directory.
