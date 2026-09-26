# Real-data case study: evidence heterogeneity across FinnGen T2D loci

**Status:** `public_source_derived_pending_human_review`

## Research question

How do fine-mapping concentration, Locus-to-Gene rankings, molecular-QTL colocalisation, and evidence availability differ across five selected FinnGen R12 type 2 diabetes loci?

## Data and scope

This descriptive case study uses the saved Open Targets 26.09 snapshot for `FINNGEN_R12_T2D` (Type 2 diabetes, definitions combined). The source contains 368 credible sets; 5 were selected by the documented p-value and 1 Mb spacing rule. It includes 65 variant rows, 11 returned L2G candidates, and 44 molecular-QTL colocalisation rows.

These are public summary-level source records, not participant-level data. PIP, L2G, H4, and CLPP are reported separately. No composite score, causal posterior, core-gene posterior, therapeutic ranking, or model-performance estimate is calculated.

## Locus-level results

| Lead variant | Top returned L2G gene | Variant rows | PIP sum | Maximum returned PIP | L2G candidates | Total QTL colocs | Candidate-matched colocs | Locus-wide max H4 | Returned-candidate max H4 | Locus-wide max CLPP | Returned-candidate max CLPP |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| [10_112998590_C_T](https://platform.opentargets.org/credible-set/83b58cde6c7d78f3211dce1b35b83674) | TCF7L2 | 1 | 0.999981 | 0.999981 | 1 | 0 | 0 | not returned | not returned | not returned | not returned |
| [12_4275678_T_G](https://platform.opentargets.org/credible-set/25c86693e33ae778acf27c3255c2925f) | CCND2 | 1 | 1.000000 | 1.000000 | 2 | 38 | 11 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| [16_53793567_A_G](https://platform.opentargets.org/credible-set/911b5c2383e73877f5d1a394a5ed728c) | FTO | 56 | 0.950564 | 0.045960 | 3 | 2 | 0 | 0.933400 | not returned | 0.012512 | not returned |
| [20_44354839_C_CA](https://platform.opentargets.org/credible-set/f26ca97da3415b5a097fdf85b538997e) | HNF4A | 1 | 1.000000 | 1.000000 | 4 | 0 | 0 | not returned | not returned | not returned | not returned |
| [6_20703721_A_G](https://platform.opentargets.org/credible-set/8b114fb0fc75e8f34eece8cee51dd41f) | CDKAL1 | 6 | 0.999189 | 0.339572 | 1 | 4 | 2 | 0.992549 | 0.985572 | 0.113979 | 0.061083 |

## What the five loci show

- **TCF7L2 / 10_112998590_C_T:** 1 returned variant row(s) carried PIP mass 0.999981; Open Targets returned 1 L2G candidate(s), with top score 0.8700; 0 molecular-QTL row(s) were returned and 0 matched a returned candidate by gene ID.
  A single returned candidate provides no within-locus ranking comparison.
  No colocalisation row was returned; this is missing source evidence, not a zero or negative biological result.
- **CCND2 / 12_4275678_T_G:** 1 returned variant row(s) carried PIP mass 1.000000; Open Targets returned 2 L2G candidate(s), with top score 0.3482; 38 molecular-QTL row(s) were returned and 11 matched a returned candidate by gene ID.
  27 colocalisation row(s) targeted genes outside the returned L2G candidate list; they were retained as locus evidence and were not reassigned.
- **FTO / 16_53793567_A_G:** 56 returned variant row(s) carried PIP mass 0.950564; Open Targets returned 3 L2G candidate(s), with top score 0.9876; 2 molecular-QTL row(s) were returned and 0 matched a returned candidate by gene ID.
  2 colocalisation row(s) targeted genes outside the returned L2G candidate list; they were retained as locus evidence and were not reassigned.
- **HNF4A / 20_44354839_C_CA:** 1 returned variant row(s) carried PIP mass 1.000000; Open Targets returned 4 L2G candidate(s), with top score 0.5327; 0 molecular-QTL row(s) were returned and 0 matched a returned candidate by gene ID.
  No colocalisation row was returned; this is missing source evidence, not a zero or negative biological result.
- **CDKAL1 / 6_20703721_A_G:** 6 returned variant row(s) carried PIP mass 0.999189; Open Targets returned 1 L2G candidate(s), with top score 0.9386; 4 molecular-QTL row(s) were returned and 2 matched a returned candidate by gene ID.
  A single returned candidate provides no within-locus ranking comparison.
  2 colocalisation row(s) targeted genes outside the returned L2G candidate list; they were retained as locus evidence and were not reassigned.

## Transformations performed

1. Read only the committed `reports/evidence.json` snapshot and validated its schema, identifiers, counts, completeness fields, probability ranges, and PIP sums.
2. Preserved source PIPs, L2G scores, H3, H4, and CLPP without normalising or combining them.
3. Ranked the returned L2G candidates by their existing source score, with target ID only as a deterministic tie-breaker.
4. Matched a colocalisation row to a candidate when either `otherStudyLocus.qtlGeneId` or `otherStudyLocus.study.target.id` equalled the candidate Ensembl ID, matching the explorer's existing rule.
5. Retained QTL type, tissue, study, method, quality-control text, completeness fields, and missing values in the long table.
6. Drew an availability chart whose counts indicate returned records, not evidence strength.

## Interpretation boundaries

- Variant PIP, L2G score, coloc H4, and eCAVIAR CLPP answer different questions and are not interchangeable probabilities.
- L2G already uses features that can include distance and colocalisation. Combining the displayed values again would risk double-counting dependent evidence.
- A maximum across tissues or studies is descriptive and subject to multiplicity; it is not independent replication.
- A missing colocalisation row remains `not_returned`. It is not converted to zero and is not evidence against a gene.
- Returned L2G candidates are an upstream-filtered set. Their scores are not renormalised, and a one-candidate locus does not validate ranking accuracy.
- Every selected locus carries the exact unresolved source warning `Study has quality control flag(s)`. The snapshot does not explain its meaning, so this analysis does not speculate about its effect.
- The five loci are a deterministic convenience sample from one Finnish study. They do not estimate performance across loci, cohorts, ancestries, or diseases.
- The analysis does not run a new GWAS, fine mapping, Mendelian randomisation, causal-gene model, CoreGene posterior, or clinical target validation.
- All derived outputs remain public-source-derived and pending human review.

## Reproduce

```bash
python analysis/build_case_study.py
python analysis/build_case_study.py --check
python -m unittest discover -s tests -v
```
