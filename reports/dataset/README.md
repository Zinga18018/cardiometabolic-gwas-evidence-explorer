# Exported real-data tables

## What “real data” means here

These CSVs contain **real public aggregate records** saved from the Open Targets Platform. The disease-association side comes from the FinnGen R12 type 2 diabetes study `FINNGEN_R12_T2D`; the molecular-QTL side contains public QTL studies integrated by Open Targets.

They are **not participant-level records**, patient data, or genotypes. They are **not a complete GWAS summary-statistics download**. A row represents a returned variant, gene candidate, or molecular-QTL colocalisation record. For example, 65 variant rows means 65 variant records, not 65 people.

The export covers five deliberately selected loci from one study:

- `variants.csv`: 65 fine-mapped variant rows.
- `l2g_candidates.csv`: 11 Open Targets Locus-to-Gene candidate rows.
- `molecular_qtl_colocalisations.csv`: 44 molecular-QTL colocalisation rows.
- `open_targets_raw_snapshot.zip`: the manifest and nine original API response JSON files behind the five-locus snapshot.

The three CSVs come deterministically from `reports/evidence.json` (SHA-256 `49e73e407ba8a062c7b0a32f145233d24505dd381e675c550e5ba76267f5c821`). No new score or probability is calculated here. Missing source values are empty CSV cells; a source value of zero remains `0` or `0.0`.

## Raw source download

`open_targets_raw_snapshot.zip` is the **raw five-locus Open Targets API snapshot used by this project**. It contains `data/manifest.json` and the nine unmodified JSON response bodies under `data/raw/`, using the same relative paths as the repository. The manifest records each request, retrieval time, byte count, and response SHA-256.

This archive is **not full FinnGen GWAS summary statistics** and contains **no participant-level data**. The discovery response lists the public study's returned credible-set metadata, while five locus responses contain the rows fetched for the selected loci. The other responses record source release metadata and the API schema. “Raw” means unmodified API response bytes within this bounded saved snapshot; it does not mean individual-level or complete FinnGen source data.

The archive uses stored, uncompressed entries with fixed timestamps and file permissions so rerunning the exporter produces the same ZIP bytes on Windows or Linux. Each archived JSON file is checked against the byte count and SHA-256 recorded in the manifest before packaging.

## Key terms in plain language

- **GWAS:** a study that looks for genetic variants associated with a trait across many people.
- **Locus:** a region of the genome around an association signal.
- **Variant:** a DNA difference at a genomic position.
- **Fine-mapping / credible set:** a method and returned group of variants used to narrow an association signal. It does not prove which variant is causal.
- **PIP (`posterior_probability`):** the source fine-mapping probability attached to a variant within its locus.
- **L2G:** Open Targets' existing model for ranking genes near a GWAS signal. `l2g_score` is the source model score, not a score created by this project.
- **QTL:** a genetic association with a measured molecular trait, such as gene expression or protein abundance.
- **Colocalisation:** a statistical check of whether a disease association and a molecular-QTL association may share a causal variant.
- **H3, H4 and CLPP:** separate statistics returned by the source colocalisation analysis. H4 represents the shared-variant hypothesis in the coloc model; CLPP is the eCAVIAR statistic. They are kept separate and are not combined.
- **SHAP value:** the source model's feature contribution to its L2G prediction. Its sign describes the model contribution, not whether a gene is biologically “good” or “bad.”

## Provenance and interpretation rules

- `record_status` remains `public_source_derived_pending_human_review`.
- `source_*` columns identify the API release, retrieval time, study, original row order, returned counts, and completeness.
- `study_locus_id`, variant IDs, Ensembl gene IDs, QTL study IDs, and biosample IDs are preserved.
- `source_quality_controls` carries the exact unresolved source warning returned for each locus: `Study has quality control flag(s)`. The snapshot does not explain the flag, so this project does not guess its meaning.
- Empty cells mean missing, absent, or not applicable in the saved source record. They must not be treated automatically as zero.
- Repeated locus and provenance fields make each CSV row understandable on its own.
- The tables are a five-locus descriptive extract. They cannot measure model accuracy or generalise across diseases, cohorts, or ancestries.

## Fields shared by all three files

| Field | Meaning |
|---|---|
| `record_status` | Review status carried from the evidence snapshot. |
| `source_platform_name` | Exact platform label returned with the snapshot. |
| `source_data_version_*` | Open Targets data-version components. A blank iteration means none was returned. |
| `source_api_version_*` | Open Targets API-version components. |
| `source_retrieved_utc` | UTC time when the source snapshot was retrieved. |
| `source_study_id` | FinnGen GWAS study identifier. |
| `source_project_id` | Source project identifier. |
| `source_trait` | Trait name supplied by the source study. |
| `coordinate_build` | Genome coordinate build (`GRCh38`). |
| `source_evidence_sha256` | SHA-256 of the exact `evidence.json` input. |
| `study_locus_id` | Open Targets identifier for the selected GWAS locus. |
| `lead_variant_id` | Lead variant for that selected locus. |
| `locus_chromosome`, `locus_position` | Lead-locus genomic coordinate. |
| `locus_p_value_mantissa`, `locus_p_value_exponent` | Source p-value stored as mantissa × 10^exponent to avoid underflow. |
| `locus_finemapping_method` | Fine-mapping method named by the source. |
| `locus_confidence` | Source credible-set confidence description. |
| `source_quality_controls` | Exact source QC warning text; multiple strings would be separated by `; `. |
| `source_url` | Open Targets page for the source credible set. |

Each file also includes `source_row_number_within_locus`, `source_reported_count`, `source_fetched_count`, and `source_collection_complete`. These preserve row order and document whether every reported row for that endpoint was fetched.

## `variants.csv`

| Field | Meaning |
|---|---|
| `variant_id` | Open Targets variant ID, formatted as chromosome_position_reference_alternate. |
| `variant_rs_ids` | Returned dbSNP rs identifier(s), separated by `; ` if needed. |
| `variant_chromosome`, `variant_position` | Variant coordinate. |
| `posterior_probability` | Source fine-mapping PIP. |
| `is_95_credible_set`, `is_99_credible_set` | Source Boolean membership flags. |
| `beta`, `standard_error` | Source GWAS effect estimate and standard error. |
| `locus_fetched_pip_sum` | Sum stored in the evidence snapshot for fetched PIPs at this locus. |
| `locus_missing_pip_count` | Number of fetched locus rows with a missing PIP. |

## `l2g_candidates.csv`

| Field | Meaning |
|---|---|
| `target_gene_id`, `target_gene_symbol`, `target_biotype` | Returned candidate-gene identifiers and type. |
| `canonical_transcript_*` | Source canonical-transcript coordinate and strand. |
| `l2g_score` | Original Open Targets L2G score. It is not renormalised. |
| `tss_distance_bp` | Absolute lead-variant distance to the returned transcript's transcription start site, derived when the evidence snapshot was built. |
| `candidate_coloc_status` | Whether a matching QTL row was returned for this candidate. |
| `candidate_coloc_rows_fetched` | Number of matching colocalisation rows found during source assembly. |
| `candidate_max_h4_fetched` | Maximum source H4 among those matching rows; blank when no matching row was returned. |
| `feature_<source name>_value` | Exact value of a named input feature returned with the source L2G row. |
| `feature_<source name>_shap_value` | Exact source SHAP contribution for that feature. |

The 31 exact source feature names are:

- `eQtlColocClppMaximum`
- `pQtlColocClppMaximum`
- `sQtlColocClppMaximum`
- `eQtlColocH4Maximum`
- `pQtlColocH4Maximum`
- `sQtlColocH4Maximum`
- `eQtlColocClppMaximumNeighbourhood`
- `pQtlColocClppMaximumNeighbourhood`
- `sQtlColocClppMaximumNeighbourhood`
- `eQtlColocH4MaximumNeighbourhood`
- `pQtlColocH4MaximumNeighbourhood`
- `sQtlColocH4MaximumNeighbourhood`
- `distanceSentinelFootprint`
- `distanceSentinelFootprintNeighbourhood`
- `distanceFootprintMean`
- `distanceFootprintMeanNeighbourhood`
- `distanceTssMean`
- `distanceTssMeanNeighbourhood`
- `distanceSentinelTss`
- `distanceSentinelTssNeighbourhood`
- `vepMaximum`
- `vepMaximumNeighbourhood`
- `vepMean`
- `vepMeanNeighbourhood`
- `e2gMean`
- `e2gMeanNeighbourhood`
- `geneCount500kb`
- `proteinGeneCount500kb`
- `credibleSetConfidence`
- `transPQtlColocH4Maximum`
- `transPQtlColocH4MaximumNeighbourhood`

## `molecular_qtl_colocalisations.csv`

| Field | Meaning |
|---|---|
| `colocalisation_study_locus_id` | GWAS-side locus ID repeated by the source colocalisation row. |
| `other_study_locus_id` | Molecular-QTL-side locus ID. |
| `qtl_gene_id` | Direct QTL gene field; blank in all 44 saved rows. The separate source target ID is retained below. |
| `qtl_study_id`, `qtl_study_type`, `right_study_type` | Identifiers and QTL types returned by the source. |
| `qtl_target_gene_id`, `qtl_target_gene_symbol` | Gene attached to the molecular-QTL study. |
| `biosample_id`, `biosample_name` | Tissue or cell-context identifier and label. |
| `colocalisation_method` | Source method name. |
| `h3`, `h4`, `clpp` | Separate source colocalisation statistics. |
| `number_colocalising_variants` | Source count of colocalising variants. |
| `beta_ratio_sign_average` | Source direction-summary field; preserved without reinterpretation. |

## Reproduce and verify

From the repository root:

```bash
python analysis/export_real_dataset.py
python analysis/export_real_dataset.py --check
```

The `--check` command rebuilds all five files in a temporary directory and fails if any committed output differs byte-for-byte.
