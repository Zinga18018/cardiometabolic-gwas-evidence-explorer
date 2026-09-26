"""Export the pinned public evidence snapshot as three readable CSV datasets.

The exporter uses only the Python standard library.  It preserves source rows,
identifiers, numeric measures, missing values, provenance, and quality-control
text.  It does not calculate a new score or convert absent values to zero.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import tempfile
import zipfile

from build_case_study import load_evidence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "reports" / "evidence.json"
DEFAULT_DATASET_DIR = ROOT / "reports" / "dataset"
DEFAULT_RAW_DATA_DIR = ROOT / "data"
PINNED_EVIDENCE_SHA256 = "49e73e407ba8a062c7b0a32f145233d24505dd381e675c550e5ba76267f5c821"
EXPECTED_COUNTS = {
    "variants": 65,
    "l2g_candidates": 11,
    "molecular_qtl_colocalisations": 44,
}
OUTPUT_FILES = (
    "variants.csv",
    "l2g_candidates.csv",
    "molecular_qtl_colocalisations.csv",
    "README.md",
    "open_targets_raw_snapshot.zip",
)
RAW_ARCHIVE_NAME = "open_targets_raw_snapshot.zip"
RAW_ARCHIVE_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
RAW_ARCHIVE_MODE = stat.S_IFREG | 0o644

L2G_FEATURE_NAMES = (
    "eQtlColocClppMaximum",
    "pQtlColocClppMaximum",
    "sQtlColocClppMaximum",
    "eQtlColocH4Maximum",
    "pQtlColocH4Maximum",
    "sQtlColocH4Maximum",
    "eQtlColocClppMaximumNeighbourhood",
    "pQtlColocClppMaximumNeighbourhood",
    "sQtlColocClppMaximumNeighbourhood",
    "eQtlColocH4MaximumNeighbourhood",
    "pQtlColocH4MaximumNeighbourhood",
    "sQtlColocH4MaximumNeighbourhood",
    "distanceSentinelFootprint",
    "distanceSentinelFootprintNeighbourhood",
    "distanceFootprintMean",
    "distanceFootprintMeanNeighbourhood",
    "distanceTssMean",
    "distanceTssMeanNeighbourhood",
    "distanceSentinelTss",
    "distanceSentinelTssNeighbourhood",
    "vepMaximum",
    "vepMaximumNeighbourhood",
    "vepMean",
    "vepMeanNeighbourhood",
    "e2gMean",
    "e2gMeanNeighbourhood",
    "geneCount500kb",
    "proteinGeneCount500kb",
    "credibleSetConfidence",
    "transPQtlColocH4Maximum",
    "transPQtlColocH4MaximumNeighbourhood",
)

COMMON_FIELDS = (
    "record_status",
    "source_platform_name",
    "source_data_version_year",
    "source_data_version_month",
    "source_data_version_iteration",
    "source_api_version_x",
    "source_api_version_y",
    "source_api_version_z",
    "source_retrieved_utc",
    "source_study_id",
    "source_project_id",
    "source_trait",
    "coordinate_build",
    "source_evidence_sha256",
    "study_locus_id",
    "lead_variant_id",
    "locus_chromosome",
    "locus_position",
    "locus_p_value_mantissa",
    "locus_p_value_exponent",
    "locus_finemapping_method",
    "locus_confidence",
    "source_quality_controls",
    "source_url",
)

VARIANT_FIELDS = COMMON_FIELDS + (
    "source_row_number_within_locus",
    "source_reported_count",
    "source_fetched_count",
    "source_collection_complete",
    "variant_id",
    "variant_rs_ids",
    "variant_chromosome",
    "variant_position",
    "posterior_probability",
    "is_95_credible_set",
    "is_99_credible_set",
    "beta",
    "standard_error",
    "locus_fetched_pip_sum",
    "locus_missing_pip_count",
)

L2G_BASE_FIELDS = COMMON_FIELDS + (
    "source_row_number_within_locus",
    "source_reported_count",
    "source_fetched_count",
    "source_collection_complete",
    "target_gene_id",
    "target_gene_symbol",
    "target_biotype",
    "canonical_transcript_chromosome",
    "canonical_transcript_start",
    "canonical_transcript_end",
    "canonical_transcript_strand",
    "l2g_score",
    "tss_distance_bp",
    "candidate_coloc_status",
    "candidate_coloc_rows_fetched",
    "candidate_max_h4_fetched",
)

L2G_FEATURE_FIELDS = tuple(
    field
    for name in L2G_FEATURE_NAMES
    for field in (f"feature_{name}_value", f"feature_{name}_shap_value")
)
L2G_FIELDS = L2G_BASE_FIELDS + L2G_FEATURE_FIELDS

COLOCALISATION_FIELDS = COMMON_FIELDS + (
    "source_row_number_within_locus",
    "source_reported_count",
    "source_fetched_count",
    "source_collection_complete",
    "colocalisation_study_locus_id",
    "other_study_locus_id",
    "qtl_gene_id",
    "qtl_study_id",
    "qtl_study_type",
    "right_study_type",
    "qtl_target_gene_id",
    "qtl_target_gene_symbol",
    "biosample_id",
    "biosample_name",
    "colocalisation_method",
    "h3",
    "h4",
    "clpp",
    "number_colocalising_variants",
    "beta_ratio_sign_average",
)


def _source_text(values) -> str:
    """Join a source string array without changing any individual value."""
    return "; ".join(values or [])


def _common_values(evidence: dict, locus: dict, evidence_hash: str) -> dict:
    platform = evidence["platform_meta"]
    data_version = platform["dataVersion"]
    api_version = platform["apiVersion"]
    study = evidence["study"]
    return {
        "record_status": locus["evidence_status"],
        "source_platform_name": platform["name"],
        "source_data_version_year": data_version.get("year"),
        "source_data_version_month": data_version.get("month"),
        "source_data_version_iteration": data_version.get("iteration"),
        "source_api_version_x": api_version.get("x"),
        "source_api_version_y": api_version.get("y"),
        "source_api_version_z": api_version.get("z"),
        "source_retrieved_utc": evidence["retrieved_utc"],
        "source_study_id": study["id"],
        "source_project_id": study.get("projectId"),
        "source_trait": study.get("traitFromSource"),
        "coordinate_build": evidence["coordinate_build"],
        "source_evidence_sha256": evidence_hash,
        "study_locus_id": locus["studyLocusId"],
        "lead_variant_id": locus["variant"]["id"],
        "locus_chromosome": locus["chromosome"],
        "locus_position": locus["position"],
        "locus_p_value_mantissa": locus["pValueMantissa"],
        "locus_p_value_exponent": locus["pValueExponent"],
        "locus_finemapping_method": locus["finemappingMethod"],
        "locus_confidence": locus["confidence"],
        "source_quality_controls": _source_text(locus.get("qualityControls")),
        "source_url": locus["source_url"],
    }


def _collection_values(collection: dict, row_number: int) -> dict:
    return {
        "source_row_number_within_locus": row_number,
        "source_reported_count": collection["reported_count"],
        "source_fetched_count": collection["fetched_count"],
        "source_collection_complete": collection["complete"],
    }


def _feature_map(candidate: dict) -> dict:
    features = candidate.get("features")
    if not isinstance(features, list):
        raise ValueError("Each L2G candidate must contain a features list")
    result = {}
    for feature in features:
        name = feature.get("name")
        if name in result:
            raise ValueError(f"Duplicate L2G feature name: {name}")
        result[name] = feature
    if set(result) != set(L2G_FEATURE_NAMES):
        missing = sorted(set(L2G_FEATURE_NAMES) - set(result))
        unexpected = sorted(set(result) - set(L2G_FEATURE_NAMES))
        raise ValueError(f"L2G feature schema changed; missing={missing}, unexpected={unexpected}")
    return result


def export_rows(evidence: dict, evidence_hash: str):
    variant_rows = []
    candidate_rows = []
    colocalisation_rows = []

    for locus in evidence["loci"]:
        common = _common_values(evidence, locus, evidence_hash)

        variants = locus["variants"]
        for row_number, source in enumerate(variants["rows"], start=1):
            variant = source["variant"]
            variant_rows.append(
                {
                    **common,
                    **_collection_values(variants, row_number),
                    "variant_id": variant["id"],
                    "variant_rs_ids": _source_text(variant.get("rsIds")),
                    "variant_chromosome": variant.get("chromosome"),
                    "variant_position": variant.get("position"),
                    "posterior_probability": source.get("posteriorProbability"),
                    "is_95_credible_set": source.get("is95CredibleSet"),
                    "is_99_credible_set": source.get("is99CredibleSet"),
                    "beta": source.get("beta"),
                    "standard_error": source.get("standardError"),
                    "locus_fetched_pip_sum": locus.get("fetched_pip_sum"),
                    "locus_missing_pip_count": locus.get("missing_pip_count"),
                }
            )

        candidates = locus["candidates"]
        for row_number, source in enumerate(candidates["rows"], start=1):
            target = source["target"]
            transcript = target.get("canonicalTranscript") or {}
            features = _feature_map(source)
            output = {
                **common,
                **_collection_values(candidates, row_number),
                "target_gene_id": target.get("id"),
                "target_gene_symbol": target.get("approvedSymbol"),
                "target_biotype": target.get("biotype"),
                "canonical_transcript_chromosome": transcript.get("chromosome"),
                "canonical_transcript_start": transcript.get("start"),
                "canonical_transcript_end": transcript.get("end"),
                "canonical_transcript_strand": transcript.get("strand"),
                "l2g_score": source.get("score"),
                "tss_distance_bp": source.get("tss_distance_bp"),
                "candidate_coloc_status": source.get("coloc_status"),
                "candidate_coloc_rows_fetched": source.get("coloc_rows_fetched"),
                "candidate_max_h4_fetched": source.get("max_h4_fetched"),
            }
            for name in L2G_FEATURE_NAMES:
                output[f"feature_{name}_value"] = features[name].get("value")
                output[f"feature_{name}_shap_value"] = features[name].get("shapValue")
            candidate_rows.append(output)

        colocalisations = locus["colocalisations"]
        for row_number, source in enumerate(colocalisations["rows"], start=1):
            other_locus = source.get("otherStudyLocus") or {}
            qtl_study = other_locus.get("study") or {}
            qtl_target = qtl_study.get("target") or {}
            biosample = qtl_study.get("biosample") or {}
            colocalisation_rows.append(
                {
                    **common,
                    **_collection_values(colocalisations, row_number),
                    "colocalisation_study_locus_id": source.get("studyLocusId"),
                    "other_study_locus_id": other_locus.get("studyLocusId"),
                    "qtl_gene_id": other_locus.get("qtlGeneId"),
                    "qtl_study_id": qtl_study.get("id"),
                    "qtl_study_type": qtl_study.get("studyType"),
                    "right_study_type": source.get("rightStudyType"),
                    "qtl_target_gene_id": qtl_target.get("id"),
                    "qtl_target_gene_symbol": qtl_target.get("approvedSymbol"),
                    "biosample_id": biosample.get("biosampleId"),
                    "biosample_name": biosample.get("biosampleName"),
                    "colocalisation_method": source.get("colocalisationMethod"),
                    "h3": source.get("h3"),
                    "h4": source.get("h4"),
                    "clpp": source.get("clpp"),
                    "number_colocalising_variants": source.get("numberColocalisingVariants"),
                    "beta_ratio_sign_average": source.get("betaRatioSignAverage"),
                }
            )

    return variant_rows, candidate_rows, colocalisation_rows


def _write_csv(path: Path, fieldnames: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def raw_source_members(raw_data_dir: Path) -> list[tuple[str, Path]]:
    """Return the exact manifest and response files allowed in the raw ZIP."""
    raw_data_dir = Path(raw_data_dir)
    manifest_path = raw_data_dir / "manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    requests = manifest.get("requests")
    if not isinstance(requests, list) or len(requests) != 9:
        raise ValueError("The pinned raw snapshot must contain exactly 9 API requests")

    listed_paths = []
    for request in requests:
        relative_text = request.get("file")
        if not isinstance(relative_text, str):
            raise ValueError("Every manifest request must identify a response file")
        relative = PurePosixPath(relative_text)
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or len(relative.parts) != 2
            or relative.parts[0] != "raw"
            or relative.suffix != ".json"
        ):
            raise ValueError(f"Unsafe or unexpected raw response path: {relative_text}")
        if relative_text in listed_paths:
            raise ValueError(f"Duplicate raw response path in manifest: {relative_text}")
        listed_paths.append(relative_text)

        source_path = raw_data_dir.joinpath(*relative.parts)
        payload = source_path.read_bytes()
        if request.get("bytes") != len(payload):
            raise ValueError(f"Manifest byte count does not match {relative_text}")
        payload_hash = hashlib.sha256(payload).hexdigest()
        if request.get("response_sha256") != payload_hash:
            raise ValueError(f"Manifest SHA-256 does not match {relative_text}")

    disk_paths = {
        path.relative_to(raw_data_dir).as_posix()
        for path in (raw_data_dir / "raw").glob("*.json")
        if path.is_file()
    }
    if disk_paths != set(listed_paths):
        missing = sorted(set(listed_paths) - disk_paths)
        unexpected = sorted(disk_paths - set(listed_paths))
        raise ValueError(f"Raw response files changed; missing={missing}, unexpected={unexpected}")

    members = [("data/manifest.json", manifest_path)]
    members.extend(
        (f"data/{relative}", raw_data_dir.joinpath(*PurePosixPath(relative).parts))
        for relative in listed_paths
    )
    return sorted(members, key=lambda item: item[0])


def _write_raw_snapshot_archive(raw_data_dir: Path, archive_path: Path) -> None:
    """Package raw bytes with fixed ZIP metadata for cross-platform stability."""
    members = raw_source_members(raw_data_dir)
    with zipfile.ZipFile(archive_path, mode="w", compression=zipfile.ZIP_STORED) as archive:
        for archive_name, source_path in members:
            info = zipfile.ZipInfo(archive_name, date_time=RAW_ARCHIVE_TIMESTAMP)
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = RAW_ARCHIVE_MODE << 16
            info.extra = b""
            info.comment = b""
            archive.writestr(info, source_path.read_bytes())


def _readme(evidence_hash: str) -> str:
    feature_names = "\n".join(f"- `{name}`" for name in L2G_FEATURE_NAMES)
    return f"""# Exported real-data tables

## What “real data” means here

These CSVs contain **real public aggregate records** saved from the Open Targets Platform. The disease-association side comes from the FinnGen R12 type 2 diabetes study `FINNGEN_R12_T2D`; the molecular-QTL side contains public QTL studies integrated by Open Targets.

They are **not participant-level records**, patient data, or genotypes. They are **not a complete GWAS summary-statistics download**. A row represents a returned variant, gene candidate, or molecular-QTL colocalisation record. For example, 65 variant rows means 65 variant records, not 65 people.

The export covers five deliberately selected loci from one study:

- `variants.csv`: 65 fine-mapped variant rows.
- `l2g_candidates.csv`: 11 Open Targets Locus-to-Gene candidate rows.
- `molecular_qtl_colocalisations.csv`: 44 molecular-QTL colocalisation rows.
- `open_targets_raw_snapshot.zip`: the manifest and nine original API response JSON files behind the five-locus snapshot.

The three CSVs come deterministically from `reports/evidence.json` (SHA-256 `{evidence_hash}`). No new score or probability is calculated here. Missing source values are empty CSV cells; a source value of zero remains `0` or `0.0`.

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

{feature_names}

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
"""


def _validate_pinned_snapshot(evidence: dict, evidence_hash: str, row_counts: dict) -> None:
    if evidence_hash != PINNED_EVIDENCE_SHA256:
        raise ValueError(
            "The evidence snapshot is not the pinned input: "
            f"expected {PINNED_EVIDENCE_SHA256}, got {evidence_hash}"
        )
    if evidence["study"]["id"] != "FINNGEN_R12_T2D":
        raise ValueError("The pinned dataset must use study FINNGEN_R12_T2D")
    if row_counts != EXPECTED_COUNTS:
        raise ValueError(f"Pinned row counts changed: expected {EXPECTED_COUNTS}, got {row_counts}")


def export_dataset(
    evidence_path: Path = DEFAULT_EVIDENCE,
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    raw_data_dir: Path = DEFAULT_RAW_DATA_DIR,
) -> dict:
    evidence_path = Path(evidence_path)
    dataset_dir = Path(dataset_dir)
    raw_data_dir = Path(raw_data_dir)
    evidence, evidence_hash = load_evidence(evidence_path)
    variants, candidates, colocalisations = export_rows(evidence, evidence_hash)
    row_counts = {
        "variants": len(variants),
        "l2g_candidates": len(candidates),
        "molecular_qtl_colocalisations": len(colocalisations),
    }
    _validate_pinned_snapshot(evidence, evidence_hash, row_counts)

    dataset_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(dataset_dir / "variants.csv", VARIANT_FIELDS, variants)
    _write_csv(dataset_dir / "l2g_candidates.csv", L2G_FIELDS, candidates)
    _write_csv(
        dataset_dir / "molecular_qtl_colocalisations.csv",
        COLOCALISATION_FIELDS,
        colocalisations,
    )
    _write_text_lf(dataset_dir / "README.md", _readme(evidence_hash))
    archive_path = dataset_dir / RAW_ARCHIVE_NAME
    _write_raw_snapshot_archive(raw_data_dir, archive_path)
    archive_bytes = archive_path.read_bytes()
    return {
        "input_sha256": evidence_hash,
        "output_files": list(OUTPUT_FILES),
        "raw_archive": {
            "bytes": len(archive_bytes),
            "members": [name for name, _ in raw_source_members(raw_data_dir)],
            "sha256": hashlib.sha256(archive_bytes).hexdigest(),
        },
        "row_counts": row_counts,
    }


def check_outputs(evidence_path: Path, dataset_dir: Path, raw_data_dir: Path) -> list[str]:
    with tempfile.TemporaryDirectory() as temporary:
        generated = Path(temporary)
        export_dataset(evidence_path, generated, raw_data_dir)
        return [
            name
            for name in OUTPUT_FILES
            if not (Path(dataset_dir) / name).exists()
            or (Path(dataset_dir) / name).read_bytes() != (generated / name).read_bytes()
        ]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--raw-data-dir", type=Path, default=DEFAULT_RAW_DATA_DIR)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that committed exports regenerate byte-for-byte",
    )
    args = parser.parse_args(argv)
    try:
        if args.check:
            changed = check_outputs(args.evidence, args.dataset_dir, args.raw_data_dir)
            if changed:
                print(json.dumps({"deterministic": False, "different_or_missing": changed}, indent=2))
                return 1
            print(json.dumps({"deterministic": True, "checked_files": list(OUTPUT_FILES)}, indent=2))
            return 0
        result = export_dataset(args.evidence, args.dataset_dir, args.raw_data_dir)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
