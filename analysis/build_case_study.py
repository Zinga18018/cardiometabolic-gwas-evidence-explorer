"""Build a deterministic descriptive case study from the saved evidence snapshot.

This module deliberately does not combine PIP, L2G, H4, or CLPP into a new
score or posterior. They are source measures with different meanings and
dependencies. The outputs describe evidence availability and preserve missing
values for later human review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "reports" / "evidence.json"
DEFAULT_REPORT_DIR = ROOT / "reports"
STATUS = "public_source_derived_pending_human_review"
ANALYSIS_TYPE = "descriptive_real_data_case_study"
OUTPUT_FILES = (
    "real_data_case_study.md",
    "locus_summary.csv",
    "candidate_evidence_long.csv",
    "evidence_availability.svg",
    "case_study_metadata.json",
)

LOCUS_FIELDS = (
    "status",
    "study_id",
    "study_locus_id",
    "lead_variant",
    "chromosome",
    "position",
    "finemapping_method",
    "confidence",
    "variant_rows",
    "fetched_pip_sum",
    "missing_pip_count",
    "maximum_returned_variant_pip",
    "source_95_credible_set_rows",
    "source_99_credible_set_rows",
    "l2g_candidates",
    "top_l2g_gene",
    "top_l2g_score",
    "second_l2g_score",
    "l2g_score_gap",
    "total_colocalisation_rows",
    "candidate_matched_colocalisation_rows",
    "unmatched_colocalisation_rows",
    "maximum_returned_h4",
    "maximum_returned_clpp",
    "candidate_matched_maximum_h4",
    "candidate_matched_maximum_clpp",
    "qtl_types_returned",
    "biosamples_returned",
    "source_quality_controls",
    "variants_complete",
    "candidates_complete",
    "colocalisations_complete",
    "source_url",
)

CANDIDATE_FIELDS = (
    "status",
    "record_type",
    "candidate_match_status",
    "study_id",
    "study_locus_id",
    "lead_variant",
    "gene_id",
    "gene_symbol",
    "l2g_rank",
    "l2g_score",
    "tss_distance_bp",
    "candidate_coloc_status",
    "matched_colocalisation_rows",
    "source_colocalisation_index",
    "qtl_gene_id",
    "qtl_target_gene_id",
    "qtl_target_gene_symbol",
    "qtl_study_id",
    "qtl_type",
    "biosample_id",
    "biosample_name",
    "colocalisation_method",
    "h3",
    "h4",
    "clpp",
    "number_colocalising_variants",
    "beta_ratio_sign_average",
    "source_quality_controls",
    "candidate_list_complete",
    "colocalisation_list_complete",
    "source_url",
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def probability(value, label: str):
    """Validate a source probability while retaining a missing value as None."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a number or null")
    result = float(value)
    if not math.isfinite(result) or not 0 <= result <= 1:
        raise ValueError(f"{label} must be finite and between zero and one")
    return result


def _rows(collection, label: str):
    if not isinstance(collection, dict) or not isinstance(collection.get("rows"), list):
        raise ValueError(f"{label} must contain a rows list")
    rows = collection["rows"]
    reported = collection.get("reported_count")
    fetched = collection.get("fetched_count")
    complete = collection.get("complete")
    if isinstance(reported, bool) or not isinstance(reported, int) or reported < 0:
        raise ValueError(f"{label} reported_count must be a non-negative integer")
    if isinstance(fetched, bool) or not isinstance(fetched, int) or fetched != len(rows):
        raise ValueError(f"{label} fetched_count does not match its rows")
    if not isinstance(complete, bool) or complete != (fetched == reported):
        raise ValueError(f"{label} completeness is inconsistent with its counts")
    return rows


def validate_evidence(evidence: dict) -> None:
    if not isinstance(evidence, dict) or evidence.get("schema_version") != 1:
        raise ValueError("Expected evidence schema_version 1")
    study = evidence.get("study") or {}
    if not isinstance(study.get("id"), str) or not study["id"]:
        raise ValueError("Evidence must identify its source study")
    loci = evidence.get("loci")
    if not isinstance(loci, list) or not loci:
        raise ValueError("Evidence must contain at least one locus")

    locus_ids = set()
    for locus in loci:
        locus_id = locus.get("studyLocusId")
        if not isinstance(locus_id, str) or not locus_id or locus_id in locus_ids:
            raise ValueError("Study-locus identifiers must be present and unique")
        locus_ids.add(locus_id)
        if locus.get("studyId") != study["id"]:
            raise ValueError("Locus study ID does not match the top-level study")
        quality_controls = locus.get("qualityControls")
        if not isinstance(quality_controls, list) or any(not isinstance(x, str) for x in quality_controls):
            raise ValueError("qualityControls must be a list of source strings")

        variants = _rows(locus.get("variants"), f"variants for {locus_id}")
        candidates = _rows(locus.get("candidates"), f"candidates for {locus_id}")
        colocs = _rows(locus.get("colocalisations"), f"colocalisations for {locus_id}")

        variant_ids = set()
        pip_sum = 0.0
        missing_pips = 0
        for variant in variants:
            variant_id = ((variant.get("variant") or {}).get("id"))
            if not isinstance(variant_id, str) or not variant_id or variant_id in variant_ids:
                raise ValueError(f"Variant identifiers must be present and unique within {locus_id}")
            variant_ids.add(variant_id)
            pip = probability(variant.get("posteriorProbability"), "variant PIP")
            if pip is None:
                missing_pips += 1
            else:
                pip_sum += pip
        if pip_sum > 1.001:
            raise ValueError(f"Variant PIPs sum above 1.001 for {locus_id}")
        if not math.isclose(float(locus.get("fetched_pip_sum")), pip_sum, rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Stored PIP sum does not match source rows for {locus_id}")
        if locus.get("missing_pip_count") != missing_pips:
            raise ValueError(f"Stored missing-PIP count does not match source rows for {locus_id}")

        target_ids = set()
        for candidate in candidates:
            target_id = ((candidate.get("target") or {}).get("id"))
            if not isinstance(target_id, str) or not target_id or target_id in target_ids:
                raise ValueError(f"Candidate target identifiers must be present and unique within {locus_id}")
            target_ids.add(target_id)
            probability(candidate.get("score"), "L2G score")
            probability(candidate.get("max_h4_fetched"), "candidate maximum H4")

        coloc_ids = set()
        for coloc in colocs:
            other = coloc.get("otherStudyLocus") or {}
            key = (other.get("studyLocusId"), coloc.get("colocalisationMethod"))
            if key in coloc_ids:
                raise ValueError(f"Colocalisation identifiers must be unique within {locus_id}")
            coloc_ids.add(key)
            for metric in ("h3", "h4", "clpp"):
                probability(coloc.get(metric), metric)


def load_evidence(path: Path):
    raw = Path(path).read_bytes()
    evidence = json.loads(raw)
    validate_evidence(evidence)
    return evidence, sha256_bytes(raw)


def quality_control_text(locus: dict) -> str:
    controls = locus.get("qualityControls") or []
    return "; ".join(controls) if controls else "none returned"


def colocalisation_matches_candidate(coloc: dict, target_id: str) -> bool:
    """Use the same two-source gene-ID matching rule as the explorer."""
    other = coloc.get("otherStudyLocus") or {}
    study_target = ((other.get("study") or {}).get("target") or {}).get("id")
    return other.get("qtlGeneId") == target_id or study_target == target_id


def _maximum(values):
    present = [value for value in values if value is not None]
    return max(present) if present else None


def _candidate_order(locus: dict):
    return sorted(
        locus["candidates"]["rows"],
        key=lambda candidate: (
            -(candidate["score"] if candidate.get("score") is not None else -1),
            candidate["target"]["id"],
        ),
    )


def summarize_locus(locus: dict) -> dict:
    variants = locus["variants"]["rows"]
    candidates = _candidate_order(locus)
    colocs = locus["colocalisations"]["rows"]
    pips = [row.get("posteriorProbability") for row in variants]
    scores = [row.get("score") for row in candidates if row.get("score") is not None]

    matched_indices = {
        index
        for index, coloc in enumerate(colocs)
        if any(colocalisation_matches_candidate(coloc, candidate["target"]["id"]) for candidate in candidates)
    }
    matched_colocs = [colocs[index] for index in sorted(matched_indices)]
    qtl_types = sorted({row.get("rightStudyType") for row in colocs if row.get("rightStudyType")})
    biosamples = sorted(
        {
            ((row.get("otherStudyLocus") or {}).get("study") or {}).get("biosample", {}).get("biosampleName")
            for row in colocs
            if (((row.get("otherStudyLocus") or {}).get("study") or {}).get("biosample") or {}).get("biosampleName")
        }
    )
    top = candidates[0] if candidates else None
    second_score = scores[1] if len(scores) > 1 else None
    top_score = scores[0] if scores else None
    return {
        "status": STATUS,
        "study_id": locus["studyId"],
        "study_locus_id": locus["studyLocusId"],
        "lead_variant": locus["variant"]["id"],
        "chromosome": locus["chromosome"],
        "position": locus["position"],
        "finemapping_method": locus.get("finemappingMethod"),
        "confidence": locus.get("confidence"),
        "variant_rows": len(variants),
        "fetched_pip_sum": locus["fetched_pip_sum"],
        "missing_pip_count": locus["missing_pip_count"],
        "maximum_returned_variant_pip": _maximum(pips),
        "source_95_credible_set_rows": sum(row.get("is95CredibleSet") is True for row in variants),
        "source_99_credible_set_rows": sum(row.get("is99CredibleSet") is True for row in variants),
        "l2g_candidates": len(candidates),
        "top_l2g_gene": (top["target"].get("approvedSymbol") if top else None),
        "top_l2g_score": top_score,
        "second_l2g_score": second_score,
        "l2g_score_gap": (top_score - second_score if top_score is not None and second_score is not None else None),
        "total_colocalisation_rows": len(colocs),
        "candidate_matched_colocalisation_rows": len(matched_indices),
        "unmatched_colocalisation_rows": len(colocs) - len(matched_indices),
        "maximum_returned_h4": _maximum(row.get("h4") for row in colocs),
        "maximum_returned_clpp": _maximum(row.get("clpp") for row in colocs),
        "candidate_matched_maximum_h4": _maximum(row.get("h4") for row in matched_colocs),
        "candidate_matched_maximum_clpp": _maximum(row.get("clpp") for row in matched_colocs),
        "qtl_types_returned": "; ".join(qtl_types),
        "biosamples_returned": "; ".join(biosamples),
        "source_quality_controls": quality_control_text(locus),
        "variants_complete": locus["variants"]["complete"],
        "candidates_complete": locus["candidates"]["complete"],
        "colocalisations_complete": locus["colocalisations"]["complete"],
        "source_url": locus["source_url"],
    }


def candidate_evidence_rows(locus: dict):
    colocs = locus["colocalisations"]["rows"]
    rows = []
    matched_source_indices = set()
    for rank, candidate in enumerate(_candidate_order(locus), start=1):
        target = candidate["target"]
        matched = [
            (index, coloc)
            for index, coloc in enumerate(colocs)
            if colocalisation_matches_candidate(coloc, target["id"])
        ]
        matched_source_indices.update(index for index, _ in matched)
        records = matched or [(None, None)]
        for index, coloc in records:
            other = (coloc or {}).get("otherStudyLocus") or {}
            study = other.get("study") or {}
            qtl_target = study.get("target") or {}
            biosample = study.get("biosample") or {}
            rows.append(
                {
                    "status": STATUS,
                    "record_type": (
                        "candidate_colocalisation_match"
                        if coloc is not None
                        else "candidate_without_matched_colocalisation"
                    ),
                    "candidate_match_status": (
                        "matched_to_returned_candidate"
                        if coloc is not None
                        else "not_returned_for_candidate"
                    ),
                    "study_id": locus["studyId"],
                    "study_locus_id": locus["studyLocusId"],
                    "lead_variant": locus["variant"]["id"],
                    "gene_id": target["id"],
                    "gene_symbol": target.get("approvedSymbol"),
                    "l2g_rank": rank,
                    "l2g_score": candidate.get("score"),
                    "tss_distance_bp": candidate.get("tss_distance_bp"),
                    "candidate_coloc_status": candidate.get("coloc_status"),
                    "matched_colocalisation_rows": len(matched),
                    "source_colocalisation_index": index,
                    "qtl_gene_id": other.get("qtlGeneId"),
                    "qtl_target_gene_id": qtl_target.get("id"),
                    "qtl_target_gene_symbol": qtl_target.get("approvedSymbol"),
                    "qtl_study_id": study.get("id"),
                    "qtl_type": (coloc or {}).get("rightStudyType"),
                    "biosample_id": biosample.get("biosampleId"),
                    "biosample_name": biosample.get("biosampleName"),
                    "colocalisation_method": (coloc or {}).get("colocalisationMethod"),
                    "h3": (coloc or {}).get("h3"),
                    "h4": (coloc or {}).get("h4"),
                    "clpp": (coloc or {}).get("clpp"),
                    "number_colocalising_variants": (coloc or {}).get("numberColocalisingVariants"),
                    "beta_ratio_sign_average": (coloc or {}).get("betaRatioSignAverage"),
                    "source_quality_controls": quality_control_text(locus),
                    "candidate_list_complete": locus["candidates"]["complete"],
                    "colocalisation_list_complete": locus["colocalisations"]["complete"],
                    "source_url": locus["source_url"],
                }
            )
    # Retain every source colocalisation that did not match a returned L2G
    # candidate. It appears once at the locus level rather than being copied to
    # each unrelated candidate.
    for index, coloc in enumerate(colocs):
        if index in matched_source_indices:
            continue
        other = coloc.get("otherStudyLocus") or {}
        study = other.get("study") or {}
        qtl_target = study.get("target") or {}
        biosample = study.get("biosample") or {}
        rows.append(
            {
                "status": STATUS,
                "record_type": "unmatched_locus_colocalisation",
                "candidate_match_status": "unmatched_to_returned_candidates",
                "study_id": locus["studyId"],
                "study_locus_id": locus["studyLocusId"],
                "lead_variant": locus["variant"]["id"],
                "gene_id": None,
                "gene_symbol": None,
                "l2g_rank": None,
                "l2g_score": None,
                "tss_distance_bp": None,
                "candidate_coloc_status": None,
                "matched_colocalisation_rows": None,
                "source_colocalisation_index": index,
                "qtl_gene_id": other.get("qtlGeneId"),
                "qtl_target_gene_id": qtl_target.get("id"),
                "qtl_target_gene_symbol": qtl_target.get("approvedSymbol"),
                "qtl_study_id": study.get("id"),
                "qtl_type": coloc.get("rightStudyType"),
                "biosample_id": biosample.get("biosampleId"),
                "biosample_name": biosample.get("biosampleName"),
                "colocalisation_method": coloc.get("colocalisationMethod"),
                "h3": coloc.get("h3"),
                "h4": coloc.get("h4"),
                "clpp": coloc.get("clpp"),
                "number_colocalising_variants": coloc.get("numberColocalisingVariants"),
                "beta_ratio_sign_average": coloc.get("betaRatioSignAverage"),
                "source_quality_controls": quality_control_text(locus),
                "candidate_list_complete": locus["candidates"]["complete"],
                "colocalisation_list_complete": locus["colocalisations"]["complete"],
                "source_url": locus["source_url"],
            }
        )
    return rows


def summarize_case(evidence: dict):
    locus_rows = [summarize_locus(locus) for locus in evidence["loci"]]
    candidate_rows = [row for locus in evidence["loci"] for row in candidate_evidence_rows(locus)]
    counts = {
        "loci": len(locus_rows),
        "variant_rows": sum(row["variant_rows"] for row in locus_rows),
        "l2g_candidates": sum(row["l2g_candidates"] for row in locus_rows),
        "colocalisation_rows": sum(row["total_colocalisation_rows"] for row in locus_rows),
        "candidate_matched_colocalisation_rows": sum(
            row["candidate_matched_colocalisation_rows"] for row in locus_rows
        ),
        "unmatched_colocalisation_rows": sum(row["unmatched_colocalisation_rows"] for row in locus_rows),
        "candidate_long_rows": len(candidate_rows),
    }
    return locus_rows, candidate_rows, counts


def _write_csv(path: Path, fields, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_text_lf(path: Path, value: str) -> None:
    """Write stable LF bytes on Windows, macOS, and Linux."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(value)


def _display(value, digits=6):
    if value is None or value == "":
        return "not returned"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def markdown_report(evidence: dict, locus_rows: list[dict], counts: dict) -> str:
    release = evidence["platform_meta"]["dataVersion"]
    selection = evidence["selection"]
    lines = [
        "# Real-data case study: evidence heterogeneity across FinnGen T2D loci",
        "",
        f"**Status:** `{STATUS}`",
        "",
        "## Research question",
        "",
        "How do fine-mapping concentration, Locus-to-Gene rankings, molecular-QTL colocalisation, and evidence availability differ across five selected FinnGen R12 type 2 diabetes loci?",
        "",
        "## Data and scope",
        "",
        f"This descriptive case study uses the saved Open Targets {release['year']}.{release['month']} snapshot for `{evidence['study']['id']}` ({evidence['study']['traitFromSource']}). The source contains {selection['discovery_count']} credible sets; {len(locus_rows)} were selected by the documented p-value and 1 Mb spacing rule. It includes {counts['variant_rows']} variant rows, {counts['l2g_candidates']} returned L2G candidates, and {counts['colocalisation_rows']} molecular-QTL colocalisation rows.",
        "",
        "These are public summary-level source records, not participant-level data. PIP, L2G, H4, and CLPP are reported separately. No composite score, causal posterior, core-gene posterior, therapeutic ranking, or model-performance estimate is calculated.",
        "",
        "## Locus-level results",
        "",
        "| Lead variant | Top returned L2G gene | Variant rows | PIP sum | Maximum returned PIP | L2G candidates | Total QTL colocs | Candidate-matched colocs | Locus-wide max H4 | Returned-candidate max H4 | Locus-wide max CLPP | Returned-candidate max CLPP |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in locus_rows:
        lines.append(
            f"| [{row['lead_variant']}]({row['source_url']}) | {row['top_l2g_gene'] or 'not returned'} | "
            f"{row['variant_rows']} | {_display(row['fetched_pip_sum'])} | {_display(row['maximum_returned_variant_pip'])} | "
            f"{row['l2g_candidates']} | {row['total_colocalisation_rows']} | "
            f"{row['candidate_matched_colocalisation_rows']} | {_display(row['maximum_returned_h4'])} | "
            f"{_display(row['candidate_matched_maximum_h4'])} | {_display(row['maximum_returned_clpp'])} | "
            f"{_display(row['candidate_matched_maximum_clpp'])} |"
        )

    lines.extend(["", "## What the five loci show", ""])
    for row in locus_rows:
        score_text = _display(row["top_l2g_score"], 4)
        statement = (
            f"- **{row['top_l2g_gene'] or row['lead_variant']} / {row['lead_variant']}:** "
            f"{row['variant_rows']} returned variant row(s) carried PIP mass {_display(row['fetched_pip_sum'])}; "
            f"Open Targets returned {row['l2g_candidates']} L2G candidate(s), with top score {score_text}; "
            f"{row['total_colocalisation_rows']} molecular-QTL row(s) were returned and "
            f"{row['candidate_matched_colocalisation_rows']} matched a returned candidate by gene ID."
        )
        lines.append(statement)
        if row["l2g_candidates"] == 1:
            lines.append("  A single returned candidate provides no within-locus ranking comparison.")
        if row["unmatched_colocalisation_rows"]:
            lines.append(
                f"  {row['unmatched_colocalisation_rows']} colocalisation row(s) targeted genes outside the returned L2G candidate list; they were retained as locus evidence and were not reassigned."
            )
        if row["total_colocalisation_rows"] == 0:
            lines.append("  No colocalisation row was returned; this is missing source evidence, not a zero or negative biological result.")

    lines.extend(
        [
            "",
            "## Transformations performed",
            "",
            "1. Read only the committed `reports/evidence.json` snapshot and validated its schema, identifiers, counts, completeness fields, probability ranges, and PIP sums.",
            "2. Preserved source PIPs, L2G scores, H3, H4, and CLPP without normalising or combining them.",
            "3. Ranked the returned L2G candidates by their existing source score, with target ID only as a deterministic tie-breaker.",
            "4. Matched a colocalisation row to a candidate when either `otherStudyLocus.qtlGeneId` or `otherStudyLocus.study.target.id` equalled the candidate Ensembl ID, matching the explorer's existing rule.",
            "5. Retained QTL type, tissue, study, method, quality-control text, completeness fields, and missing values in the long table.",
            "6. Drew an availability chart whose counts indicate returned records, not evidence strength.",
            "",
            "## Interpretation boundaries",
            "",
            "- Variant PIP, L2G score, coloc H4, and eCAVIAR CLPP answer different questions and are not interchangeable probabilities.",
            "- L2G already uses features that can include distance and colocalisation. Combining the displayed values again would risk double-counting dependent evidence.",
            "- A maximum across tissues or studies is descriptive and subject to multiplicity; it is not independent replication.",
            "- A missing colocalisation row remains `not_returned`. It is not converted to zero and is not evidence against a gene.",
            "- Returned L2G candidates are an upstream-filtered set. Their scores are not renormalised, and a one-candidate locus does not validate ranking accuracy.",
            "- Every selected locus carries the exact unresolved source warning `Study has quality control flag(s)`. The snapshot does not explain its meaning, so this analysis does not speculate about its effect.",
            "- The five loci are a deterministic convenience sample from one Finnish study. They do not estimate performance across loci, cohorts, ancestries, or diseases.",
            "- The analysis does not run a new GWAS, fine mapping, Mendelian randomisation, causal-gene model, CoreGene posterior, or clinical target validation.",
            "- All derived outputs remain public-source-derived and pending human review.",
            "",
            "## Reproduce",
            "",
            "```bash",
            "python analysis/build_case_study.py",
            "python analysis/build_case_study.py --check",
            "python -m unittest discover -s tests -v",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def availability_svg(locus_rows: list[dict]) -> str:
    columns = (
        ("variant_rows", "Variant rows"),
        ("l2g_candidates", "L2G candidates"),
        ("total_colocalisation_rows", "All QTL colocs"),
        ("candidate_matched_colocalisation_rows", "Matched QTL colocs"),
        ("source_quality_controls", "Source QC flags"),
    )
    width, left, top, cell_w, row_h = 1120, 275, 150, 155, 54
    height = top + row_h * len(locus_rows) + 125
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Evidence availability across five FinnGen type 2 diabetes loci</title>',
        '<desc id="desc">A table-like chart of returned variant, L2G, total colocalisation, candidate-matched colocalisation, and source quality-control records. Counts describe availability, not evidence strength.</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111827}.title{font-size:22px;font-weight:700}.subtitle{font-size:13px;fill:#4b5563}.head{font-size:12px;font-weight:700}.label{font-size:13px;font-weight:700}.small{font-size:11px;fill:#4b5563}.value{font-size:16px;font-weight:700}.foot{font-size:12px;fill:#374151}</style>',
        '<text id="title-text" x="36" y="38" class="title">Evidence availability across five FinnGen T2D loci</text>',
        '<text x="36" y="64" class="subtitle">Counts show returned public-source records. They do not measure biological support or causal strength.</text>',
    ]
    for column_index, (_, label) in enumerate(columns):
        x = left + column_index * cell_w + cell_w / 2
        pieces.append(f'<text x="{x}" y="112" text-anchor="middle" class="head">{html.escape(label)}</text>')
    for row_index, row in enumerate(locus_rows):
        y = top + row_index * row_h
        label = f"{row['top_l2g_gene'] or 'No L2G'} · {row['lead_variant']}"
        pieces.append(f'<text x="36" y="{y + 21}" class="label">{html.escape(label)}</text>')
        pieces.append(f'<text x="36" y="{y + 39}" class="small">study locus {html.escape(row["study_locus_id"][:10])}…</text>')
        values = (
            row["variant_rows"],
            row["l2g_candidates"],
            row["total_colocalisation_rows"],
            row["candidate_matched_colocalisation_rows"],
            0 if row["source_quality_controls"] == "none returned" else len(row["source_quality_controls"].split("; ")),
        )
        for column_index, value in enumerate(values):
            x = left + column_index * cell_w
            fill = "#dbeafe" if value else "#f3f4f6"
            stroke = "#60a5fa" if value else "#d1d5db"
            if column_index == 4 and value:
                fill, stroke = "#fef3c7", "#f59e0b"
            pieces.append(f'<rect x="{x + 7}" y="{y}" width="{cell_w - 14}" height="42" rx="7" fill="{fill}" stroke="{stroke}"/>')
            pieces.append(f'<text x="{x + cell_w / 2}" y="{y + 27}" text-anchor="middle" class="value">{value}</text>')
    foot_y = top + row_h * len(locus_rows) + 20
    pieces.extend(
        [
            f'<rect x="36" y="{foot_y}" width="16" height="16" rx="3" fill="#dbeafe" stroke="#60a5fa"/>',
            f'<text x="60" y="{foot_y + 13}" class="foot">one or more records returned</text>',
            f'<rect x="250" y="{foot_y}" width="16" height="16" rx="3" fill="#f3f4f6" stroke="#d1d5db"/>',
            f'<text x="274" y="{foot_y + 13}" class="foot">no record returned</text>',
            f'<rect x="430" y="{foot_y}" width="16" height="16" rx="3" fill="#fef3c7" stroke="#f59e0b"/>',
            f'<text x="454" y="{foot_y + 13}" class="foot">unresolved source QC flag</text>',
            f'<text x="36" y="{foot_y + 50}" class="foot">Status: {STATUS}</text>',
            "</svg>",
        ]
    )
    return "\n".join(pieces) + "\n"


def build_case_study(evidence_path: Path = DEFAULT_EVIDENCE, report_dir: Path = DEFAULT_REPORT_DIR):
    evidence_path, report_dir = Path(evidence_path), Path(report_dir)
    evidence, input_hash = load_evidence(evidence_path)
    locus_rows, candidate_rows, counts = summarize_case(evidence)
    report_dir.mkdir(parents=True, exist_ok=True)

    _write_csv(report_dir / "locus_summary.csv", LOCUS_FIELDS, locus_rows)
    _write_csv(report_dir / "candidate_evidence_long.csv", CANDIDATE_FIELDS, candidate_rows)
    _write_text_lf(
        report_dir / "real_data_case_study.md",
        markdown_report(evidence, locus_rows, counts),
    )
    _write_text_lf(
        report_dir / "evidence_availability.svg",
        availability_svg(locus_rows),
    )

    output_hashes = {
        name: sha256_bytes((report_dir / name).read_bytes())
        for name in OUTPUT_FILES
        if name != "case_study_metadata.json"
    }
    release = evidence["platform_meta"]
    metadata = {
        "analysis_type": ANALYSIS_TYPE,
        "candidate_match_rule": (
            "otherStudyLocus.qtlGeneId OR otherStudyLocus.study.target.id equals candidate target.id"
        ),
        "counts": counts,
        "input": {
            "path": str(evidence_path.relative_to(ROOT)).replace("\\", "/")
            if evidence_path.is_relative_to(ROOT)
            else str(evidence_path),
            "sha256": input_hash,
        },
        "output_sha256": output_hashes,
        "principles": {
            "composite_score_created": False,
            "missing_values_converted_to_zero": False,
            "posterior_created": False,
            "source_metrics_kept_separate": ["PIP", "L2G", "H3", "H4", "CLPP"],
        },
        "schema_version": 1,
        "source": {
            "api_version": release["apiVersion"],
            "coordinate_build": evidence["coordinate_build"],
            "data_version": release["dataVersion"],
            "retrieved_utc": evidence["retrieved_utc"],
            "study_id": evidence["study"]["id"],
        },
        "status": STATUS,
    }
    _write_text_lf(
        report_dir / "case_study_metadata.json",
        json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )
    return metadata


def check_outputs(evidence_path: Path, report_dir: Path) -> list[str]:
    with tempfile.TemporaryDirectory() as temporary:
        generated = Path(temporary)
        build_case_study(evidence_path, generated)
        return [
            name
            for name in OUTPUT_FILES
            if not (Path(report_dir) / name).exists()
            or (Path(report_dir) / name).read_bytes() != (generated / name).read_bytes()
        ]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--check", action="store_true", help="verify that committed outputs regenerate byte-for-byte")
    args = parser.parse_args(argv)
    try:
        if args.check:
            changed = check_outputs(args.evidence, args.report_dir)
            if changed:
                print(json.dumps({"deterministic": False, "different_or_missing": changed}, indent=2))
                return 1
            print(json.dumps({"deterministic": True, "checked_files": list(OUTPUT_FILES)}, indent=2))
            return 0
        metadata = build_case_study(args.evidence, args.report_dir)
        print(json.dumps({"status": metadata["status"], **metadata["counts"]}, indent=2))
        return 0
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
