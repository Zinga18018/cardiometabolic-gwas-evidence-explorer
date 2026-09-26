"""A bounded, provenance-preserving explorer of existing public GWAS evidence.

Python 3.11+, standard library only. No GWAS, fine mapping, or causal model is fit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import platform
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
STUDY_ID = "FINNGEN_R12_T2D"
PAGE_SIZE = 500
MAX_PAGES = 8
META_QUERY = "{meta{name dataVersion{year month iteration} apiVersion{x y z}}}"
SCHEMA_QUERY = """{__schema{types{kind name fields{name type{kind name ofType{kind name ofType{kind name}}}} enumValues{name}}}}"""
DISCOVERY_QUERY = """query($study:String!,$index:Int!){study(studyId:$study){
 id traitFromSource nSamples nCases nControls projectId hasSumstats
 diseases{id name} discoverySamples{sampleSize ancestry}
 credibleSets(page:{size:500,index:$index}){count rows{
 studyLocusId chromosome position pValueMantissa pValueExponent
 finemappingMethod confidence variant{id rsIds}}}}}"""
DETAIL_QUERY = """query($id:String!,$index:Int!){credibleSet(studyLocusId:$id){
 studyLocusId studyId chromosome position pValueMantissa pValueExponent
 confidence finemappingMethod qualityControls variant{id rsIds}
 l2GPredictions(page:{size:500,index:$index}){count rows{score
 target{id approvedSymbol biotype canonicalTranscript{start end strand chromosome}}
 features{name value shapValue}}}
 locus(page:{size:500,index:$index}){count rows{posteriorProbability is95CredibleSet
 is99CredibleSet beta standardError variant{id rsIds chromosome position}}}
 colocalisation(studyTypes:[eqtl,pqtl,sqtl,sceqtl,scpqtl,scsqtl,tuqtl,sctuqtl],
 page:{size:500,index:$index}){count rows{studyLocusId rightStudyType h3 h4 clpp
 colocalisationMethod numberColocalisingVariants betaRatioSignAverage
 otherStudyLocus{studyLocusId qtlGeneId study{id studyType
 target{id approvedSymbol} biosample{biosampleId biosampleName}}}}}}}"""


def utc_now():
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def parse_graphql(raw):
    value = json.loads(raw)
    if not isinstance(value, dict) or value.get("errors"):
        raise ValueError("GraphQL returned errors; partial data must not be analysed: " + str(value.get("errors") if isinstance(value, dict) else value))
    if not isinstance(value.get("data"), dict):
        raise ValueError("GraphQL response has no data object")
    return value["data"]


class Snapshot:
    """Raw response bytes and the exact request are saved before analysis."""
    def __init__(self, root, create=False):
        self.root = Path(root)
        self.manifest_path = self.root / "manifest.json"
        if create:
            if self.manifest_path.exists():
                raise ValueError("Snapshot already exists. Build offline, or choose a NEW --data-dir for a fresh fetch.")
            self.root.mkdir(parents=True, exist_ok=True)
            self.manifest = {"schema_version": 1, "created_utc": utc_now(), "endpoint": ENDPOINT,
                             "python_version": platform.python_version(), "requests": []}
            self.save()
        else:
            self.manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))

    def save(self):
        write_json(self.manifest_path, self.manifest)

    def fetch(self, key, query, variables=None):
        payload = {"query": query, "variables": variables or {}}
        request_bytes = json.dumps(payload, sort_keys=True).encode("utf-8")
        request = urllib.request.Request(ENDPOINT, data=request_bytes,
                                        headers={"Content-Type": "application/json", "User-Agent": "CardiometabolicEvidenceExplorer/0.1"})
        started = utc_now()
        with urllib.request.urlopen(request, timeout=60) as response:
            raw = response.read()
            status = response.status
        filename = f"raw/{key}.json"
        output = self.root / filename
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        record = {"key": key, "file": filename, "request": payload, "request_sha256": digest(request_bytes),
                  "response_sha256": digest(raw), "bytes": len(raw), "http_status": status,
                  "requested_utc": started, "retrieved_utc": utc_now(), "source_url": ENDPOINT}
        self.manifest["requests"].append(record)
        self.save()
        return parse_graphql(raw)

    def load(self, key):
        matches = [r for r in self.manifest["requests"] if r["key"] == key]
        if len(matches) != 1:
            raise ValueError(f"Expected one cached response for {key}; found {len(matches)}")
        record = matches[0]
        raw = (self.root / record["file"]).read_bytes()
        if digest(raw) != record["response_sha256"]:
            raise ValueError(f"SHA-256 mismatch: {record['file']}")
        request_bytes = json.dumps(record["request"], sort_keys=True).encode("utf-8")
        if digest(request_bytes) != record["request_sha256"]:
            raise ValueError(f"Request SHA-256 mismatch: {key}")
        return parse_graphql(raw)

    def verify(self):
        keys = [r["key"] for r in self.manifest["requests"]]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate cache keys")
        for key in keys:
            self.load(key)
        return len(keys)


def probability(value, label):
    """Preserve null as missing. A numeric zero remains a measured zero."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"Invalid {label}: {value!r}")
    return float(value)


def quality_controls(value):
    """Preserve source quality-control text exactly and reject malformed values."""
    if value is None:
        return None
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError(f"Invalid source qualityControls field: {value!r}")
    return value


def quality_control_text(value):
    controls = quality_controls(value)
    if controls is None:
        return "not returned"
    return "; ".join(controls) if controls else "none returned"


def choose_loci(rows, limit=5, min_separation=1_000_000):
    ids = [r["studyLocusId"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate study-locus identifiers in discovery response")
    for row in rows:
        m, e = row["pValueMantissa"], row["pValueExponent"]
        if not 1 <= m < 10 or not isinstance(e, int):
            raise ValueError("Invalid scientific-notation p-value")
    # Comparing scientific notation avoids converting very small p-values to zero.
    ordered = sorted(rows, key=lambda r: (r["pValueExponent"], r["pValueMantissa"], r["studyLocusId"]))
    selected = []
    for row in ordered:
        if all(row["chromosome"] != previous["chromosome"] or
               abs(row["position"] - previous["position"]) >= min_separation for previous in selected):
            selected.append(row)
        if len(selected) == limit:
            break
    return selected


def tss_distance(target, chromosome, position):
    transcript = target.get("canonicalTranscript")
    if not transcript or transcript.get("chromosome") != chromosome:
        return None
    if transcript.get("strand") == "POSITIVE":
        coordinate = transcript.get("start")
    elif transcript.get("strand") == "NEGATIVE":
        coordinate = transcript.get("end")
    else:
        return None
    return abs(position - coordinate) if isinstance(coordinate, int) else None


def merge_pages(pages, field, identity):
    expected = pages[0][field]["count"]
    merged, seen = [], set()
    for page in pages:
        if page[field]["count"] != expected:
            raise ValueError(f"Count changed during pagination: {field}")
        for row in page[field]["rows"]:
            key = identity(row)
            if key in seen:
                raise ValueError(f"Duplicate across pages: {field}: {key}")
            seen.add(key)
            merged.append(row)
    if len(merged) > expected:
        raise ValueError(f"More rows than advertised: {field}")
    return {"rows": merged, "reported_count": expected, "fetched_count": len(merged), "complete": len(merged) == expected}


def fetch_snapshot(data_dir, limit=5):
    snapshot = Snapshot(data_dir, create=True)
    snapshot.fetch("meta_start", META_QUERY)
    snapshot.fetch("schema", SCHEMA_QUERY)
    discovery_pages = []
    for index in range(MAX_PAGES):
        study = snapshot.fetch(f"discovery_{index}", DISCOVERY_QUERY, {"study": STUDY_ID, "index": index})["study"]
        if not study:
            raise ValueError("Configured study was not returned")
        discovery_pages.append(study)
        if (index + 1) * PAGE_SIZE >= study["credibleSets"]["count"]:
            break
    discovery = merge_pages(discovery_pages, "credibleSets", lambda r: r["studyLocusId"])
    if not discovery["complete"]:
        raise ValueError("Discovery exceeds bounded pagination; selection would be biased by truncation")
    selected = choose_loci(discovery["rows"], limit=limit)
    if len(selected) < limit:
        raise ValueError("Fewer eligible loci than requested")
    snapshot.manifest["study_id"] = STUDY_ID
    snapshot.manifest["selection"] = {"count": limit, "rule": "ascending p-value; ties studyLocusId; same-chromosome lead variants at least 1000000 bp apart",
                                       "min_separation_bp": 1_000_000, "selected_ids": [r["studyLocusId"] for r in selected],
                                       "discovery_count": discovery["reported_count"]}
    snapshot.save()
    for row in selected:
        locus_id = row["studyLocusId"]
        for index in range(MAX_PAGES):
            detail = snapshot.fetch(f"locus_{locus_id}_{index}", DETAIL_QUERY, {"id": locus_id, "index": index})["credibleSet"]
            if not detail or detail["studyLocusId"] != locus_id or detail["studyId"] != STUDY_ID:
                raise ValueError("Study-locus identity mismatch")
            if (index + 1) * PAGE_SIZE >= max(detail[k]["count"] for k in ("locus", "l2GPredictions", "colocalisation")):
                break
    snapshot.fetch("meta_end", META_QUERY)
    if snapshot.load("meta_start")["meta"] != snapshot.load("meta_end")["meta"]:
        raise ValueError("Platform release changed during collection; keep this snapshot for inspection and fetch a new one")
    snapshot.manifest["completed_utc"] = utc_now()
    snapshot.save()
    return snapshot


def compile_evidence(snapshot):
    snapshot.verify()
    if not snapshot.manifest.get("completed_utc"):
        raise ValueError("Snapshot is incomplete")
    study = dict(snapshot.load("discovery_0")["study"])
    study.pop("credibleSets")
    records = []
    for locus_id in snapshot.manifest["selection"]["selected_ids"]:
        keys = sorted((r["key"] for r in snapshot.manifest["requests"] if r["key"].startswith(f"locus_{locus_id}_")),
                      key=lambda key: int(key.rsplit("_", 1)[1]))
        pages = [snapshot.load(key)["credibleSet"] for key in keys]
        row = {k: v for k, v in pages[0].items() if k not in ("locus", "l2GPredictions", "colocalisation")}
        row["qualityControls"] = quality_controls(row.get("qualityControls"))
        for page in pages:
            if page["studyLocusId"] != locus_id or page["studyId"] != study["id"]:
                raise ValueError("Cached locus identity mismatch")
        variants = merge_pages(pages, "locus", lambda r: r["variant"]["id"])
        candidates = merge_pages(pages, "l2GPredictions", lambda r: r["target"]["id"])
        colocs = merge_pages(pages, "colocalisation", lambda r: (r["otherStudyLocus"]["studyLocusId"], r["colocalisationMethod"]))
        for variant in variants["rows"]:
            variant["posteriorProbability"] = probability(variant.get("posteriorProbability"), "variant PIP")
        known_pips = [v["posteriorProbability"] for v in variants["rows"] if v["posteriorProbability"] is not None]
        row["fetched_pip_sum"] = sum(known_pips)
        row["missing_pip_count"] = len(variants["rows"]) - len(known_pips)
        if row["fetched_pip_sum"] > 1.001:
            raise ValueError("Credible-set PIPs sum above 1.001")
        for candidate in candidates["rows"]:
            candidate["score"] = probability(candidate.get("score"), "L2G score")
            candidate["tss_distance_bp"] = tss_distance(candidate["target"], row["chromosome"], row["position"])
        candidates["rows"].sort(key=lambda x: (-(x["score"] if x["score"] is not None else -1), x["target"]["id"]))
        for coloc in colocs["rows"]:
            for metric in ("h3", "h4", "clpp"):
                coloc[metric] = probability(coloc.get(metric), metric)
        for candidate in candidates["rows"]:
            target_id = candidate["target"]["id"]
            matched = [c for c in colocs["rows"] if c["otherStudyLocus"].get("qtlGeneId") == target_id or
                       ((c["otherStudyLocus"].get("study") or {}).get("target") or {}).get("id") == target_id]
            h4s = [c["h4"] for c in matched if c["h4"] is not None]
            candidate["coloc_rows_fetched"] = len(matched)
            candidate["max_h4_fetched"] = max(h4s) if h4s else None
            candidate["coloc_status"] = ("present" if matched else "not_returned") if colocs["complete"] else "partial_fetch"
        distance_rows = [c for c in candidates["rows"] if c["tss_distance_bp"] is not None]
        min_distance = min((c["tss_distance_bp"] for c in distance_rows), default=None)
        row["nearest_returned_candidate_ids"] = [c["target"]["id"] for c in distance_rows if c["tss_distance_bp"] == min_distance]
        row["top_l2g_target"] = candidates["rows"][0]["target"]["approvedSymbol"] if candidates["rows"] else None
        row["candidates"] = candidates
        row["variants"] = variants
        row["colocalisations"] = colocs
        row["source_url"] = f"https://platform.opentargets.org/credible-set/{locus_id}"
        row["evidence_status"] = "public_source_derived_pending_human_review"
        records.append(row)
    return {"schema_version": 1, "study": study, "platform_meta": snapshot.load("meta_start")["meta"],
            "retrieved_utc": snapshot.manifest["completed_utc"], "selection": snapshot.manifest["selection"],
            "coordinate_build": "GRCh38", "loci": records}


def evidence_csv(evidence, path):
    fields = ["study_id", "study_locus_id", "lead_variant", "gene_id", "gene_symbol", "l2g_score",
              "tss_distance_bp", "nearest_among_returned_candidates", "coloc_rows_fetched", "max_h4_fetched",
              "coloc_status", "candidate_list_complete", "coloc_list_complete", "source_quality_controls",
              "source_url", "status"]
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        for locus in evidence["loci"]:
            for candidate in locus["candidates"]["rows"]:
                target = candidate["target"]
                writer.writerow(dict(zip(fields, [locus["studyId"], locus["studyLocusId"], locus["variant"]["id"],
                    target["id"], target["approvedSymbol"], candidate["score"], candidate["tss_distance_bp"],
                    target["id"] in locus["nearest_returned_candidate_ids"], candidate["coloc_rows_fetched"],
                    candidate["max_h4_fetched"], candidate["coloc_status"], locus["candidates"]["complete"],
                    locus["colocalisations"]["complete"], quality_control_text(locus.get("qualityControls")),
                    locus["source_url"], locus["evidence_status"]])))


def report_markdown(evidence):
    loci = evidence["loci"]
    version = evidence["platform_meta"]["dataVersion"]
    rows = ["# Cardiometabolic GWAS evidence audit", "", f"Snapshot: {evidence['retrieved_utc']}. Open Targets data release {version['year']}.{version['month']}.",
            "", f"Study: [{STUDY_ID}](https://platform.opentargets.org/study/{STUDY_ID}) — {evidence['study']['traitFromSource']}.",
            f"Selected {len(loci)} of {evidence['selection']['discovery_count']} credible sets by the smallest reported p-values, requiring lead variants on the same chromosome to be at least 1 Mb apart. This is a convenience sample, not an independent-locus proof or performance benchmark.",
            "", "| Lead variant | Top returned L2G gene | Variants fetched | PIP sum | L2G candidates | Molecular-QTL colocalisation rows | Source quality controls |", "|---|---|---:|---:|---:|---:|---|"]
    for locus in loci:
        rows.append(f"| [{locus['variant']['id']}]({locus['source_url']}) | {locus['top_l2g_target'] or 'not returned'} | {locus['variants']['fetched_count']} | {locus['fetched_pip_sum']:.6f} | {locus['candidates']['fetched_count']} | {locus['colocalisations']['fetched_count']} | {quality_control_text(locus.get('qualityControls'))} |")
    rows.extend(["", "## What was computed", "", "Parsed and validated source probabilities, checked pagination and cached-response hashes, preserved missing values, and calculated absolute distances from lead variants to canonical transcription start sites for returned candidates. Negative-strand transcripts use their end coordinate. L2G and colocalisation values were supplied by Open Targets; they were not estimated by this project.",
                 "", "## Interpretation boundaries", "", "- A missing colocalisation row means no row was returned by this query and release; it is not evidence against a gene.",
                 "- L2G scores are existing model outputs; they are not this project's causal posterior or clinical target probability. Colocalisation H4 and eCAVIAR CLPP are distinct fields, not interchangeable scores.",
                 "- The nearest comparison covers only returned L2G candidates, which are filtered upstream; it does not identify the nearest gene across the genome. Agreement is descriptive, not accuracy.",
                 "- Maximum H4 is the maximum among fetched molecular-QTL records. It may reflect a tissue unrelated to the disease mechanism. Multiple QTL signals and tissues are not independent replications.",
                 "- Source quality-control text is reproduced exactly from the `qualityControls` field. The cached Open Targets API response provides no additional detail about a returned flag, so this project does not infer its meaning or effect.",
                 "- No new GWAS, fine-mapping, Mendelian randomisation, drug-target validation, wet-lab experiment or AHBA analysis was performed.",
                 "- Study-level nSamples and discoverySamples are retained as separate source annotations and can differ. No denominator was silently changed.",
                 "- Every derived record remains public_source_derived_pending_human_review.", "", "## Reproduce", "", "Run `python gwas_explorer.py verify`, `python gwas_explorer.py build`, and `python -m unittest discover -s tests -v` from the repository root. A new download requires a new snapshot directory."])
    return "\n".join(rows) + "\n"


def build(snapshot, report_dir):
    evidence = compile_evidence(snapshot)
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    write_json(report_dir / "evidence.json", evidence)
    evidence_csv(evidence, report_dir / "candidate_evidence.csv")
    (report_dir / "audit_report.md").write_text(report_markdown(evidence), encoding="utf-8")
    template = Path(__file__).with_name("explorer_template.html").read_text(encoding="utf-8")
    encoded = json.dumps(evidence, ensure_ascii=True, allow_nan=False).replace("<", "\\u003c")
    (report_dir / "index.html").write_text(template.replace("__EVIDENCE_JSON__", encoded), encoding="utf-8")
    summary = {"study": evidence["study"]["id"], "loci": len(evidence["loci"]),
               "variants": sum(x["variants"]["fetched_count"] for x in evidence["loci"]),
               "l2g_candidates": sum(x["candidates"]["fetched_count"] for x in evidence["loci"]),
               "molecular_qtl_colocalisations": sum(x["colocalisations"]["fetched_count"] for x in evidence["loci"]),
               "all_pages_complete": all(x[k]["complete"] for x in evidence["loci"] for k in ("variants", "candidates", "colocalisations")),
               "cached_responses_verified": len(snapshot.manifest["requests"])}
    write_json(report_dir / "run_summary.json", summary)
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["fetch", "build", "verify"])
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--report-dir", type=Path, default=Path("reports"))
    parser.add_argument("--loci", type=int, default=5)
    args = parser.parse_args()
    if not 3 <= args.loci <= 8:
        parser.error("--loci must be between 3 and 8")
    try:
        if args.command == "fetch":
            snapshot = fetch_snapshot(args.data_dir, args.loci)
            build(snapshot, args.report_dir)
        elif args.command == "verify":
            print(json.dumps({"verified_responses": Snapshot(args.data_dir).verify()}))
        else:
            build(Snapshot(args.data_dir), args.report_dir)
    except (ValueError, KeyError, OSError, urllib.error.URLError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
