import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gwas_explorer import Snapshot, build, choose_loci, digest, merge_pages, parse_graphql, probability, tss_distance


class EvidenceValidationTests(unittest.TestCase):
    def test_graphql_partial_errors_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "partial data"):
            parse_graphql('{"data":{"study":{}},"errors":[{"message":"resolver failed"}]}')

    def test_missing_probability_is_not_zero(self):
        self.assertIsNone(probability(None, "h4"))
        self.assertEqual(probability(0.0, "h4"), 0.0)

    def test_invalid_probabilities_rejected(self):
        for value in [1.01, -0.1, float("nan"), float("inf"), True, "0.9"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                probability(value, "PIP")

    def test_negative_strand_tss_uses_end(self):
        target = {"canonicalTranscript":{"start":100,"end":200,"strand":"NEGATIVE","chromosome":"2"}}
        self.assertEqual(tss_distance(target,"2",180),20)

    def test_positive_strand_tss_uses_start(self):
        target = {"canonicalTranscript":{"start":100,"end":200,"strand":"POSITIVE","chromosome":"2"}}
        self.assertEqual(tss_distance(target,"2",180),80)

    def test_unknown_or_mismatched_coordinates_do_not_get_distance(self):
        self.assertIsNone(tss_distance({},"2",180))
        target = {"canonicalTranscript":{"start":100,"end":200,"strand":"UNKNOWN","chromosome":"2"}}
        self.assertIsNone(tss_distance(target,"2",180))
        target["canonicalTranscript"]["strand"] = "POSITIVE"
        self.assertIsNone(tss_distance(target,"3",180))

    def test_selection_avoids_pvalue_underflow(self):
        rows = [{"studyLocusId":name,"chromosome":chrom,"position":100,"pValueMantissa":1.5,"pValueExponent":exp}
                for name,chrom,exp in [("weaker","1",-400),("stronger","2",-500)]]
        self.assertEqual(choose_loci(rows,limit=1)[0]["studyLocusId"],"stronger")

    def test_separation_is_chromosome_aware(self):
        rows = [{"studyLocusId":name,"chromosome":chrom,"position":pos,"pValueMantissa":1.2,"pValueExponent":exp}
                for name,chrom,pos,exp in [("first","1",100,-30),("nearby","1",150,-29),("otherchr","2",150,-28),("distant","1",1000100,-27)]]
        self.assertEqual([r["studyLocusId"] for r in choose_loci(rows,limit=3)],["first","otherchr","distant"])

    def test_duplicate_discovery_ids_rejected(self):
        r={"studyLocusId":"same","chromosome":"1","position":100,"pValueMantissa":1.2,"pValueExponent":-10}
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            choose_loci([r,r])

    def test_partial_pagination_is_visible(self):
        pages=[{"locus":{"count":2,"rows":[{"id":"one"}]}}]
        result=merge_pages(pages,"locus",lambda r:r["id"])
        self.assertFalse(result["complete"])
        self.assertEqual(result["fetched_count"],1)

    def test_duplicate_pages_rejected(self):
        page={"locus":{"count":2,"rows":[{"id":"one"}]}}
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            merge_pages([page,page],"locus",lambda r:r["id"])

    def test_changed_pagination_counts_rejected(self):
        with self.assertRaisesRegex(ValueError,"Count changed"):
            merge_pages([{"x":{"count":2,"rows":[]}},{"x":{"count":3,"rows":[]}}],"x",lambda r:r["id"])

    def test_response_hash_detects_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=Snapshot(tmp,create=True)
            payload={"query":"{meta{name}}","variables":{}}
            raw=b'{"data":{"meta":{"name":"source"}}}'
            p=Path(tmp)/"raw.json";p.write_bytes(raw)
            s.manifest["requests"]=[{"key":"one","file":"raw.json","request":payload,
                "request_sha256":digest(json.dumps(payload,sort_keys=True).encode()),"response_sha256":digest(raw)}]
            self.assertEqual(s.load("one")["meta"]["name"],"source")
            p.write_bytes(raw+b" ")
            with self.assertRaisesRegex(ValueError,"SHA-256 mismatch"):
                s.load("one")

    def test_existing_snapshot_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            Snapshot(tmp,create=True)
            with self.assertRaisesRegex(ValueError,"already exists"):
                Snapshot(tmp,create=True)

    def test_source_quality_controls_surface_in_generated_outputs(self):
        root = Path(__file__).resolve().parents[1]
        expected = "Study has quality control flag(s)"
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            build(Snapshot(root / "data"), tmp)
            output = Path(tmp)
            report = (output / "audit_report.md").read_text(encoding="utf-8")
            page = (output / "index.html").read_text(encoding="utf-8")
            evidence = json.loads((output / "evidence.json").read_text(encoding="utf-8"))
            with (output / "candidate_evidence.csv").open(newline="", encoding="utf-8") as handle:
                csv_rows = list(csv.DictReader(handle))

        self.assertEqual([locus["qualityControls"] for locus in evidence["loci"]], [[expected]] * 5)
        self.assertEqual(report.count(f"| {expected} |"), 5)
        self.assertIn("provides no additional detail", report)
        self.assertIn("<strong>Source quality controls:</strong>", page)
        self.assertIn("provides no additional detail", page)
        self.assertEqual({row["source_quality_controls"] for row in csv_rows}, {expected})


if __name__ == "__main__":
    unittest.main()
