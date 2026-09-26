import copy
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import build_case_study as case_study


class RealDataCaseStudyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence_path = ROOT / "reports" / "evidence.json"
        cls.evidence, cls.input_hash = case_study.load_evidence(cls.evidence_path)
        cls.locus_rows, cls.candidate_rows, cls.counts = case_study.summarize_case(cls.evidence)

    def test_pinned_snapshot_and_counts(self):
        self.assertEqual(
            self.input_hash,
            "49e73e407ba8a062c7b0a32f145233d24505dd381e675c550e5ba76267f5c821",
        )
        self.assertEqual(self.evidence["study"]["id"], "FINNGEN_R12_T2D")
        self.assertEqual(
            self.counts,
            {
                "loci": 5,
                "variant_rows": 65,
                "l2g_candidates": 11,
                "colocalisation_rows": 44,
                "candidate_matched_colocalisation_rows": 13,
                "unmatched_colocalisation_rows": 31,
                "candidate_long_rows": 53,
            },
        )

    def test_candidate_matching_uses_either_source_gene_field(self):
        target = "ENSG_TEST"
        by_qtl_gene = {"otherStudyLocus": {"qtlGeneId": target, "study": {"target": {"id": "other"}}}}
        by_study_target = {"otherStudyLocus": {"qtlGeneId": None, "study": {"target": {"id": target}}}}
        unrelated = {"otherStudyLocus": {"qtlGeneId": "other", "study": {"target": {"id": "also-other"}}}}
        self.assertTrue(case_study.colocalisation_matches_candidate(by_qtl_gene, target))
        self.assertTrue(case_study.colocalisation_matches_candidate(by_study_target, target))
        self.assertFalse(case_study.colocalisation_matches_candidate(unrelated, target))

    def test_missing_colocalisation_metrics_remain_missing(self):
        tcf7l2 = next(row for row in self.candidate_rows if row["gene_symbol"] == "TCF7L2")
        self.assertEqual(tcf7l2["candidate_coloc_status"], "not_returned")
        self.assertEqual(tcf7l2["matched_colocalisation_rows"], 0)
        self.assertIsNone(tcf7l2["h4"])
        self.assertIsNone(tcf7l2["clpp"])

    def test_l2g_scores_are_preserved_without_normalisation(self):
        ccnd2_locus = next(
            locus for locus in self.evidence["loci"] if locus["variant"]["id"] == "12_4275678_T_G"
        )
        source_scores = {row["target"]["approvedSymbol"]: row["score"] for row in ccnd2_locus["candidates"]["rows"]}
        output_scores = {
            row["gene_symbol"]: row["l2g_score"]
            for row in self.candidate_rows
            if row["lead_variant"] == "12_4275678_T_G" and row["gene_symbol"] is not None
        }
        self.assertEqual(output_scores, source_scores)
        self.assertNotAlmostEqual(sum(output_scores.values()), 1.0)

    def test_h4_and_clpp_remain_distinct_source_fields(self):
        cdkal1 = [row for row in self.candidate_rows if row["gene_symbol"] == "CDKAL1"]
        self.assertEqual(len(cdkal1), 2)
        single_cell = next(row for row in cdkal1 if row["qtl_type"] == "sceqtl")
        self.assertAlmostEqual(single_cell["h4"], 0.9855721010213954)
        self.assertAlmostEqual(single_cell["clpp"], 0.06108280581803706)
        self.assertNotEqual(single_cell["h4"], single_cell["clpp"])

        locus = next(row for row in self.locus_rows if row["top_l2g_gene"] == "CDKAL1")
        self.assertAlmostEqual(locus["maximum_returned_h4"], 0.9925491069838929)
        self.assertAlmostEqual(locus["candidate_matched_maximum_h4"], 0.9855721010213954)
        report = case_study.markdown_report(self.evidence, self.locus_rows, self.counts)
        self.assertIn("Locus-wide max H4", report)
        self.assertIn("Returned-candidate max H4", report)

    def test_fto_locus_keeps_unmatched_qtl_rows_visible(self):
        fto = next(row for row in self.locus_rows if row["top_l2g_gene"] == "FTO")
        self.assertEqual(fto["variant_rows"], 56)
        self.assertAlmostEqual(fto["fetched_pip_sum"], 0.9505637925173493)
        self.assertAlmostEqual(fto["maximum_returned_variant_pip"], 0.0459602012564266)
        self.assertEqual(fto["total_colocalisation_rows"], 2)
        self.assertEqual(fto["candidate_matched_colocalisation_rows"], 0)
        self.assertEqual(fto["unmatched_colocalisation_rows"], 2)

    def test_long_table_retains_every_source_colocalisation(self):
        coloc_records = [
            row
            for row in self.candidate_rows
            if row["record_type"] in {"candidate_colocalisation_match", "unmatched_locus_colocalisation"}
        ]
        self.assertEqual(len(coloc_records), 44)
        self.assertEqual(
            sum(row["record_type"] == "candidate_colocalisation_match" for row in coloc_records),
            13,
        )
        self.assertEqual(
            sum(row["record_type"] == "unmatched_locus_colocalisation" for row in coloc_records),
            31,
        )

    def test_every_output_row_retains_status_and_qc_warning(self):
        all_rows = self.locus_rows + self.candidate_rows
        self.assertTrue(all(row["status"] == case_study.STATUS for row in all_rows))
        self.assertEqual(
            {row["source_quality_controls"] for row in all_rows},
            {"Study has quality control flag(s)"},
        )

    def test_output_schema_has_no_composite_or_new_posterior(self):
        fields = set(case_study.LOCUS_FIELDS) | set(case_study.CANDIDATE_FIELDS)
        forbidden = {
            "composite_score",
            "combined_score",
            "causal_probability",
            "core_probability",
            "actionability_probability",
        }
        self.assertTrue(fields.isdisjoint(forbidden))

    def test_probability_validation_distinguishes_zero_and_missing(self):
        self.assertIsNone(case_study.probability(None, "metric"))
        self.assertEqual(case_study.probability(0.0, "metric"), 0.0)
        for invalid in (True, -0.01, 1.01, float("nan"), float("inf"), "0.5"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                case_study.probability(invalid, "metric")

    def test_validation_rejects_schema_and_count_corruption(self):
        bad_schema = copy.deepcopy(self.evidence)
        bad_schema["schema_version"] = 2
        with self.assertRaisesRegex(ValueError, "schema_version 1"):
            case_study.validate_evidence(bad_schema)

        bad_count = copy.deepcopy(self.evidence)
        bad_count["loci"][0]["variants"]["fetched_count"] += 1
        with self.assertRaisesRegex(ValueError, "fetched_count"):
            case_study.validate_evidence(bad_count)

    def test_regeneration_is_byte_deterministic_and_does_not_mutate_input(self):
        before = case_study.sha256_bytes(self.evidence_path.read_bytes())
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            case_study.build_case_study(self.evidence_path, first)
            case_study.build_case_study(self.evidence_path, second)
            for name in case_study.OUTPUT_FILES:
                first_bytes = (Path(first) / name).read_bytes()
                second_bytes = (Path(second) / name).read_bytes()
                committed_bytes = (ROOT / "reports" / name).read_bytes()
                self.assertEqual(first_bytes, second_bytes, name)
                self.assertEqual(first_bytes, committed_bytes, name)
                self.assertNotIn(b"\r\n", first_bytes, name)
        self.assertEqual(before, case_study.sha256_bytes(self.evidence_path.read_bytes()))

    def test_csv_blank_cells_represent_missing_colocalisation(self):
        with (ROOT / "reports" / "candidate_evidence_long.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        tcf7l2 = next(row for row in rows if row["gene_symbol"] == "TCF7L2")
        self.assertEqual(tcf7l2["h4"], "")
        self.assertEqual(tcf7l2["clpp"], "")
        self.assertNotEqual(tcf7l2["h4"], "0")

    def test_metadata_declares_descriptive_boundaries(self):
        metadata = json.loads((ROOT / "reports" / "case_study_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["analysis_type"], "descriptive_real_data_case_study")
        self.assertFalse(metadata["principles"]["composite_score_created"])
        self.assertFalse(metadata["principles"]["posterior_created"])
        self.assertFalse(metadata["principles"]["missing_values_converted_to_zero"])


if __name__ == "__main__":
    unittest.main()
