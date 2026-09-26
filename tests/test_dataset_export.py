import csv
import hashlib
from pathlib import Path
import stat
import sys
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import export_real_dataset as dataset_export


def csv_text(value):
    """Match csv.DictWriter's representation of scalar source values."""
    return "" if value is None else str(value)


class RealDatasetExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence, cls.evidence_hash = dataset_export.load_evidence(
            ROOT / "reports" / "evidence.json"
        )
        cls.dataset_dir = ROOT / "reports" / "dataset"
        cls.rows = {}
        for name in (
            "variants",
            "l2g_candidates",
            "molecular_qtl_colocalisations",
        ):
            with (cls.dataset_dir / f"{name}.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                cls.rows[name] = list(csv.DictReader(handle))

    def test_pinned_source_and_exact_row_counts(self):
        self.assertEqual(self.evidence_hash, dataset_export.PINNED_EVIDENCE_SHA256)
        self.assertEqual(
            {name: len(rows) for name, rows in self.rows.items()},
            dataset_export.EXPECTED_COUNTS,
        )

    def test_source_order_and_primary_identifiers_are_preserved(self):
        expected_variants = [
            source["variant"]["id"]
            for locus in self.evidence["loci"]
            for source in locus["variants"]["rows"]
        ]
        expected_candidates = [
            source["target"]["id"]
            for locus in self.evidence["loci"]
            for source in locus["candidates"]["rows"]
        ]
        expected_colocalisations = [
            source["otherStudyLocus"]["studyLocusId"]
            for locus in self.evidence["loci"]
            for source in locus["colocalisations"]["rows"]
        ]
        self.assertEqual(
            [row["variant_id"] for row in self.rows["variants"]], expected_variants
        )
        self.assertEqual(
            [row["target_gene_id"] for row in self.rows["l2g_candidates"]],
            expected_candidates,
        )
        self.assertEqual(
            [row["other_study_locus_id"] for row in self.rows["molecular_qtl_colocalisations"]],
            expected_colocalisations,
        )

    def test_variant_source_metrics_are_not_changed(self):
        exported = iter(self.rows["variants"])
        for locus in self.evidence["loci"]:
            for source in locus["variants"]["rows"]:
                row = next(exported)
                for output_name, source_name in (
                    ("posterior_probability", "posteriorProbability"),
                    ("is_95_credible_set", "is95CredibleSet"),
                    ("is_99_credible_set", "is99CredibleSet"),
                    ("beta", "beta"),
                    ("standard_error", "standardError"),
                ):
                    self.assertEqual(row[output_name], csv_text(source.get(source_name)))

    def test_l2g_scores_and_every_feature_metric_are_not_changed(self):
        exported = iter(self.rows["l2g_candidates"])
        for locus in self.evidence["loci"]:
            for source in locus["candidates"]["rows"]:
                row = next(exported)
                self.assertEqual(row["l2g_score"], csv_text(source.get("score")))
                features = {feature["name"]: feature for feature in source["features"]}
                self.assertEqual(set(features), set(dataset_export.L2G_FEATURE_NAMES))
                for name, feature in features.items():
                    self.assertEqual(
                        row[f"feature_{name}_value"], csv_text(feature.get("value"))
                    )
                    self.assertEqual(
                        row[f"feature_{name}_shap_value"],
                        csv_text(feature.get("shapValue")),
                    )

    def test_colocalisation_source_metrics_are_not_changed(self):
        exported = iter(self.rows["molecular_qtl_colocalisations"])
        for locus in self.evidence["loci"]:
            for source in locus["colocalisations"]["rows"]:
                row = next(exported)
                for output_name, source_name in (
                    ("h3", "h3"),
                    ("h4", "h4"),
                    ("clpp", "clpp"),
                    ("number_colocalising_variants", "numberColocalisingVariants"),
                    ("beta_ratio_sign_average", "betaRatioSignAverage"),
                ):
                    self.assertEqual(row[output_name], csv_text(source.get(source_name)))

    def test_missing_values_remain_blank_and_zero_remains_zero(self):
        candidates = self.rows["l2g_candidates"]
        colocalisations = self.rows["molecular_qtl_colocalisations"]
        self.assertEqual(sum(row["candidate_max_h4_fetched"] == "" for row in candidates), 9)
        self.assertEqual(sum(row["qtl_gene_id"] == "" for row in colocalisations), 44)
        self.assertTrue(any(row["h3"] == "0.0" for row in colocalisations))
        tcf7l2 = next(row for row in candidates if row["target_gene_symbol"] == "TCF7L2")
        self.assertEqual(tcf7l2["candidate_max_h4_fetched"], "")
        self.assertEqual(tcf7l2["feature_eQtlColocH4Maximum_value"], "0")

    def test_provenance_qc_and_completeness_are_present_on_every_row(self):
        all_rows = [row for rows in self.rows.values() for row in rows]
        self.assertEqual(
            {row["record_status"] for row in all_rows},
            {"public_source_derived_pending_human_review"},
        )
        self.assertEqual(
            {row["source_quality_controls"] for row in all_rows},
            {"Study has quality control flag(s)"},
        )
        self.assertEqual(
            {row["source_evidence_sha256"] for row in all_rows},
            {dataset_export.PINNED_EVIDENCE_SHA256},
        )
        self.assertEqual({row["source_collection_complete"] for row in all_rows}, {"True"})
        self.assertTrue(all(row["source_url"].startswith("https://platform.opentargets.org/") for row in all_rows))

    def test_schema_contains_no_new_combined_or_causal_score(self):
        fields = (
            set(dataset_export.VARIANT_FIELDS)
            | set(dataset_export.L2G_FIELDS)
            | set(dataset_export.COLOCALISATION_FIELDS)
        )
        forbidden = {
            "combined_score",
            "composite_score",
            "causal_probability",
            "core_probability",
            "actionability_probability",
        }
        self.assertTrue(fields.isdisjoint(forbidden))
        self.assertEqual(len(dataset_export.L2G_FEATURE_FIELDS), 62)

    def test_readme_states_real_data_boundaries_and_defines_terms(self):
        readme = (self.dataset_dir / "README.md").read_text(encoding="utf-8")
        for required in (
            "real public aggregate records",
            "not participant-level records",
            "not a complete GWAS summary-statistics download",
            "65 variant rows means 65 variant records, not 65 people",
            "Missing source values are empty CSV cells",
            "## Key terms in plain language",
            "**GWAS:**",
            "**L2G:**",
            "**Colocalisation:**",
            "raw five-locus Open Targets API snapshot used by this project",
            "not full FinnGen GWAS summary statistics",
            "no participant-level data",
            "nine unmodified JSON response bodies",
        ):
            self.assertIn(required.casefold(), readme.casefold())

    def test_raw_archive_has_exact_unmodified_members_and_fixed_metadata(self):
        archive_path = self.dataset_dir / dataset_export.RAW_ARCHIVE_NAME
        source_members = dataset_export.raw_source_members(ROOT / "data")
        expected_names = [name for name, _ in source_members]
        with zipfile.ZipFile(archive_path) as archive:
            self.assertIsNone(archive.testzip())
            self.assertEqual(archive.namelist(), expected_names)
            self.assertEqual(len(archive.infolist()), 10)
            for info, (expected_name, source_path) in zip(
                archive.infolist(), source_members, strict=True
            ):
                self.assertEqual(info.filename, expected_name)
                self.assertEqual(info.date_time, dataset_export.RAW_ARCHIVE_TIMESTAMP)
                self.assertEqual(info.create_system, 3)
                self.assertEqual(info.compress_type, zipfile.ZIP_STORED)
                self.assertEqual(info.external_attr >> 16, stat.S_IFREG | 0o644)
                self.assertEqual(info.extra, b"")
                self.assertEqual(info.comment, b"")
                self.assertEqual(archive.read(info), source_path.read_bytes())

        archive_bytes = archive_path.read_bytes()
        self.assertEqual(len(archive_bytes), 306502)
        self.assertEqual(
            hashlib.sha256(archive_bytes).hexdigest(),
            "9418e42f6c82d04976b310a4fd558df3d5da6a2bfcc0f3123c339f567c0ad994",
        )

    def test_regeneration_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            dataset_export.export_dataset(ROOT / "reports" / "evidence.json", first)
            dataset_export.export_dataset(ROOT / "reports" / "evidence.json", second)
            for name in dataset_export.OUTPUT_FILES:
                first_bytes = (Path(first) / name).read_bytes()
                self.assertEqual(first_bytes, (Path(second) / name).read_bytes(), name)
                self.assertEqual(first_bytes, (self.dataset_dir / name).read_bytes(), name)
                if name != dataset_export.RAW_ARCHIVE_NAME:
                    self.assertNotIn(b"\r\n", first_bytes, name)


if __name__ == "__main__":
    unittest.main()
