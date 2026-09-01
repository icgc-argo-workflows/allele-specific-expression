import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEMVER_PATTERN = re.compile(
    r"^(0|[1-9]\d*)\."
    r"(0|[1-9]\d*)\."
    r"(0|[1-9]\d*)"
    r"(?:-[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?"
    r"(?:\+[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?$"
)


class VersionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
        cls.nextflow_config = (ROOT / "nextflow.config").read_text(
            encoding="utf-8"
        )

    def test_version_uses_semantic_versioning(self):
        self.assertRegex(self.version, SEMVER_PATTERN)

    def test_nextflow_manifest_matches_version_file(self):
        manifest = re.search(
            r"manifest\s*\{(?P<body>.*?)\}", self.nextflow_config, re.DOTALL
        )
        self.assertIsNotNone(manifest, "nextflow.config has no manifest block")

        manifest_version = re.search(
            r"^\s*version\s*=\s*['\"](?P<version>[^'\"]+)['\"]",
            manifest.group("body"),
            re.MULTILINE,
        )
        self.assertIsNotNone(
            manifest_version, "Nextflow manifest has no literal version"
        )
        self.assertEqual(manifest_version.group("version"), self.version)

    def test_changelog_contains_current_version(self):
        changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
        self.assertIn(f"## [{self.version}]", changelog)

    def test_readme_badge_matches_current_version(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            f"img.shields.io/badge/version-{self.version}-blue.svg", readme
        )

    def test_readme_has_license_badge(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "img.shields.io/github/license/"
            "icgc-argo-workflows/allele-specific-expression",
            readme,
        )

    def test_citation_matches_current_version(self):
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertRegex(
            citation,
            rf"(?m)^version:\s*['\"]?{re.escape(self.version)}['\"]?\s*$",
        )

    def test_license_is_populated(self):
        license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
        self.assertIn("MIT License", license_text)
        self.assertIn("Permission is hereby granted", license_text)

    def test_conda_lock_tracks_environment_file(self):
        lock_text = (ROOT / "conda-lock.yml").read_text(encoding="utf-8")
        self.assertIn("linux-64:", lock_text)
        self.assertIn("- ase-environment.yaml", lock_text)


if __name__ == "__main__":
    unittest.main()
