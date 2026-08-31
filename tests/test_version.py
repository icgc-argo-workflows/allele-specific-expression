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


if __name__ == "__main__":
    unittest.main()
