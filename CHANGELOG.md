# Changelog

All notable changes to this project are documented in this file. The project
uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Input preflight validation for alignment, VCF, and reference contig compatibility and VCF reference alleles.
- A lightweight end-to-end Nextflow smoke test in CI.
- Automatic conversion and coordinate sorting of SAM inputs.
- Phase-set-aware haplotype aggregation, VCF/BAM sample identity checks, and explicit biallelic SNV filtering.
- Run provenance outputs, reproducible Linux dependency locking, citation metadata, and a populated MIT license.
- Release container publication to GitHub Container Registry.
- Automated Markdown linting for the README and changelog.

### Fixed

- Updated the workflow to pass Nextflow strict-syntax linting.
- Consolidated smoke and integration data under `tests/data`, corrected the optional integration settings, and fixed documentation and log-message typos.
- Hardened reference downloads with retries, resume support, and SHA-256 integrity checks.

## [1.0.0] - 2026-08-31

### Added

- Formal semantic version metadata for the pipeline.
- Automated Python, Nextflow, and container checks with GitHub Actions.
- Tag-driven GitHub releases with version consistency checks.

### Fixed

- Deterministic PyEnsembl database publication from `createGtfDB`.
- Compatibility with current scientific Python and bioinformatics packages.
- VCF sample selection, contig matching, and several ASE analysis edge cases.

[Unreleased]: https://github.com/schwarzlab-ccb/pipeline-ase/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/schwarzlab-ccb/pipeline-ase/releases/tag/v1.0.0
