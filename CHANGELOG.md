# Changelog

## [Unreleased]

## [1.0.1] - 2026-08-15

- Cached codon scores and finite-state constraint transitions within each exact optimization while preserving optimality and deterministic tie-breaking.
- Added a deterministic end-to-end constrained-optimization benchmark with exact result checksums.


## [1.0.0] - 2026-08-12

First stable release.

- Exact finite-state synonymous codon optimization under DNA constraints.
- Added a deterministic property-based suite for translation, optimality, and reverse-complement invariants.
- Documented mutation testing: 768 of 807 mutants killed (95.17%); all 39 survivors reviewed as behavior-equivalent.
- Adopted strict mypy checking and shipped typed-package metadata.
- Expanded CI across Linux and macOS on Python 3.11–3.13.
