# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [2025-04-24]

Initial release of nf-core/abotyper, created with the [nf-core](https://nf-co.re/) template.

### Added

- Initial ABO\*A1/A2 differentiation algorithm (WORK IN PROGRESS)
- Metadata propagation for exonic information
- Improved meta.yaml documentation for local modules
- Updated documenting and staging of versions.yaml files for MultiQC
- Testing docker and singularity containers for requisite profiles support

### Changed

- Removed modules aliasing for cleaner code structure
- Improved staging for MultiQC files
- Updated MultiQC custom config to improve reporting
- Updated modules config file to match removed modules aliasing
- Moved index preparation to main workflow (streamlining processes)
- Upgraded all python scripts to class-style declaration
- Refactored all Python scripts to remove stale imports and improve maintainability

### Fixed

- Fixed MultiQC configuration errors
- Improved file staging for MultiQC compatibility
- Fixed issues with metadata propagation in the pipeline to reduce modules aliasing

### Known Issues

- MultiQC-1.30 treats the `\.` at the end of the `multiqc module command` as an illegal character causing failure. After long discussions and testing on different platforms, it was determined that this is a potential `AVX2` issue. Some python libraries may have been compiled with AVX2 instructions under containerised environments.

## [1.1.0] - [2025-07-31]

### New Features

- Initial release of ABO typing pipeline
- Comprehensive test suite with `test` profile support
- Debug mode compatibility for enhanced troubleshooting
- Updated usage documentation in `docs/usage.md`
- Updated output documentation in `docs/output.md`
- Tool citations and contributor acknowledgments in `README.md`

### Improvements

- Pipeline now fully compatible with nf-core standards
- Enhanced linting compliance (`nf-core pipelines lint`)
- Improved test dataset integration

### Bug Fixes

- All pipeline tests now pass successfully
- Debug mode warnings resolved
- Linting issues addressed

### Documentation Updates

- Updated `README.md` with new tool citations and contributors
- Comprehensive `docs/usage.md` documentation
- Detailed `docs/output.md` with pipeline outputs
- This `CHANGELOG.md` updated with all changes
