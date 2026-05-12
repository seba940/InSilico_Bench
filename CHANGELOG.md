# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-05-12
### Added
- **Multi-species Gene Lookup**: Support for multiple species beyond SGD (Saccharomyces cerevisiae), including Human, Mouse, PomBase, FlyBase, and WormBase.
- **Provider Registry Pattern**: Architectural foundation for easily adding new genomic databases.
- **Annotation Metadata**: Added `species` and `provider` fields to Annotation References.
- **Dynamic UI**: Gene Lookup UI now adapts its database buttons and links based on the selected species.
- **"No GFF" Mode**: Directly open external databases using gene names even without a local GFF file.

## [1.0.0] - 2026-05-12
### Added
- **SnapGene-style Gene Viewer**: A unified sequence visualization module across all tabs in InSilico_Bench.
- **PCR Filtering**: Only selected primers are visualized in the PCR execution tab for better clarity.
- **Interactive Features**: Clicking on a primer or feature now automatically selects the sequence range and displays Tm/GC% information.
- **Adjustable Layout**: Added a PanedWindow in the Primer Tab to allow resizing between Forward and Reverse tables.
- **Search & Find**: Added Ctrl+F sequence search functionality for templates.
- **Project Documentation**: Created a comprehensive README.md and this Changelog.
- **GitHub Integration**: Initialized repository and pushed to GitHub.

### Fixed
- Resolved `NameError` in multiple tabs regarding `SnapGeneViewer` imports.
- Fixed `AttributeError` in `find_primer_bindings` when handling filtered primer lists.
- Corrected table layout alignment issues in the Primer tab.
