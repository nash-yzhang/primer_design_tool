# Scope

A interest driven spontaneous project for browser-based PCR primer design tool. Designed by me ([nash-yzhang @ github](https://github.com/nash-yzhang)) and implemented with Claude's aid. Identifies conserved primer binding sites in a multiple sequence alignment (MSA) using windowed Jensen-Shannon Divergence (JSD) scoring ([Capra & Singh, Bioinformatics 2009](https://academic.oup.com/bioinformatics/article/23/15/1875/203579)). Open-source under MIT License, with no dependencies or backend - runs entirely in the browser.

## Quick Start

1. Open primer_design_tool_standalone.html in any modern browser
2. Drop or click to load a BLAST-aligned multi-FASTA file (.fasta, .fa, .txt, .aln)
3. Adjust parameters if needed, then click **Run Analysis**
4. Click any candidate primer to inspect per-locus conservation, GC%, Tm, and haplotype variants

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **k** - score threshold | 0.60 | Minimum 3-end-weighted JSD score to retain a candidate |
| **Min / Max length** | 18 / 24 nt | Primer length search range |
| **c** - locus threshold | 0.70 | Per-position JSD cutoff for conserved label |

## Algorithm

For each alignment column, computes JSD between the observed base frequency and a uniform background (lambda = 0.5). Applies a +-3 nt sliding-window smooth, then scores each candidate window with a 3-end-weighted sum. Greedy overlap removal retains the highest-scoring non-redundant set.

## Project layout

`
primer_design_tool_standalone.html  <- single-file app, open directly
primerscope/
  tool.html          <- componentised entry point
  primerscope.css    <- all styles
  core.js            <- JSD algorithm, FASTA parser
  ui.js              <- rendering, events, state
  config.json        <- default parameter reference
  example.fasta      <- test alignment (20 seqs x 120 nt)
`

primerscope/ is identical in behaviour - HTML only loads the two JS files and the CSS.

## Input format

Standard multi-FASTA, pre-aligned to equal length.

## License

MIT - see [LICENSE](./LICENSE)
