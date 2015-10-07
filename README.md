superFreq
==========

A CNA, SNV and clonality tracker, taking full advantage of multiple samples from the same individual.

The main selling points are
- Extensive output, both plots and tables, diagnostic and analytical. Typically both the raw data and identified interesting features are shown, making it easier for the user to spot patterns and identify issues with the data or analysis.
- SuperFreq uses a pool of (potentially unrelated) normal samples to filter SNVs and estimate variance in coverage for CNA calls.
- All variants (SNVs and CNAs) are cross-checked between samples of the same individual, allowing detection of very low clonality mutations if the mutation is present at higher clonality in a related sample.
- The CNA calls are allele sensitive and use both coverage and SNP frequency to segment the genome. Does not require a matched normal sample, but will use a matched normal if available.
- The clonality tracking uses both SNVs and CNAs to group mutations into clones.
- The phylogenetic tree is checked for self-consistency.

Best way to run it
==================
The instructions are in the manual.

Basically it involves downloading an example data set, run the pipeline on the example, and if that works, point the pipeline to your data and rerun.
