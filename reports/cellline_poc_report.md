# Cell-Line Proof-of-Concept Report

## Completion State

- Fine-tuned generated directory: not found
- POC generation complete: no
- Enformer validation: not run by this script; run the GPU Enformer pipeline after generation.

## Sequence QC Summary

| group | cell_type | n | unique | duplicate_rate | gc_mean | cpg_mean | novelty_vs_endogenous | novelty_vs_pretrained |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| fine_tuned | K562 | 0 | 0 | nan | nan | nan | nan | nan |
| fine_tuned | HepG2 | 0 | 0 | nan | nan | nan | nan | nan |
| fine_tuned | GM12878 | 0 | 0 | nan | nan | nan | nan | nan |
| pretrained | K562 | 1000 | 1000 | 0.000 | 0.499 | 4.196 | 1.000 | 0.000 |
| pretrained | HepG2 | 1000 | 1000 | 0.000 | 0.464 | 2.720 | 1.000 | 0.000 |
| pretrained | GM12878 | 1000 | 1000 | 0.000 | 0.432 | 2.406 | 1.000 | 0.000 |
| endogenous | K562 | 11968 | 11968 | 0.000 | 0.496 | 3.329 | 0.000 | 1.000 |
| endogenous | HepG2 | 11968 | 11968 | 0.000 | 0.462 | 2.518 | 0.000 | 1.000 |
| endogenous | GM12878 | 11968 | 11968 | 0.000 | 0.453 | 2.542 | 0.000 | 1.000 |

## Go/No-Go Notes

- Go to hematopoietic implementation only after fine-tuned sequences exist for all three cell types.
- Keep pretrained and endogenous references as controls for all downstream comparisons.
- If fine-tuned sequences are not novel or drift strongly from endogenous GC/k-mer distributions, revisit training length and pretrained initialization.
