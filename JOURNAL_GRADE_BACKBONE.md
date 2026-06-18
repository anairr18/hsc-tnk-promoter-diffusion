# Computational Backbone Upgrades

This project now has the components needed for a stronger computational-journal proof of concept:

- explicit random, chromosome-heldout, and gene-heldout split bundles
- exact and approximate leakage/novelty auditing
- pretrained fine-tuning patch support for DNA-Diffusion
- continuous promoter-profile conditioning patch support
- ensemble predictor input/merge interface for Enformer, Borzoi, ChromBPNet-style scores
- explicit guided reward ranking for target-high/offtarget-low candidates
- active-learning round design from MPRA measurements
- MPRA library and barcode-count analysis scripts

The strongest POC should report results across at least:

1. from-scratch DNA-Diffusion baseline
2. pretrained fine-tuned DNA-Diffusion
3. random split
4. chromosome-heldout split
5. gene-heldout split
6. endogenous and pretrained generated controls
7. leakage/novelty audit
8. predictor-guided ranking with target/off-target reward

This is the computational evidence package most likely to survive reviewer questions before wet-lab validation.
