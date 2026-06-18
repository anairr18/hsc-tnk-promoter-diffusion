# DNA-Diffusion HSC Project — Export for Codex

START HERE: read `PROJECT_CONTEXT.md` first. It contains the full project
history, decisions, current state, gotchas, and next steps so you (or an
agent) can resume work without re-deriving anything.

## Folder structure

```
PROJECT_CONTEXT.md       <- read this first
pipeline/                <- current pipeline scripts (dataset + training + generation)
scripts/                 <- validation, HSC/T-NK promoter, ranking, and MPRA scripts
reports/                 <- generated POC reports after validation
reference_data/          <- files Giacomo sent (pretrained-model generations + endogenous data)
```

## How to use this with Codex

1. Open a new Codex session (CLI or IDE extension).
2. Point it at this folder, or paste `PROJECT_CONTEXT.md` into the chat as
   the first message so it has full context before doing anything.
3. Ask Codex to adapt `pipeline/cell1_dataset_and_training.py` for your
   actual compute environment (remove the `google.colab.drive` mount call,
   replace the cache path with wherever you want data stored — local disk,
   S3, a persistent volume, etc.).
4. Have Codex fix the output-buffering issue from Colab (unbuffered
   subprocess output / streamed logs) so training progress is visible.
5. Proceed through the "Immediate Next Steps" section in
   `PROJECT_CONTEXT.md`.

## Implemented plan

See `IMPLEMENTATION_STATUS.md` for the concrete two-stage execution flow:

1. Finish the K562/HepG2/GM12878 cell-line proof of concept.
2. Build the HSC-to-T/NK CAR promoter generation and MPRA validation pipeline.
