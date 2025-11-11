# AlphaGenome Structural Variant Loop Experiment

This experiment extends the enhancer stacking study by introducing simple
structural variant operations around the 10× HS2 enhancer array. We probe
whether AlphaGenome responds to topological rewiring or only to sequence
composition by adding convergent CTCF loop anchors, inverting anchors,
deleting the enhancer–promoter spacer, and relocating the enhancer block to a
distal site.

## TL;DR

Quick results and interpretation live in `summary.md`. Highlights:

- Promoter accessibility remains ~0.0049 regardless of loop polarity, implying AlphaGenome is largely orientation agnostic.
- Removing the spacer sharpens enhancer accessibility locally but does not lift promoter signal.
- Relocating the enhancer to 800 kb raises the enhancer peak (0.2695) while collapsing signal within the anchored loop, indicating the model treats the distal block as an isolated hotspot.

## Constructs Produced

Four constructs are generated automatically by
`build_structural_variant_constructs.py`:

| Construct | Description |
|-----------|-------------|
| `LoopIntact_10x` | 10× HS2 array bracketed by convergent CTCF anchors at 350 kb and 450 kb; promoter at 500 kb. |
| `LoopInverted_10x` | Same as above but the right CTCF anchor is flipped (non-convergent orientation). |
| `LoopDeleted_10x` | Convergent anchors; spacer between right anchor and promoter is deleted so the promoter abuts the loop. |
| `LoopRelocated_10x` | Convergent anchors remain, but the 10× HS2 array is moved to 800 kb, outside the anchored loop. |

All constructs remain 1,048,576 bp long so they satisfy AlphaGenome's input
length requirement.

## Pipeline

1. **Build constructs**: `python build_structural_variant_constructs.py`
   - Creates FASTA files under `sequences/`
   - Writes a manifest with coordinates and operations
2. **Run AlphaGenome**: `python run_structural_variant_predictions.py`
   - Requires `ALPHA_GENOME_KEY` (or `ALPHA_GENOME_API_KEY`) in the environment
   - Saves predictions to `alphagenome_outputs/` and execution logs to `logs/`

Optional: the metadata manifest can be imported into notebooks or downstream
analysis scripts to align prediction outputs with the structural edits that
were applied.

## Requirements

- Same dependencies as the root enhancer stacking project
- AlphaGenome API access configured via `.env` or exported environment
  variable

## Next Steps

- Expand the variant library with additional operations (duplications,
  inversions of the enhancer block, insulation with CTCF boundaries, etc.)
- Add automated plotting/analysis similar to the main experiment once
  predictions are available.
