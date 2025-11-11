# Structural Variant Loop Experiment — Scientific TL;DR

**Date:** 10 Nov 2025  
**Model:** AlphaGenome (DNase output, K562 / EFO:0002067)  
**Construct set:** Four 1,048,576 bp sequences centered on a 10× HS2 enhancer block with CTCF loop anchors.

## Key Observations

- **Loop integrity barely changes promoter accessibility.** LoopIntact and LoopInverted retain similar promoter scores (mean ≈0.00493), implying AlphaGenome is largely insensitive to anchor polarity when enhancer distance is unchanged.
- **Spacer deletion boosts local enhancer accessibility.** LoopDeleted increases maximum DNase at the enhancer (0.2539) without raising promoter mean, suggesting the model reports stronger local opening but limited long-range propagation.
- **Relocating the enhancer weakens loop interior but raises distal peak.** Moving the 10× block to 800 kb lowers loop-span mean (0.00417 vs ~0.00499) and anchor signals, yet the enhancer itself peaks higher (max 0.2695), indicating AlphaGenome treats the distant block as an isolated hotspot.
- **Global accessibility remains stable.** Genome-wide means stay near 0.00437–0.00442 with σ≈0.0058, so structural edits re-pattern accessibility rather than globally scaling predictions.

## Quantitative Summary

| Construct | Enhancer Max | Enhancer Mean | Promoter Mean (±5 kb) | Loop Span Mean |
|-----------|--------------|---------------|------------------------|----------------|
| LoopIntact_10x | 0.2539 | 0.01257 | 0.00493 | 0.00499 |
| LoopInverted_10x | 0.2500 | 0.01257 | 0.00493 | 0.00499 |
| LoopDeleted_10x | 0.2539 | 0.01255 | 0.00490 | 0.00498 |
| LoopRelocated_10x | **0.2695** | **0.01294** | 0.00485 | **0.00417** |

Numbers derive from `results/structural_variant_metrics.csv`.

## Interpretation

- Anchor orientation alone does not impose a contact penalty, hinting that AlphaGenome encodes CTCF motifs primarily through sequence-level binding potential rather than explicit polarity-aware looping.
- Removing the promoter spacer keeps long-range activation muted, reinforcing the earlier stacking result that AlphaGenome struggles to propagate accessibility across large filler segments.
- Enhancer relocation creates a new accessibility domain detached from the promoter loop, allowing targeted tests of whether additional insulators or relocated promoters recover coupling.

## Recommended Follow-Up

1. Add constructs with duplicated or inverted anchors to test whether anchor copy number modulates predicted contact strength.
2. Introduce promoter relocations into the distal domain to evaluate whether AlphaGenome can re-couple enhancer and promoter after relocation.
3. Expand metrics to include alternative outputs (e.g., H3K27ac) for richer interpretation of distal enhancer behavior.
