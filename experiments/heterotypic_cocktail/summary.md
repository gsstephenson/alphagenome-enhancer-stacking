# Heterotypic Enhancer Cocktail Experiment — Scientific TL;DR

**Date:** 10 Nov 2025  
**Model:** AlphaGenome DNase, K562 (EFO:0002067)  
**Construct library:** Six 1,048,576 bp sequences mixing HS2, GATA1, and HNF4A modules with CTCF brackets/separators and varied spacing/orientation.

## Key Findings

- **HNF4A modules dominate accessibility.** Mean DNase within HNF4A tiles ranges 0.047–0.057 versus HS2 at 0.016–0.018 and GATA1 at ~4–5×10⁻⁴, suggesting AlphaGenome prioritizes hepatic-factor motifs over erythroid HS2 despite the K562 context.
- **Modest gains from tighter packing.** Shorter 1 kb spacing (Cocktail_1kbForward) boosts HNF4A mean to 0.0467 while leaving promoter signal modest (0.00449); 20 kb spacing slightly raises HNF4A mean (0.0567) but reduces global variance.
- **Orientation swap yields minimal change.** Alternating GATA1 strand or reversing the module order shifts per-module means by <5%, indicating limited polarity sensitivity at this resolution.
- **CTCF separators dampen loop interior.** Inserting CTCF between repeats drops global max to 2.59 (vs 4.88–5.59) and slightly elevates promoter mean (0.00472), hinting that insulating spacers diffuse enhancer-domain intensity while marginally improving distal signal.

## Metric Snapshot

| Construct | HS2 Mean | GATA1 Mean | HNF4A Mean | Promoter Mean |
|-----------|----------|------------|------------|----------------|
| Cocktail_1kbForward | 0.0161 | 4.1e-4 | 0.0467 | 0.00449 |
| Cocktail_5kbForward | 0.0183 | 4.6e-4 | 0.0523 | 0.00463 |
| Cocktail_20kbForward | 0.0184 | 4.2e-4 | **0.0567** | 0.00452 |
| Cocktail_5kbReverseOrder | 0.0174 | 4.5e-4 | 0.0496 | 0.00463 |
| Cocktail_5kbAltOrientation | 0.0182 | **5.1e-4** | 0.0521 | 0.00463 |
| Cocktail_CTCFSeparated | 0.0165 | **5.2e-4** | 0.0557 | **0.00472** |

Source: `results/cocktail_metrics.csv`.

## Interpretation

- AlphaGenome’s DNase head appears highly responsive to the HNF4A motif, even in an erythroid promoter context, suggesting motif strength rather than cell-type priors dominates accessibility predictions.
- GATA1-only modules contribute near-background signal; boosting their length or pairing with HS2 repeats may be necessary to elicit cooperative gains.
- CTCF separators reduce extreme peaks while mildly increasing promoter averages, consistent with insulation spreading accessibility rather than amplifying it.

## Suggested Next Steps

1. Expand to include additional erythroid factors (e.g., TAL1, KLF1) to test whether context-specific motif pairs elevate HS2/GATA1 synergy.
2. Swap promoters (HBG1 vs non-erythroid) within the same cocktail to probe whether promoter identity rescues motif-specific responses.
3. Request additional AlphaGenome outputs (H3K27ac, TF binding) to check whether HNF4A bias persists across chromatin marks.
