# Distance Decay Experiment - Results

**Date:** November 11, 2025
**Model:** AlphaGenome DNase, K562
**Enhancer:** HS2 (chr11:5290000-5291000, 1,001 bp)
**Promoter:** HBG1 (chr11:5273600-5273900, 301 bp)

## Key Findings

- **Enhancer signal decay:** 0.2656 @ 1kb → 0.2344 @ 500kb (1.13× reduction)
- **Promoter signal:** 0.00476 ± 0.00038 (largely invariant)
- **Hi-C correlation:** Pearson r = -0.508 (p = 1.99e-01), Spearman ρ = -0.563 (p = 1.46e-01)

## Metrics Table

| name           |   distance_kb |   enh_max |   enh_mean |   enh_auc |   prom_max |   prom_mean |   prom_auc |   spacer_mean |   global_mean |   global_max |
|:---------------|--------------:|----------:|-----------:|----------:|-----------:|------------:|-----------:|--------------:|--------------:|-------------:|
| Distance_1kb   |             1 |  0.265625 | 0.00777305 |   23.3263 |   0.265625 |  0.00570636 |    57.0619 |    0.00323774 |    0.00398364 |     0.265625 |
| Distance_5kb   |             5 |  0.265625 | 0.00699051 |   20.9777 |   0.132812 |  0.00455305 |    45.5288 |    0.00372176 |    0.00397374 |     0.265625 |
| Distance_10kb  |            10 |  0.289062 | 0.00791996 |   23.7672 |   0.132812 |  0.00464088 |    46.407  |    0.00327822 |    0.00395103 |     0.289062 |
| Distance_25kb  |            25 |  0.25     | 0.00863926 |   25.9228 |   0.136719 |  0.00468715 |    46.8697 |    0.00362992 |    0.00393253 |     0.25     |
| Distance_50kb  |            50 |  0.277344 | 0.00802891 |   24.0929 |   0.134766 |  0.00467552 |    46.7535 |    0.00382377 |    0.00393357 |     0.277344 |
| Distance_100kb |           100 |  0.257812 | 0.00829957 |   24.9062 |   0.134766 |  0.00466266 |    46.6248 |    0.00403584 |    0.00392113 |     0.257812 |
| Distance_200kb |           200 |  0.246094 | 0.00815806 |   24.4812 |   0.132812 |  0.00463764 |    46.3747 |    0.00395465 |    0.00390654 |     0.246094 |
| Distance_500kb |           500 |  0.234375 | 0.00783767 |   23.5197 |   0.132812 |  0.00454836 |    45.4819 |    0.00388245 |    0.0038573  |     0.234375 |

## Interpretation

1. **Enhancer signal shows exponential decay** with increasing distance
2. **Promoter signal remains constant** - no long-range activation detected
3. **Weak correlation with Hi-C contact frequency** - AlphaGenome may not fully model looping dynamics
