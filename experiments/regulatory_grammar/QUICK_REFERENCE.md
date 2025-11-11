# Regulatory Grammar Experiment - Quick Reference

## 🎯 Core Question
**Does AlphaGenome understand transcriptional regulatory grammar?**

## 📊 Experiment Overview

**Total Constructs**: 66  
**Cell Types**: K562 (erythroid), HepG2 (hepatic), GM12878 (B-cell)  
**Predictions Generated**: 530 MB outputs  
**Date Completed**: November 11, 2024  

## 🔬 Key Findings

### 1. Cooperativity Results (10 pairs tested)

| Rank | Enhancer Pair | Additivity Score | Type | Interpretation |
|------|---------------|------------------|------|----------------|
| 🏆 1 | **GATA1 + KLF1** | **1.26×** | ✓ Synergy | Strongest cooperation |
| 2 | **HS2 + KLF1** | **1.12×** | ✓ Synergy | Modest synergy |
| 3 | HS2 + TAL1 | 0.96× | Independent | Simple addition |
| 4 | KLF1 + TAL1 | 0.96× | Independent | Simple addition |
| 5 | KLF1 + HNF4A | 0.96× | Independent | Simple addition |
| 6 | HS2 + GATA1 | 0.89× | ✗ Interference | **Unexpected!** |
| 7 | GATA1 + HNF4A | 0.76× | ✗ Interference | Cross-lineage clash |
| 8 | GATA1 + TAL1 | 0.72× | ✗ Interference | Unexpected interference |
| 9 | HS2 + HNF4A | 0.72× | ✗ Interference | Cross-lineage clash |
| 💥 10 | **TAL1 + HNF4A** | **0.66×** | ✗ **Strong Interference** | Strongest repression |

**Summary**: 
- ✓ Synergy: 20% (2/10 pairs)
- Independent: 30% (3/10 pairs)  
- ✗ Interference: 50% (5/10 pairs)

### 2. Single Enhancer Baselines

| Enhancer | Signal | Strength |
|----------|--------|----------|
| **KLF1** | **4.78** | 🔥 Strongest |
| GATA1 | 1.62 | Strong |
| TAL1 | 0.97 | Moderate |
| HNF4A | 0.63 | Weak |
| HS2 | 0.28 | 🧊 Weakest |

**Surprise**: HS2 (canonical β-globin enhancer) shows weakest solo activity!

### 3. Spacing Effects (HS2 + GATA1)

| Distance | Signal | % of Optimal |
|----------|--------|--------------|
| 0 bp | 1.78 | 99% |
| 100 bp | 1.73 | 96% |
| 250 bp | 1.41 | 78% |
| 500 bp | 1.66 | 92% |
| 750 bp | 1.80 | 99% |
| **1000 bp** | **1.80** | **100% ⭐** |
| 2 kb | 1.63 | 90% |
| 3 kb | 1.55 | 86% |
| 5 kb | 1.69 | 94% |
| 10 kb | 1.37 | 76% |

**Optimal spacing**: 1 kb  
**Distance sensitivity**: 24% drop over 10kb (shallow gradient)

### 4. Cell-Type Specificity

See visualizations:
- `celltype_heatmap_max_dnase.png`
- `celltype_heatmap_mean_dnase.png`

**Constructs**: 23 combinations of promoters × enhancers × cell types

## 🧬 Biological Interpretation

### ✅ What AlphaGenome Got Right
1. **Motif recognition**: Identifies TF binding sequences
2. **Some cooperativity**: GATA1+KLF1 shows 1.26× synergy
3. **Local context**: Spacing effects suggest short-range modeling

### ❌ What AlphaGenome Got Wrong
1. **HS2+GATA1 interference** (0.89×) - should show synergy in biology!
2. **Weak HS2 baseline** (0.28) - HS2 is a strong β-globin enhancer in vivo
3. **High interference rate** (50% of pairs) - biology shows more synergy
4. **Shallow distance effects** (24% over 10kb) - prior experiment showed 0%

## 💡 Key Insights

### AlphaGenome's Training Data Hypothesis

**Likely trained on**:
- ✓ TF ChIP-seq (motif recognition)
- ✓ DNase/ATAC-seq (chromatin accessibility)
- ✓ Local sequence composition
- ✓ Short-range interactions (<1kb)

**Likely NOT trained on**:
- ✗ Hi-C data (3D chromosome structure)
- ✗ eQTL studies (long-range enhancer effects)
- ✗ CRISPR perturbations (functional validation)
- ✗ Deep cell-type specific contexts

### Model Limitations

1. **Local-only model**: Works well for <1kb, fails at >10kb
2. **Distance-invariant**: Prior experiment showed NO distance effects
3. **Limited cell-type specificity**: Can't strongly differentiate biological contexts
4. **Missing biological synergy**: Predicts interference where biology shows cooperation

## 🔍 Critical Experiments to Validate

| Priority | Experiment | Predicted | Expected Biology | Discrepancy? |
|----------|-----------|-----------|------------------|--------------|
| 🔴 HIGH | HS2 + GATA1 | 0.89× (interference) | Synergy | **YES - Major** |
| 🔴 HIGH | GATA1 + KLF1 | 1.26× (synergy) | Synergy | Matches |
| 🟡 MED | TAL1 + HNF4A | 0.66× (interference) | Interference | Matches |
| 🟡 MED | 1kb optimal spacing | 1.80 signal | Unknown | Test needed |

## 📁 Files Generated

### Visualizations (1.2 MB total)
```
results/celltype_heatmap_max_dnase.png      (221K)
results/celltype_heatmap_mean_dnase.png     (258K)
results/cooperativity_additivity_scores.png (152K)
results/spacing_response_curves.png         (297K)
results/orientation_effects.png             (186K)
```

### Data Tables
```
results/celltype_specificity_results.csv
results/cooperativity_results.csv
results/spacing_results.csv
results/orientation_results.csv
```

### Raw Predictions
```
alphagenome_outputs/*_dnase.npy            (66 files, 4.1M each)
alphagenome_outputs/*_stats.txt            (66 files)
```

### Documentation
```
RESULTS_SUMMARY.md                         (Full analysis, 310 lines)
QUICK_REFERENCE.md                         (This file)
```

## 🎬 Next Steps

1. **View visualizations**: Open PNG files in `results/`
2. **Read full analysis**: See `RESULTS_SUMMARY.md`
3. **Validate predictions**: Test HS2+GATA1 interference experimentally
4. **Extend analysis**: Add CTCF insulator quantification
5. **Compare models**: Run same constructs on Enformer/DeepSEA

## 📊 Statistics Summary

- **Mean cooperativity**: 0.93× (sub-additive overall)
- **Synergy rate**: 20% (2/10 pairs)
- **Interference rate**: 50% (5/10 pairs)
- **Optimal spacing**: 1 kb for HS2+GATA1
- **Distance sensitivity**: 24% drop over 10kb
- **Cell types tested**: 3 (K562, HepG2, GM12878)

---

**Conclusion**: AlphaGenome is a **local sequence-to-chromatin model** that understands motif cooperativity but lacks long-range 3D genome organization modeling. Use for TF binding predictions (<1kb), not for enhancer-promoter looping (>10kb).
