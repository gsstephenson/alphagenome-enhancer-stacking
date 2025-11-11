# Regulatory Grammar Experiment - COMPLETE ✓

## 📋 Status: ALL ANALYSIS COMPLETE

**Date**: November 11, 2024  
**Total Constructs**: 66  
**Total Predictions**: 66 (100% complete)  
**Data Generated**: 530 MB  
**Analysis Status**: ✅ Complete  

---

## 🎯 Quick Access

### 📖 Start Here
1. **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Best overview, key findings in tables
2. **[RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)** - Full 310-line detailed analysis
3. **[results/COMPREHENSIVE_SUMMARY.png](results/COMPREHENSIVE_SUMMARY.png)** - Single-page visual summary

### 📊 Key Visualizations

**Core Findings**:
- `results/COMPREHENSIVE_SUMMARY.png` (751K) - All cooperativity results in one figure
- `results/SPACING_SUMMARY.png` (324K) - Distance-activity relationships

**Individual Analyses**:
- `results/cooperativity_additivity_scores.png` (152K)
- `results/spacing_response_curves.png` (297K)
- `results/celltype_heatmap_max_dnase.png` (221K)
- `results/celltype_heatmap_mean_dnase.png` (258K)
- `results/orientation_effects.png` (186K)

### 📈 Data Tables (CSV)
- `results/cooperativity_results.csv` - All 10 pairs with additivity scores
- `results/spacing_results.csv` - Distance vs signal data
- `results/celltype_specificity_results.csv` - Cell-type matching results
- `results/orientation_results.csv` - Strand orientation effects

---

## 🔬 Top-Line Results

### Main Finding
**AlphaGenome shows local motif cooperativity but lacks biological context.**

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Synergistic pairs** | 2/10 (20%) | Only GATA1+KLF1, HS2+KLF1 show synergy |
| **Interference pairs** | 5/10 (50%) | Most pairs are sub-additive |
| **Mean additivity** | 0.90× | Overall sub-additive (not synergistic) |
| **Optimal spacing** | 1000 bp | Local cooperativity peak |
| **Distance sensitivity** | 24% drop over 10kb | Shallow gradient (limited range) |

### Critical Mismatch
⚠️ **HS2+GATA1: Predicted 0.89× (interference), Biology expects synergy!**

This is a **major discrepancy** - HS2 and GATA1 are canonical β-globin regulators that should cooperate.

---

## 📁 File Structure

```
experiments/regulatory_grammar/
├── INDEX.md                              ← You are here
├── QUICK_REFERENCE.md                    ← Best quick overview
├── RESULTS_SUMMARY.md                    ← Full detailed analysis
├── README.md                             ← Original experiment design
│
├── sequences/                            ← Input constructs (66 FASTA files)
│   ├── CellType_*.fa (23 files)
│   ├── Cooperativity_*.fa (17 files)
│   ├── Spacing_*.fa (10 files)
│   └── Orientation_*.fa (16 files)
│
├── alphagenome_outputs/                  ← Raw predictions (530 MB)
│   ├── *_dnase.npy (66 files, 4.1M each)
│   ├── *_dnase_mean.npy (66 files)
│   ├── *_dnase.txt (66 files)
│   └── *_stats.txt (66 files)
│
├── results/                              ← Analysis outputs (2.4 MB)
│   ├── COMPREHENSIVE_SUMMARY.png         ← ⭐ Main summary figure
│   ├── SPACING_SUMMARY.png               ← ⭐ Distance effects
│   ├── cooperativity_additivity_scores.png
│   ├── spacing_response_curves.png
│   ├── celltype_heatmap_max_dnase.png
│   ├── celltype_heatmap_mean_dnase.png
│   ├── orientation_effects.png
│   ├── cooperativity_results.csv
│   ├── spacing_results.csv
│   ├── celltype_specificity_results.csv
│   └── orientation_results.csv
│
├── build_regulatory_grammar_constructs.py  ← Construct builder
├── run_regulatory_grammar_predictions.py   ← AlphaGenome runner
├── analyze_regulatory_grammar.py          ← Analysis script
└── create_summary_figures.py              ← Summary visualizations
```

---

## 🧬 Top 10 Cooperativity Results (Ranked)

| Rank | Pair | Additivity | Type | Biology Match? |
|------|------|------------|------|----------------|
| 🏆 1 | **GATA1+KLF1** | **1.26×** | ✓ Synergy | ✅ Yes |
| 2 | HS2+KLF1 | 1.12× | ✓ Synergy | ✅ Yes |
| 3 | HS2+TAL1 | 0.96× | Independent | ⚠️ Expected synergy |
| 4 | KLF1+TAL1 | 0.96× | Independent | ⚠️ Expected synergy |
| 5 | KLF1+HNF4A | 0.96× | Independent | ✅ Yes (cross-lineage) |
| 6 | **HS2+GATA1** | **0.89×** | ✗ Interference | ❌ **NO - expects synergy!** |
| 7 | GATA1+HNF4A | 0.76× | ✗ Interference | ✅ Yes (cross-lineage) |
| 8 | GATA1+TAL1 | 0.72× | ✗ Interference | ❌ NO - same pathway |
| 9 | HS2+HNF4A | 0.72× | ✗ Interference | ✅ Yes (cross-lineage) |
| 💥 10 | TAL1+HNF4A | 0.66× | ✗ Strong Interference | ✅ Yes (cross-lineage) |

**Legend**: ✅ Prediction matches biology | ⚠️ Unexpected | ❌ Major mismatch

---

## 📊 Spacing Results Summary

| Distance | Signal | % of Optimal | Notes |
|----------|--------|--------------|-------|
| 0 bp | 1.78 | 99% | Nearly optimal (direct contact) |
| 100 bp | 1.73 | 96% | Minimal drop |
| 250 bp | 1.41 | 78% | Local minimum |
| 500 bp | 1.66 | 92% | Recovery |
| 750 bp | 1.80 | 99% | Near-optimal |
| **1000 bp** | **1.80** | **100%** | ⭐ **Optimal** |
| 2 kb | 1.63 | 90% | Gradual decay begins |
| 3 kb | 1.55 | 86% | Continued decay |
| 5 kb | 1.69 | 94% | Slight recovery |
| 10 kb | 1.37 | 76% | 24% drop from optimal |

**Interpretation**: Shallow distance gradient suggests local modeling, not true 3D looping.

---

## 💡 Key Biological Insights

### ✅ What AlphaGenome Does Well
1. **Recognizes TF motifs** - All enhancers produce signals
2. **Local cooperativity** - GATA1+KLF1 shows 1.26× synergy
3. **Short-range effects** - Spacing matters within 1-10kb
4. **Some cross-lineage interference** - TAL1+HNF4A shows 0.66× (correct)

### ❌ What AlphaGenome Gets Wrong
1. **HS2+GATA1 interference** - Biology shows strong synergy
2. **Weak HS2 baseline** - HS2 is strong enhancer in vivo
3. **60% interference rate** - Biology shows more synergy
4. **Limited cell-type specificity** - Can't strongly differentiate contexts

### 🔬 Model Training Hypothesis

**Likely trained on**:
- ✓ TF ChIP-seq (motif patterns)
- ✓ DNase/ATAC-seq (chromatin accessibility)
- ✓ Local sequence composition
- ✓ Short-range (<1kb) interactions

**Likely NOT trained on**:
- ✗ Hi-C (3D structure)
- ✗ eQTL (long-range effects)
- ✗ CRISPR screens (functional validation)
- ✗ Deep cell-type contexts

---

## 🎯 Recommended Next Steps

### For Experimentalists
1. **Validate HS2+GATA1** - Test if true biology shows interference or synergy
2. **Test GATA1+KLF1 synergy** - Validate 1.26× prediction
3. **Clone spacing constructs** - Test 1kb optimal spacing hypothesis
4. **Measure cell-type specificity** - Does HNF4A truly dominate in K562?

### For Computational
1. **Add CTCF analysis** - Quantify insulator effects
2. **Extend spacing** - Test 50kb, 100kb, 1Mb distances
3. **Compare to Enformer** - Does it show same patterns?
4. **Test synthetic motifs** - Isolate pure TF cooperativity

---

## 📚 Related Experiments

### Previous Work
1. **Enhancer Stacking** (`../distance_decay/`)
   - Result: Linear 1-10×, saturation at 40×
   - Consistent with 0.90× mean additivity

2. **Distance Decay** (`../distance_decay/`)
   - Result: NO distance effects (distance-invariant)
   - Current: 24% drop over 10kb (slightly inconsistent)

3. **Heterotypic Cocktail** (`../heterotypic_cocktail/`)
   - Result: HNF4A dominated in K562 (wrong cell)
   - Consistent with poor cell-type specificity

---

## 🔗 Quick Links

### Documentation
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Summary tables and key findings
- [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) - Full detailed analysis

### Key Figures
- [COMPREHENSIVE_SUMMARY.png](results/COMPREHENSIVE_SUMMARY.png) - Main results
- [SPACING_SUMMARY.png](results/SPACING_SUMMARY.png) - Distance effects

### Raw Data
- [cooperativity_results.csv](results/cooperativity_results.csv)
- [spacing_results.csv](results/spacing_results.csv)
- [celltype_specificity_results.csv](results/celltype_specificity_results.csv)
- [orientation_results.csv](results/orientation_results.csv)

---

## ⚡ TL;DR

**Question**: Does AlphaGenome understand transcriptional regulatory grammar?

**Answer**: **Partially.** AlphaGenome recognizes TF motifs and shows local cooperativity within 1kb, but:
- ❌ 50% of pairs show interference (not synergy)
- ❌ HS2+GATA1 predicted as interference (biology shows synergy)
- ❌ Limited cell-type specificity
- ❌ Shallow distance effects (24% over 10kb)

**Conclusion**: AlphaGenome is a **local sequence-to-chromatin model** (good for <1kb), not a full 3D genome simulator (poor for >10kb).

---

**Status**: ✅ **ALL ANALYSIS COMPLETE**  
**Next Action**: Review visualizations and validate key predictions experimentally  
**Contact**: See lab documentation for questions
