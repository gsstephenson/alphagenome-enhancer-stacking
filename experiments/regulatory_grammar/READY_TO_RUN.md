# Regulatory Grammar Experiment - Ready to Run!

**Date:** November 11, 2025  
**Status:** ✅ Constructs Built - Ready for Predictions

---

## Summary

You've successfully created a **comprehensive regulatory grammar test suite** with **66 synthetic constructs** that systematically probe AlphaGenome's understanding of transcriptional regulation.

### What We Built

#### ✅ Part 1: Cell-Type Specificity (23 constructs)
Tests whether AlphaGenome understands context-dependent enhancer activity:
- 3 promoters (HBG1, ALB, CD19) 
- 3 enhancers (HS2, GATA1, KLF1, TAL1, HNF4A)
- 3 cell types (K562, HepG2, GM12878)
- **Key question:** Does HNF4A+ALB score higher in HepG2 than K562?

#### ✅ Part 2: Motif Cooperativity (17 constructs)
Tests TF-TF interactions:
- 5 single enhancers (controls)
- 10 pairwise combinations (6 expected synergy, 4 expected interference)
- 2 CTCF separator tests
- **Key question:** Do erythroid TFs (HS2+GATA1) show super-additive signal?

#### ✅ Part 3: Short-Range Spacing (10 constructs)
Tests optimal TF spacing within model's receptive field:
- HS2 + GATA1 at 10 different spacings (0-10 kb)
- **Key question:** Is there an optimal spacing for cooperation (~200-500 bp)?

#### ✅ Part 4: Orientation Effects (16 constructs)  
Tests strand sensitivity:
- 4 TF pairs × 4 orientation combinations
- Includes CTCF (known to be orientation-dependent)
- **Key question:** Does orientation matter for cooperativity?

---

## File Structure

```
experiments/regulatory_grammar/
├── README.md                                      # Experiment design
├── build_regulatory_grammar_constructs.py         # ✅ COMPLETE
├── run_regulatory_grammar_predictions.py          # 🔴 READY TO RUN
├── analyze_regulatory_grammar.py                  # ⏳ Run after predictions
├── construct_manifest.json                        # ✅ Generated (66 constructs)
├── sequences/                                     # ✅ 66 FASTA files
│   ├── CellType_HBG1_HS2_K562.fa
│   ├── CellType_ALB_HNF4A_HepG2.fa
│   ├── Pair_HS2_GATA1.fa
│   ├── Spacing_500bp_HS2_GATA1.fa
│   ├── Orient_HS2+_CTCF-.fa
│   └── ... (61 more)
├── alphagenome_outputs/                           # ⏳ Will contain predictions
├── results/                                       # ⏳ Will contain analysis
└── logs/                                          # ⏳ Will contain run logs
```

---

## Next Steps

### 1. Set API Key (if not already set)
```bash
export ALPHA_GENOME_KEY="your_api_key_here"
```

### 2. Run Predictions (~1.5 hours, ~$7)
```bash
cd /mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/regulatory_grammar

conda activate alphagenome-env

python run_regulatory_grammar_predictions.py 2>&1 | tee logs/predictions_$(date +%Y%m%d_%H%M%S).log
```

**What this does:**
- Loads each of the 66 constructs
- Calls AlphaGenome API (handles cell-type-specific predictions automatically)
- Saves predictions as `.npy` files
- Generates summary statistics
- Estimated cost: 66 constructs × $0.10 = $6.60
- Estimated time: 66 constructs × 90 seconds = 99 minutes

### 3. Analyze Results (~30 minutes)
```bash
python analyze_regulatory_grammar.py
```

**What this does:**
- Generates cell-type specificity heatmaps
- Calculates cooperativity additivity scores
- Plots spacing response curves
- Compares orientation effects
- Saves all figures and data tables to `results/`

---

## Expected Findings

### Hypothesis 1: AlphaGenome Understands Cell-Type Context
**Test:** Compare HNF4A+ALB in HepG2 vs K562

**Prediction A (Model is smart):**
- Signal higher in HepG2 (correct context)

**Prediction B (Motif-driven, from your cocktail data):**
- Signal similar in both (HNF4A dominates regardless)

**Impact:** Reveals whether cell-type embeddings are meaningful

---

### Hypothesis 2: Erythroid TFs Show Synergy
**Test:** HS2 + GATA1 vs (HS2 alone + GATA1 alone)

**Prediction A (Learned cooperativity):**
- Additivity score > 1.1 (super-additive)

**Prediction B (Additive only, from your data):**
- Additivity score ≈ 1.0 (independent)

**Impact:** Tests if model learned TF-TF interactions

---

### Hypothesis 3: Optimal Spacing Exists
**Test:** HS2 + GATA1 at 0-10 kb spacings

**Prediction A (Short-range interaction):**
- Peak signal at ~200-500 bp

**Prediction B (Distance-invariant, from your data):**
- Flat response across all spacings

**Impact:** Defines model's "interaction radius"

---

### Hypothesis 4: Orientation Matters for CTCF
**Test:** HS2 + CTCF in (+,+) vs (+,-) vs (-,+) vs (-,-)

**Prediction A (Directional insulation):**
- Signal differs by orientation

**Prediction B (Sequence-only model):**
- No difference (reverse complement equivalent)

**Impact:** Tests if model encodes strand-specific biology

---

## Comparison to Your Previous Experiments

| Experiment | Previous | This Experiment | Improvement |
|------------|----------|----------------|-------------|
| **Enhancer stacking** | ✅ Done (9 constructs) | N/A | Baseline established |
| **Distance decay** | ✅ Done (24 constructs) | N/A | No effect found |
| **Cocktails** | ✅ Done (6 constructs) | ✅ Expanded (17 cooperativity) | Systematic TF pairs |
| **Cell-type specificity** | ❌ Not tested | ✅ New (23 constructs) | Tests context |
| **Short-range spacing** | ❌ Not tested | ✅ New (10 constructs) | Within receptive field |
| **Orientation** | 🟡 Partial (5 constructs) | ✅ Expanded (16 constructs) | Systematic test |

---

## Cost & Time Breakdown

### Computational
- **API calls:** 66 constructs × $0.10 = **$6.60**
- **Runtime:** 66 × 90 seconds = **99 minutes (~1.7 hours)**
- **Storage:** 66 × 250 MB = **16.5 GB**

### Your Time
- **Building constructs:** ✅ DONE (30 minutes)
- **Running predictions:** 1.7 hours (mostly waiting)
- **Analysis:** 30 minutes (automated)
- **Writing:** 1-2 days (manuscript draft)

**Total:** ~2 days to complete analysis + 1 week to manuscript

---

## Key Advantages of This Design

1. **Systematic:** Every combination tested, not cherry-picked
2. **Controlled:** Matched construct lengths, identical filler
3. **Quantitative:** Additivity scores, not qualitative descriptions
4. **Falsifiable:** Clear predictions for each hypothesis
5. **Publishable:** Positive OR negative results are valuable

---

## Potential Manuscript Title

> **"Dissecting AlphaGenome's Regulatory Grammar: Cell-Type Specificity, Motif Cooperativity, and the Limits of Sequence-Based Models"**

### Abstract Preview
*We systematically tested AlphaGenome's understanding of transcriptional grammar using 66 synthetic constructs. AlphaGenome showed [motif-driven / context-aware] predictions with [additive / synergistic] TF interactions. Cell-type specificity was [present / limited], and optimal TF spacing was [detected / absent]. These findings define appropriate applications for sequence-based epigenome models and highlight the need for 3D genome integration.*

---

## When to Run This

### Option A: Run Now (Recommended)
- Get results while constructs are fresh in mind
- Parallel track with other experiments
- Results in ~2 hours

### Option B: Run in Batch Later
- Combine with future experiments
- Save on API warm-up time
- Run overnight

**Recommendation:** Run now - the constructs are built and validated, and results will inform your overall strategy.

---

## Success Criteria

### Minimum Viable Result
- ✅ All 66 constructs successfully predicted
- ✅ Cell-type effects quantified (even if weak)
- ✅ Cooperativity scores calculated
- ✅ Spacing curve generated
- ✅ Orientation effects measured

### Bonus Achievements
- 🌟 Clear cell-type specificity detected
- 🌟 TF synergy vs interference patterns
- 🌟 Optimal spacing identified
- 🌟 CTCF orientation dependence
- 🌟 Technical replicates added (if time)

---

## Troubleshooting

### If predictions fail:
```bash
# Check API key
echo $ALPHA_GENOME_KEY

# Test single construct
conda run -n alphagenome-env python -c "
from alphagenome import AlphaGenomeClient
client = AlphaGenomeClient()
print('API key valid!')
"
```

### If analysis fails:
```bash
# Check predictions exist
ls alphagenome_outputs/*.npy | wc -l
# Should see 66 files

# Run analysis on subset
# Edit analyze_regulatory_grammar.py to skip missing files
```

---

## Ready to Proceed?

Run this command to start predictions:

```bash
cd /mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/regulatory_grammar

conda activate alphagenome-env

python run_regulatory_grammar_predictions.py 2>&1 | tee logs/predictions_$(date +%Y%m%d_%H%M%S).log
```

**This will take ~1.7 hours and cost ~$7. You can monitor progress in real-time.**

---

🚀 **Ready to test AlphaGenome's regulatory grammar understanding!**
