# Regulatory Grammar Experiment: Comprehensive Test of AlphaGenome

**Date:** November 11, 2025  
**Investigator:** Grant Stephenson  
**Lab:** Layer Laboratory, CU Boulder  
**Status:** 🚀 In Progress

---

## Objective

Systematically test whether AlphaGenome understands **transcriptional grammar rules**:

1. **Cell-type specificity** - Do promoter-enhancer pairs show context-dependent activity?
2. **Motif cooperativity** - Do TF pairs show synergy vs interference?
3. **Spacing sensitivity** - What's the optimal distance for cooperation?
4. **Orientation dependence** - Does strand matter for TF binding?

---

## Experimental Design

### Part 1: Cell-Type Specificity Matrix (27 constructs)

Test 3 promoters × 3 enhancers × 3 cell types = 27 predictions

| Promoter | Cell Type | Enhancer | Expected Result |
|----------|-----------|----------|-----------------|
| HBG1 (erythroid) | K562 | HS2 (erythroid) | ✅ HIGH - correct context |
| HBG1 (erythroid) | K562 | HNF4A (hepatic) | ❌ LOW - wrong enhancer |
| HBG1 (erythroid) | K562 | GATA1 (erythroid) | ✅ HIGH - correct context |
| ALB (hepatic) | HepG2 | HS2 (erythroid) | ❌ LOW - wrong enhancer |
| ALB (hepatic) | HepG2 | HNF4A (hepatic) | ✅ HIGH - correct context |
| ALB (hepatic) | HepG2 | GATA1 (erythroid) | ❌ LOW - wrong enhancer |
| CD19 (B-cell) | GM12878 | HS2 (erythroid) | ❌ LOW - wrong enhancer |
| CD19 (B-cell) | GM12878 | HNF4A (hepatic) | ❌ LOW - wrong enhancer |
| CD19 (B-cell) | GM12878 | GATA1 (erythroid) | ❌ LOW - wrong enhancer |

**Hypothesis:** Cell-type-specific combinations should show higher signal than mismatched combinations.

**Control from previous cocktail experiment:** HNF4A dominated even in K562 (wrong cell type), suggesting motif strength > context.

---

### Part 2: Pairwise Motif Cooperativity (15 constructs)

Test all pairwise combinations of erythroid TFs at 5 kb spacing

| Pair | Expected Interaction | Rationale |
|------|---------------------|-----------|
| HS2 + GATA1 | Synergy | Both erythroid, known co-binding |
| HS2 + KLF1 | Synergy | Both activate β-globin locus |
| HS2 + TAL1 | Synergy | Both erythroid master regulators |
| GATA1 + KLF1 | Synergy | Cooperative at β-globin |
| GATA1 + TAL1 | Synergy | Part of erythroid TF complex |
| KLF1 + TAL1 | Synergy | Co-activate erythroid genes |
| HS2 + HNF4A | Interference | Different cell types |
| GATA1 + HNF4A | Interference | Competing contexts |
| KLF1 + HNF4A | Interference | Unrelated pathways |
| TAL1 + HNF4A | Interference | Unrelated pathways |

**Also test with CTCF separator:**
- HS2 + [CTCF] + GATA1 - Does insulator block cooperation?
- HS2 + [CTCF] + HNF4A - Control (already independent)

**Plus individual controls (5 constructs):**
- HS2 alone
- GATA1 alone
- KLF1 alone
- TAL1 alone
- HNF4A alone

**Additivity Score Calculation:**
```python
score = signal(E1+E2) / (signal(E1_alone) + signal(E2_alone))
> 1.1: Synergy (super-additive)
0.9-1.1: Independent (additive)
< 0.9: Interference (sub-additive)
```

---

### Part 3: Short-Range Spacing Screen (20 constructs)

Test HS2 + GATA1 at fine-grained spacings (within model's receptive field)

**Spacings:** 0, 100, 250, 500, 750, 1000, 2000, 3000, 5000, 10000 bp

**Design:**
```
[Filler] - HS2 - [X bp] - GATA1 - [100kb - X] - HBG1_Promoter - [Filler]
```

**Expected:** Peak cooperation at ~200-500 bp (typical TF spacing)

**Contrast with distance decay result:** No effect at 1-500 kb (long-range), but this tests <10 kb (short-range).

---

### Part 4: Orientation Effects (16 constructs)

Test 4 TF pairs × 4 orientation combinations

**Test pairs:**
- HS2 + GATA1 (expect cooperativity)
- HS2 + HNF4A (expect independence)
- GATA1 + KLF1 (expect cooperativity)
- HS2 + CTCF (expect CTCF orientation matters)

**Orientations:**
- (+, +) Both forward
- (+, -) E1 forward, E2 reverse
- (-, +) E1 reverse, E2 forward
- (-, -) Both reverse

**Hypothesis:** CTCF orientation will matter (known directional insulator), others less so.

---

## Total Constructs

| Experiment | Count | Status |
|------------|-------|--------|
| Cell-type matrix | 27 | 🔴 To build |
| Pairwise cooperativity | 15 | 🔴 To build |
| Short-range spacing | 20 | 🔴 To build |
| Orientation effects | 16 | 🔴 To build |
| **TOTAL** | **78** | **Ready to start** |

---

## Pipeline

### Step 1: Build Constructs
```bash
cd experiments/regulatory_grammar
python build_regulatory_grammar_constructs.py
```

### Step 2: Run Predictions
```bash
python run_regulatory_grammar_predictions.py
```

### Step 3: Analyze Results
```bash
python analyze_regulatory_grammar.py
```

---

## Expected Findings

### If AlphaGenome Understands Grammar:

✅ **Cell-type specificity:** HBG1+HS2 in K562 > ALB+HS2 in K562  
✅ **Erythroid synergy:** HS2+GATA1 shows super-additive signal  
✅ **Optimal spacing:** Peak at ~200-500 bp for HS2+GATA1  
✅ **CTCF orientation:** Signal differs for (+) vs (-) CTCF  

### If AlphaGenome Is Motif-Driven Only:

❌ **Cell-type invariance:** HNF4A dominates everywhere (seen in cocktail experiment)  
❌ **Additivity only:** All pairs sum linearly (no synergy/interference)  
❌ **Spacing-invariant:** Flat response across all distances  
❌ **Orientation-insensitive:** No difference between (+) and (-) orientations  

---

## Key Questions

1. **Does cell-type embedding matter?** Or just motif strength?
2. **Did AlphaGenome learn TF cooperativity?** Or just additive rules?
3. **What's the model's interaction radius?** (Previous: no effect 1-500 kb)
4. **Does orientation encode biology?** Or just sequence features?

---

## Success Metrics

**Publishable findings require:**
- ✅ Clear positive or negative results for each question
- ✅ Technical replicates (lesson from distance decay)
- ✅ Statistical tests (ANOVA, correlations, effect sizes)
- ✅ Comparison to known biology (ChIP-seq, published cooperativity)

---

## Files

```
experiments/regulatory_grammar/
├── README.md                              # This file
├── build_regulatory_grammar_constructs.py # Construct builder
├── run_regulatory_grammar_predictions.py  # Batch predictor
├── analyze_regulatory_grammar.py          # Statistical analysis
├── construct_manifest.json                # Metadata
├── sequences/                             # FASTA files
├── alphagenome_outputs/                   # Predictions
├── results/                               # Analysis outputs
└── logs/                                  # Run logs
```

---

## Timeline

- **Day 1:** Build all 78 constructs (~2 hours)
- **Day 2:** Run predictions (~2 hours API time, $8 cost)
- **Day 3:** Analyze cell-type specificity
- **Day 4:** Analyze cooperativity + spacing
- **Day 5:** Analyze orientation effects + write up

**Total:** 1 week to complete analysis

---

## Next Command

```bash
python build_regulatory_grammar_constructs.py
```
