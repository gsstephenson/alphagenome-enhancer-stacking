# Quick Start: Running Regulatory Grammar Predictions

## Fixed Import Issue ✅

The script now uses the correct AlphaGenome API:
```python
from alphagenome.models import dna_client  # Correct
# NOT: from alphagenome import AlphaGenomeClient  # Wrong
```

## Step-by-Step Instructions

### 1. Test API Setup
```bash
cd /mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/regulatory_grammar

# Make sure you're in the environment
conda activate alphagenome-env

# Test everything works
python test_api.py
```

**Expected output:**
```
✓ AlphaGenome package imported successfully
✓ API key found (length: XX)
✓ AlphaGenome client created successfully
✓ Found manifest with 66 constructs
✓ Found 66 FASTA files

Testing API with a short sequence...
✓ API test successful! Shape: (16384,)
  Max DNase: 0.XXXXXX
  Mean DNase: 0.XXXXXX

✅ All checks passed! Ready to run predictions.
```

If the test fails, check:
- API key is set: `echo $ALPHA_GENOME_KEY`
- Environment is active: `which python` should show `alphagenome-env`

### 2. Run Full Predictions

```bash
# This will take ~1.7 hours and cost ~$7
python run_regulatory_grammar_predictions.py 2>&1 | tee logs/predictions_$(date +%Y%m%d_%H%M%S).log
```

**What happens:**
- Processes all 66 constructs
- For cell-type experiments: runs specific cell type (K562, HepG2, or GM12878)
- For all others: runs K562 only
- Saves predictions to `alphagenome_outputs/`
- Saves summary statistics
- Creates log file

**Progress tracking:**
The script prints status for each construct:
```
[1/66] CellType_HBG1_HS2_K562
  Experiment: cell_type_specificity
  Sequence length: 1,048,576 bp
  Cell types: K562
  Predicting CellType_HBG1_HS2_K562 in K562...
    ✓ Max DNase: 0.XXXX
    ✓ Mean DNase: 0.XXXXXX
```

**If a prediction fails:**
The script will continue with the next construct and report failures at the end.

### 3. Analyze Results

```bash
python analyze_regulatory_grammar.py
```

**Generates:**
- `results/celltype_heatmap_max_dnase.png` - Cell-type specificity
- `results/cooperativity_additivity_scores.png` - TF synergy/interference
- `results/spacing_response_curves.png` - Optimal spacing
- `results/orientation_effects.png` - Orientation dependence
- `results/*.csv` - All data tables

---

## Troubleshooting

### Import Error
```
ImportError: cannot import name 'AlphaGenomeClient'
```
**Fixed!** Script now uses correct import: `from alphagenome.models import dna_client`

### API Key Error
```
ValueError: ALPHA_GENOME_API_KEY or ALPHA_GENOME_KEY environment variable not set
```
**Solution:**
```bash
export ALPHA_GENOME_KEY="your_key_here"
# Or add to ~/.bashrc for persistence
```

### Sequence Length Error
```
Error: Sequence length must be power of 2
```
**Solution:** All constructs are already 1,048,576 bp (2^20). This shouldn't happen.

### API Rate Limit
```
Error: Rate limit exceeded
```
**Solution:** Script has 1-second delays between requests. If still hitting limits, increase sleep time in script.

---

## Expected Runtime

| Stage | Time | Cost |
|-------|------|------|
| Test API | 30 sec | $0 |
| Cell-type (23) | ~35 min | $2.30 |
| Cooperativity (17) | ~26 min | $1.70 |
| Spacing (10) | ~15 min | $1.00 |
| Orientation (16) | ~24 min | $1.60 |
| **TOTAL** | **~1.7 hours** | **~$6.60** |

---

## What You'll Learn

### Question 1: Cell-Type Specificity
Does HNF4A + ALB show higher signal in HepG2 (hepatocyte) than K562 (erythroid)?

**Your cocktail data suggests:** HNF4A dominates everywhere (motif strength > context)

### Question 2: Motif Cooperativity
Do erythroid TFs (HS2, GATA1, KLF1, TAL1) show synergy?

**Your cocktail data suggests:** Near-additive, limited synergy

### Question 3: Optimal Spacing
Is there a peak in signal at ~200-500 bp for HS2 + GATA1?

**Your distance decay data suggests:** Flat response (distance-invariant)

### Question 4: Orientation Matters?
Does CTCF show different signal in (+) vs (-) orientation?

**Expected:** Yes, CTCF is directional insulator

---

## Ready to Run!

```bash
# Test first
python test_api.py

# If successful, run predictions
python run_regulatory_grammar_predictions.py 2>&1 | tee logs/predictions_$(date +%Y%m%d_%H%M%S).log
```

**You can monitor progress in real-time or let it run overnight.**
