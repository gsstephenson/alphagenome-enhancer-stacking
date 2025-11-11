# Logic Gates Experiment - Quick Start Guide

## 🚀 Overview

This experiment tests whether AlphaGenome understands Boolean logic in transcriptional regulation by systematically testing AND/OR/NOT/XOR gate behavior.

## 📋 Prerequisites

1. **Python packages**: numpy, pandas, matplotlib, seaborn, scikit-learn
2. **AlphaGenome access**: Update `ALPHAGENOME_SCRIPT` path in `run_logic_gate_predictions.py`
3. **Sequence files**: Enhancers and promoters in `../../sequences/`
4. **Filler sequence**: `../../filler/1M_filler.txt`

## ⚡ Quick Start (3 Commands)

```bash
# 1. Build all constructs (~52 FASTA files)
python build_logic_gate_constructs.py

# 2. Run predictions (may take hours)
python run_logic_gate_predictions.py

# 3. Analyze and visualize
python analyze_logic_gates.py
python create_logic_gate_figures.py
```

**Or run everything at once:**
```bash
chmod +x run_full_pipeline.sh
./run_full_pipeline.sh
```

## 📁 Output Structure

```
logic_gates/
├── sequences/              # FASTA files for all constructs
│   └── LogicGate_AND_GATA1_KLF1_11_K562.fasta
├── alphagenome_outputs/    # Prediction results
│   └── LogicGate_AND_GATA1_KLF1_11_K562_dnase.npy
├── results/
│   ├── figures/           # All visualizations
│   │   ├── summary_dashboard.png
│   │   ├── all_truth_tables.png
│   │   ├── confusion_matrix.png
│   │   └── logic_gate_*.png
│   └── logic_gate_summary.csv
├── logic_gate_manifest.json   # All construct metadata
└── logic_gate_analysis.json   # Detailed results
```

## 🔬 What Gets Tested

### 13 TF Pairs × 4 Conditions = 52 Constructs

**AND Gates (5 pairs)**
- GATA1 × KLF1 (your 1.26× synergy finding)
- HS2 × GATA1 (your 0.89× interference - retest)
- HS2 × KLF1
- GATA1 × TAL1
- OCT4 × SOX2 (if available)

**OR Gates (4 pairs)**
- GATA1 × GATA2
- KLF1 × TAL1
- HS2 × TAL1
- GATA1 × KLF1 (compare to AND behavior)

**XOR Gates (4 pairs)**
- GATA1 × PU.1 (THE canonical XOR)
- GATA1 × HNF4A (your 0.76× interference)
- TAL1 × HNF4A (your 0.66× interference)
- HS2 × HNF4A

### 4 Truth Table Conditions Per Pair

| Code | TF_A | TF_B | Description |
|------|------|------|-------------|
| 00 | ✗ | ✗ | Baseline (neither TF) |
| 01 | ✗ | ✓ | Only TF_B |
| 10 | ✓ | ✗ | Only TF_A |
| 11 | ✓ | ✓ | Both TFs |

## 📊 Key Metrics

### 1. Logic Scores (R² with ideal pattern)
- **1.0** = Perfect match to expected gate
- **0.5** = Random
- **0.0** = No correlation

### 2. Synergy Metrics
- **Additivity ratio** = signal(A+B) / [signal(A) + signal(B)]
  - \>1.1 = Synergy
  - 0.9-1.1 = Additive
  - <0.9 = Interference

### 3. Classification Accuracy
- % of pairs correctly matching expected gate type

## 🎯 Expected Outcomes

### Success Scenario
- AND gates: High R² for AND pattern, synergy >1.1
- OR gates: High R² for OR pattern, saturation
- XOR gates: High R² for XOR pattern, interference when both present

### Null Scenario
- All pairs show simple additivity regardless of gate type
- Random classification (25% accuracy)
- No biological context understanding

## 🔍 Interpreting Results

### Check Summary Dashboard First
```
results/figures/summary_dashboard.png
```

Key plots:
1. **Classification accuracy** - Are expected gates predicted correctly?
2. **Mean logic scores** - How well do gates match ideal patterns?
3. **Synergy distribution** - Are most pairs additive/synergistic/interfering?
4. **Example truth tables** - Visual inspection of best examples

### Dive Deeper
- `confusion_matrix.png` - Where do misclassifications occur?
- `logic_scores_comparison.png` - Which pairs have highest scores?
- `synergy_analysis.png` - Does synergy predict logic score?

## 🐛 Troubleshooting

### "Module not found: OCT4, PU1"
Some TFs may not be in your enhancer library. The script will skip missing modules and continue with available ones.

**Solution**: Either add sequences or comment out those gate definitions in `build_logic_gate_constructs.py`

### "AlphaGenome script not found"
You need to configure AlphaGenome API access.

**Solution**: Edit `run_logic_gate_predictions.py` line 16:
```python
ALPHAGENOME_SCRIPT = Path("/your/path/to/alphagenome/predict.py")
```

Or use dry-run mode:
```bash
python run_logic_gate_predictions.py --dry-run
```

### "Missing prediction data"
Predictions haven't completed yet.

**Solution**: 
```bash
# Check how many predictions exist
ls alphagenome_outputs/*.npy | wc -l

# Should be ~52 files (one per construct)
```

## 📈 Pilot Experiment (Fast Test)

Test with just 3 pairs (12 constructs):

```bash
# Edit build_logic_gate_constructs.py
# In define_logic_gates(), keep only first AND, OR, XOR pair

python build_logic_gate_constructs.py
python run_logic_gate_predictions.py --limit 12
python analyze_logic_gates.py
python create_logic_gate_figures.py
```

Should complete in <1 hour depending on AlphaGenome speed.

## 📝 Manuscript Outline

### Title
"Boolean Logic Gates in Deep Learning Models of Gene Regulation"

### Key Findings to Report
1. Classification accuracy vs random (25%)
2. Which gate types work best (AND? OR?)
3. Correlation with synergy metrics
4. Specific TF pairs that match/violate expectations
5. Comparison to biological literature

### Figures
- Figure 1: Experimental design schematic
- Figure 2: Summary dashboard (main results)
- Figure 3: Selected truth tables (best examples)
- Figure 4: Confusion matrix + synergy analysis
- Supplementary: All individual pairs

## 🔗 Related Experiments

This builds on your regulatory grammar work:
- Your GATA1+KLF1 synergy → Test as AND gate
- Your TAL1+HNF4A interference → Test as XOR gate
- Your HS2+GATA1 paradox → Retest in logic framework

## 💡 Extensions

1. **Test in multiple cell types** - Does logic change by context?
2. **Vary spacing** - Do gates work at 1kb? 10kb?
3. **Add more gate types** - NAND, NOR, XNOR
4. **Compare to Enformer** - Do different models learn same logic?

## 🆘 Need Help?

Check existing regulatory grammar experiment for examples:
```bash
cd ../regulatory_grammar
cat README.md
cat RESULTS_SUMMARY.md
```

---

**Questions?** Check the detailed README.md or review existing analysis scripts.
