# Logic Gates Experiment - Setup Complete! 🎉

## ✅ What's Been Created

Your complete Boolean logic gates experiment framework is ready to run!

### 📁 File Structure
```
logic_gates/
├── README.md                           # Full experiment documentation
├── QUICKSTART.md                       # Quick start guide (READ THIS FIRST!)
├── run_full_pipeline.sh                # Master pipeline script
├── build_logic_gate_constructs.py      # Step 1: Generate FASTA sequences
├── run_logic_gate_predictions.py       # Step 2: Run AlphaGenome
├── analyze_logic_gates.py              # Step 3: Compute logic scores
└── create_logic_gate_figures.py        # Step 4: Generate visualizations
```

## 🔧 Fixed Issues

**Problem:** Module names didn't match loaded sequences  
**Solution:** Updated all TF names to match your actual files:
- `GATA1` → `GATA1_MODULE`
- `KLF1` → `KLF1_MODULE`
- `TAL1` → `TAL1_MODULE`
- `HS2` → `HS2_ENHANCER`
- `HNF4A` → `HNF4A_MODULE`
- `HBG1` → `HBG1_PROMOTER`

## 🚀 Quick Start (3 Steps)

### 1. Build Constructs
```bash
python build_logic_gate_constructs.py
```
**Expected output:** ~44 FASTA files (11 pairs × 4 conditions, minus missing TFs)
- AND gates: 5 pairs
- OR gates: 3 pairs (GATA2 not available)
- XOR gates: 3 pairs (PU1 not available)

### 2. Configure AlphaGenome
Edit `run_logic_gate_predictions.py` line 16:
```python
ALPHAGENOME_SCRIPT = Path("/your/actual/path/to/alphagenome_predict.py")
```

Or look at your regulatory grammar experiment to see how you ran predictions.

### 3. Run Pipeline
```bash
# Option A: Run everything
./run_full_pipeline.sh

# Option B: Step by step
python build_logic_gate_constructs.py
python run_logic_gate_predictions.py
python analyze_logic_gates.py
python create_logic_gate_figures.py
```

## 📊 What You'll Test

### Available TF Pairs (11 total)

**AND Gates (5 pairs):**
1. ✅ GATA1 + KLF1 (your 1.26× synergy)
2. ✅ GATA1 + TAL1 (erythroid cooperation)
3. ✅ HS2 + GATA1 (your 0.89× paradox - retest!)
4. ✅ HS2 + KLF1 (your 1.12× synergy)
5. ⚠️ OCT4 + SOX2 (will skip - OCT4/SOX2 sequences not available)

**OR Gates (3 pairs):**
6. ✅ KLF1 + TAL1
7. ✅ GATA1 + KLF1 (compare AND vs OR behavior!)
8. ✅ HS2 + TAL1
9. ⚠️ GATA1 + GATA2 (will skip - GATA2 not available)

**XOR Gates (3 pairs):**
10. ✅ GATA1 + HNF4A (your 0.76× interference)
11. ✅ TAL1 + HNF4A (your 0.66× strong interference)
12. ✅ HS2 + HNF4A (cross-lineage)
13. ⚠️ GATA1 + PU1 (will skip - PU1 not available)

**Total: 44 constructs** (11 pairs × 4 conditions)

## 🎯 Expected Results

### Success Indicators
- **AND gates show synergy**: GATA1+KLF1, HS2+KLF1 should have signal(11) >> signal(01) + signal(10)
- **OR gates show saturation**: signal(11) ≈ max(signal(01), signal(10))
- **XOR gates show interference**: signal(11) < signal(01) AND signal(11) < signal(10)

### Key Comparisons
1. **GATA1+KLF1 as AND vs OR**: Same pair, different expected behavior
2. **HS2+GATA1 retest**: Your 0.89× result was surprising - does it fit XOR better?
3. **Cross-lineage pairs**: Do TAL1+HNF4A and GATA1+HNF4A behave like XOR gates?

## 📈 Output Files

After running the full pipeline:

```
logic_gates/
├── sequences/                          # Your 44 FASTA files
│   └── LogicGate_AND_GATA1_MODULE_KLF1_MODULE_11_K562.fasta
├── alphagenome_outputs/                # Prediction arrays
│   └── LogicGate_AND_GATA1_MODULE_KLF1_MODULE_11_K562_dnase.npy
├── results/
│   ├── logic_gate_summary.csv         # Table of all results
│   └── figures/
│       ├── summary_dashboard.png      # ⭐ START HERE
│       ├── all_truth_tables.png       # All 11 pairs at a glance
│       ├── confusion_matrix.png       # Classification accuracy
│       ├── synergy_analysis.png       # Additivity vs logic score
│       └── logic_gate_*.png           # Individual pair details
└── logic_gate_analysis.json           # Detailed metrics
```

## 🔍 Analysis Metrics

For each TF pair, you'll get:

1. **Truth Table**: 2×2 heatmap of signals [00, 01, 10, 11]
2. **Logic Scores**: R² fit to ideal AND/OR/NOT/XOR patterns
3. **Best Fit Gate**: Which pattern matches best?
4. **Additivity Ratio**: signal(11) / [signal(01) + signal(10)]
5. **Synergy Class**: Synergy, Additive, or Interference

## 💡 Novel Findings This Could Reveal

1. **Does AlphaGenome understand Boolean logic?**
   - Random = 25% accuracy (4 gate types)
   - Good = >50% accuracy
   - Excellent = >75% accuracy

2. **Which gates work best?**
   - Hypothesis: AND/OR easier than XOR/NOT

3. **Does synergy predict logic?**
   - Your prior finding: GATA1+KLF1 synergy, TAL1+HNF4A interference
   - Test: Does synergy → AND gate, interference → XOR gate?

4. **Context dependence**
   - Test same pairs (GATA1+HNF4A) in K562 vs HepG2
   - Does cell type change the logic?

## 🐛 Troubleshooting

### "Missing module: OCT4, SOX2, PU1, GATA2"
**Expected!** These aren't in your enhancer library. The script will skip them and continue with the 11 available pairs.

To add them: Place FASTA files in `../../sequences/enhancers/` with format `TF_NAME.fasta`

### "AlphaGenome script not found"
Update the path in `run_logic_gate_predictions.py` or copy your prediction code from the regulatory_grammar experiment.

## 📝 Next Steps

1. **Run the build script** (fixed now!)
   ```bash
   python build_logic_gate_constructs.py
   ```

2. **Check the output**
   ```bash
   ls sequences/  # Should see ~44 FASTA files
   cat logic_gate_manifest.json | jq '.[0]'  # Inspect first construct
   ```

3. **Set up AlphaGenome** 
   - Update prediction script path
   - Or adapt your existing prediction pipeline

4. **Run pilot test** (optional - just 1 pair = 4 constructs)
   ```bash
   python run_logic_gate_predictions.py --limit 4
   ```

5. **Full run**
   ```bash
   ./run_full_pipeline.sh
   ```

## 📚 References

- Your regulatory grammar experiment: `../regulatory_grammar/`
- GATA1+KLF1 synergy: Your finding (1.26× additivity)
- TAL1+HNF4A interference: Your finding (0.66× additivity)
- Boolean logic in biology: Buchler et al., PNAS 2003

---

**Ready to discover if AlphaGenome learned Boolean logic! 🚀**

Questions? Check QUICKSTART.md or README.md for detailed documentation.
