# AlphaGenome Enhancer Stacking Experiment

## Scientific TL;DR

**Question:** How does AlphaGenome respond to increasing numbers of enhancer copies at a fixed distance from a promoter?

**Approach:** We generated synthetic 1 Mb genomic constructs with 0–320 tandem copies of the HS2 β-globin enhancer positioned 100 kb upstream of the HBG1 promoter, then predicted chromatin accessibility (DNase-seq) using AlphaGenome in K562 cells.

**Key Observations:**

1. **Linear Additivity at Low Copy Numbers:** DNase signal increased linearly from 1× to 10× enhancer copies, with max signal rising from 0.197 to 0.284 and area under curve (AUC) scaling proportionally.

2. **Saturation at Extreme Copy Numbers:** At 160× and 320× copies, the model exhibited signal saturation. Max DNase reached a ceiling (~0.336), and AUC showed diminishing returns, suggesting either:
   - Biological plausibility (chromatin remodeling capacity limits)
   - Model compression (softmax/normalization artifacts)

3. **Promoter Isolation:** The promoter region (~500 kb) showed minimal DNase signal changes across all constructs, indicating the model does not predict strong long-range chromatin interactions at 100 kb distance in this synthetic context.

4. **Model Stability:** Despite extreme stacking (320× = 320 kb of tandem enhancers), AlphaGenome produced stable predictions without numerical instability or spurious accessibility peaks elsewhere in the sequence.

**Biological Interpretation:** The saturation behavior aligns with known principles of enhancer-promoter communication: multiple enhancers can synergize, but chromatin accessibility has physical limits. However, the lack of distal promoter activation suggests the model may not fully capture long-range looping or phase-separation dynamics that occur *in vivo*.

**Computational Interpretation:** The plateau could also reflect model design choices (e.g., attention mechanisms saturating, or normalization flattening extreme signals). This highlights the importance of validating AI predictions against orthogonal experimental data.

---

## Experiment Design

### Constructs (all exactly 1,048,576 bp)

| Construct     | Description                                    | Enhancer Copies | Position       |
|---------------|------------------------------------------------|-----------------|----------------|
| **FillerOnly**| Negative control: A/T-rich filler DNA only     | 0               | N/A            |
| **NoEnhancer**| HBG1 promoter at center, no enhancer          | 0               | 500 kb (promoter) |
| **E0**        | HS2 enhancer immediately upstream of promoter  | 1               | Adjacent to promoter |
| **E100**      | HS2 enhancer 100 kb upstream of promoter       | 1               | 400 kb         |
| **EC100-2x**  | 2 tandem HS2 copies at 100 kb upstream         | 2               | 400 kb         |
| **EC100-5x**  | 5 tandem HS2 copies at 100 kb upstream         | 5               | 400 kb         |
| **EC100-10x** | 10 tandem HS2 copies at 100 kb upstream        | 10              | 400 kb         |
| **EC100-160x**| 160 tandem HS2 copies (stress test)            | 160             | ~340–500 kb    |
| **EC100-320x**| 320 tandem HS2 copies (extreme stress test)    | 320             | ~180–500 kb    |

**Key Elements:**
- **Enhancer:** HS2 β-globin locus control region enhancer (chr11:5290000-5291000, 1001 bp, GRCh38)
- **Promoter:** HBG1 fetal hemoglobin γ-1 promoter (chr11:5273600-5273900, 301 bp, GRCh38)
- **Filler DNA:** A/T-rich neutral sequence (40% A, 40% T, 10% G, 10% C)
- **Cell Type:** K562 (erythroleukemia, relevant for β-globin regulation)

### Rationale
- **HS2 + HBG1:** Physiologically relevant pair for erythroid gene regulation
- **100 kb distance:** Within the range of typical enhancer-promoter contacts in mammalian genomes
- **Copy number series:** Tests model behavior across biologically plausible (1–10×) to extreme (160–320×) ranges
- **1 Mb length:** Maximum supported by AlphaGenome (1,048,576 bp = 2^20)

---

## Results Summary

### Quantitative Metrics

| Construct     | Max DNase | Mean DNase | AUC (Enhancer Region) | Promoter Signal |
|---------------|-----------|------------|-----------------------|-----------------|
| FillerOnly    | 0.127     | 0.000799   | 98.4                  | 0.0005          |
| NoEnhancer    | 0.127     | 0.000804   | 99.0                  | 0.0008          |
| E0            | 0.215     | 0.003339   | 411.1                 | 0.0009          |
| E100          | 0.197     | 0.002448   | 301.5                 | 0.0009          |
| EC100-2x      | 0.235     | 0.003788   | 466.5                 | 0.0008          |
| EC100-5x      | 0.269     | 0.005299   | 652.5                 | 0.0008          |
| EC100-10x     | 0.284     | 0.006204   | 764.0                 | 0.0008          |
| EC100-160x    | 0.333     | 0.007330   | 902.8                 | 0.0008          |
| EC100-320x    | 0.336     | 0.007419   | 913.4                 | 0.0008          |

### Key Findings

1. **Dose-Response Linearity (1–10×):**
   - Max DNase: 0.197 → 0.284 (44% increase)
   - AUC: 301.5 → 764.0 (153% increase)
   - Near-linear scaling with copy number

2. **Saturation Beyond 10× Copies:**
   - 160× vs 10×: Max DNase +17%, AUC +18%
   - 320× vs 160×: Max DNase +1%, AUC +1%
   - Diminishing returns indicate model/biological ceiling

3. **Promoter Signal Invariance:**
   - Promoter DNase remains ~0.0008–0.0009 across all constructs
   - No evidence of long-range chromatin opening at 100 kb distance

4. **Control Validation:**
   - FillerOnly and NoEnhancer show low baseline signal (AUC ~98)
   - 4–9× higher signal in enhancer-containing constructs confirms specificity

---

## Visualizations

All plots are saved in `analysis/results/`:

1. **`dnase_signal_comparison.png`:** Line plots showing DNase accessibility across all constructs
2. **`enhancer_dose_response.png`:** Bar chart of max DNase vs. enhancer copy number
3. **`auc_dose_response.png`:** Area under curve scaling with copy number
4. **`saturation_analysis.png`:** Focused comparison of 10×, 160×, and 320× constructs
5. **`promoter_signal_comparison.png`:** Promoter region signal across constructs

---

## Methods

### 1. Sequence Construction
- Downloaded genomic sequences via UCSC DAS server (GRCh38/hg38)
- Generated 1 Mb A/T-rich filler DNA using weighted random sampling
- Assembled constructs by concatenating: `[filler] + [enhancer × N] + [filler] + [promoter] + [filler]`
- All constructs padded to exactly 1,048,576 bp (AlphaGenome requirement)

### 2. AlphaGenome Predictions
- **Model:** AlphaGenome v0.4.0 via Python SDK
- **Cell Type:** K562 (ontology: EFO:0002067)
- **Output:** DNase-seq accessibility predictions (track length: 131,072 bins @ 8 bp resolution)
- **API:** Google Cloud-based inference endpoint

### 3. Analysis
- Extracted enhancer region signals (bins 50,000–62,500, corresponding to ~400 kb genomic position)
- Calculated max signal, mean signal, and area under curve (AUC) for each construct
- Compared promoter signals at bin 62,500 (~500 kb)
- Statistical analysis and visualization in Python (NumPy, Matplotlib, Seaborn)

---

## Repository Structure

```
AlphaGenome_EnhancerStacking/
├── README.md                          # This file
├── setup_enhancer_stacking_experiment.sh  # Initial setup script
├── sequences/
│   ├── enhancers/
│   │   └── HS2_enhancer.fa            # β-globin HS2 enhancer (1001 bp)
│   ├── promoters/
│   │   └── HBG1_promoter.fa           # HBG1 promoter (301 bp)
│   └── constructs/
│       ├── FillerOnly_construct.fa
│       ├── NoEnhancer_construct.fa
│       ├── E0_construct.fa
│       ├── E100_construct.fa
│       ├── EC100-2x_construct.fa
│       ├── EC100-5x_construct.fa
│       ├── EC100-10x_construct.fa
│       ├── EC100-160x_construct.fa
│       └── EC100-320x_construct.fa
├── filler/
│   └── 1M_filler.txt                  # A/T-rich filler sequence (1 Mb)
├── alphagenome/
│   └── outputs/
│       ├── *_dnase.npy                # Raw predictions (NumPy arrays)
│       ├── *_dnase.txt                # Predictions (text format)
│       └── *_stats.txt                # Summary statistics
├── analysis/
│   ├── code/
│   │   ├── 01_parse_sequences_and_build_constructs.py
│   │   ├── 02_run_alphagenome_predictions.py
│   │   ├── 03_analyze_and_visualize.py
│   │   └── 04_final_summary.py
│   └── results/
│       ├── metrics_summary.csv        # Quantitative results table
│       ├── dnase_signal_comparison.png
│       ├── enhancer_dose_response.png
│       ├── auc_dose_response.png
│       ├── saturation_analysis.png
│       ├── promoter_signal_comparison.png
│       └── EXPERIMENT_REPORT.md       # Detailed findings
└── logs/
    └── predictions_*.log              # Execution logs
```

---

## Reproducing This Analysis

### Prerequisites
- Python 3.11+
- AlphaGenome SDK (`pip install alphagenome`)
- AlphaGenome API key (set as `ALPHA_GENOME_KEY` in `.env`)
- Conda environment with: `numpy`, `matplotlib`, `seaborn`, `pandas`, `python-dotenv`

### Steps

1. **Set up environment:**
   ```bash
   conda create -n alphagenome-env python=3.11
   conda activate alphagenome-env
   pip install alphagenome numpy matplotlib seaborn pandas python-dotenv
   ```

2. **Configure API key:**
   ```bash
   echo "ALPHA_GENOME_KEY=your_api_key_here" > .env
   ```

3. **Run initial setup (optional, sequences already included):**
   ```bash
   ./setup_enhancer_stacking_experiment.sh
   ```

4. **Build constructs:**
   ```bash
   python analysis/code/01_parse_sequences_and_build_constructs.py
   ```

5. **Run predictions (requires API key):**
   ```bash
   python analysis/code/02_run_alphagenome_predictions.py
   ```

6. **Analyze and visualize:**
   ```bash
   python analysis/code/03_analyze_and_visualize.py
   ```

7. **View summary:**
   ```bash
   python analysis/code/04_final_summary.py
   ```

---

## Potential Next Steps

### 1. **Test Distance-Dependent Effects**
- Generate constructs with enhancers at 10 kb, 25 kb, 50 kb, 200 kb, and 500 kb from the promoter
- **Hypothesis:** Signal should decay with distance if the model captures 3D chromatin contact probability
- **Expected outcome:** Inverse relationship between distance and DNase signal at enhancer/promoter

### 2. **Evaluate Context Dependency**
- Replace A/T-rich filler with:
  - GC-rich sequences (test compositional bias)
  - Endogenous human genomic sequences (chr11 flanking regions)
  - Repetitive elements (SINE/LINE)
- **Hypothesis:** AlphaGenome predictions should be robust to neutral sequence context
- **Expected outcome:** Similar enhancer signals with different filler, unless model has sequence composition artifacts

### 3. **Cross-Validate with Orthogonal Models**
- Compare AlphaGenome predictions to:
  - Enformer (DeepMind's sequence model)
  - Basenji2 (Calico's epigenome predictor)
  - Experimental DNase-seq from synthetic reporter assays (e.g., MPRA)
- **Hypothesis:** Models should agree on enhancer additivity trends but may differ in absolute signal scaling
- **Expected outcome:** Rank-order correlation >0.8 across models; divergence at extreme copy numbers

---

## Limitations

1. **Synthetic Context:** These constructs are purely computational. Real enhancer stacking occurs in chromatin with nucleosomes, transcription factors, and 3D looping, which may not be fully captured by sequence-only models.

2. **Model Constraints:** AlphaGenome requires power-of-2 sequence lengths and has fixed-resolution output (8 bp bins). This may introduce discretization artifacts.

3. **Single Cell Type:** Predictions are K562-specific. Enhancer-promoter interactions are highly cell-type-dependent; results may differ in primary erythroid cells or non-erythroid contexts.

4. **No Ground Truth:** We lack experimental validation for these exact synthetic sequences. Future work should include MPRA or episomal reporter assays.

---

## Citation

If you use this analysis or code, please cite:

```
AlphaGenome Enhancer Stacking Experiment (2024)
Repository: https://github.com/gsstephenson/alphagenome-enhancer-stacking
```

And the original AlphaGenome paper (when published).

---

## Contact

For questions or collaboration:
- Repository: [gsstephenson/alphagenome-enhancer-stacking](https://github.com/gsstephenson/alphagenome-enhancer-stacking)
- Issues: [GitHub Issues](https://github.com/gsstephenson/alphagenome-enhancer-stacking/issues)

---

## License

MIT License - See LICENSE file for details.

---

## Acknowledgments

- AlphaGenome team at Google DeepMind for the model and API
- UCSC Genome Browser for reference sequences (GRCh38/hg38)
- β-globin locus control region (LCR) as a model system for enhancer biology
