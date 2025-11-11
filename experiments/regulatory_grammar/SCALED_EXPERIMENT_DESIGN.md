# Scaled Regulatory Grammar Experiment - Computational Design

**Date:** November 11, 2025  
**Investigator:** Grant Stephenson  
**Lab:** Layer Laboratory, CU Boulder  
**Timeline:** 4-6 weeks (rotation-compatible)  
**Budget:** $150-200 (computational only)

---

## 🎯 Objectives

1. **Expand TF pair library** from 10 → 30 pairs (covering multiple tissue types)
2. **Compare AI models** (AlphaGenome vs Enformer vs DeepSEA)
3. **Validate with ChIP-seq** co-occupancy data as ground truth
4. **Test PWM constructs** to separate sequence from chromatin effects
5. **Deliver publication-ready** computational manuscript

---

## 📊 Experimental Design Overview

### Current Status (Completed)
- ✅ **10 TF pairs** tested (erythroid + hepatic)
- ✅ **66 constructs** across 4 dimensions
- ✅ AlphaGenome predictions complete
- ✅ Key finding: 60% interference, HS2+GATA1 mismatch

### Proposed Expansion
- 🔄 **30 additional TF pairs** (tissue-specific)
- 🔄 **150-200 total constructs**
- 🔄 Multi-model predictions (Enformer, DeepSEA)
- 🔄 ChIP-seq correlation analysis
- 🔄 Synthetic PWM constructs

---

## 🧬 Part 1: TF Pair Selection (20-30 New Pairs)

### Selection Criteria
1. **Known cooperativity** from literature (positive controls)
2. **Tissue-specific** expression patterns (test cell-type specificity)
3. **Different TF families** (diversity in binding modes)
4. **ENCODE data availability** (ChIP-seq in multiple cell types)
5. **Clinical relevance** (disease-associated TFs)

### Proposed TF Pairs by Category

#### A. **Muscle/Cardiac (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **MyoD + MEF2A** | ✓ Synergy | Canonical muscle differentiation | ✓ C2C12 |
| **MyoD + E12** | ✓ Synergy | bHLH heterodimer | ✓ Myotubes |
| **GATA4 + NKX2-5** | ✓ Synergy | Cardiac specification | ✓ H9 |
| **SRF + MEF2C** | ✓ Synergy | Smooth muscle | ✓ HSMM |
| **TBX5 + NKX2-5** | ✓ Synergy | Heart development | ✓ Literature |

#### B. **Immune/Myeloid (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **PU.1 + C/EBPα** | ✓ Synergy | Myeloid differentiation | ✓ THP1, K562 |
| **PU.1 + IRF8** | ✓ Synergy | Dendritic cell fate | ✓ GM12878 |
| **NF-κB + AP-1** | ✓ Synergy | Inflammatory response | ✓ Multiple |
| **RUNX1 + GATA1** | Competitive | Erythroid-megakaryocyte switch | ✓ K562 |
| **PU.1 + GATA1** | ✗ Interference | Cross-lineage antagonism | ✓ K562 |

#### C. **Pluripotency/Stem Cell (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **Oct4 + Sox2** | ✓ Synergy | Pluripotency core | ✓ H1-hESC |
| **Oct4 + Nanog** | ✓ Synergy | Self-renewal | ✓ H1-hESC |
| **Sox2 + Brn2** | ✓ Synergy | Neural progenitors | ✓ Literature |
| **Klf4 + c-Myc** | ✓ Synergy | Reprogramming factors | ✓ Multiple |
| **Oct4 + GATA1** | ✗ Interference | Cross-lineage block | ✓ K562, H1 |

#### D. **Neural/Brain (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **NeuroD1 + Ascl1** | ✓ Synergy | Neurogenesis | ✓ Literature |
| **Pax6 + Sox2** | ✓ Synergy | Eye/brain development | ✓ Literature |
| **REST + CoREST** | ✓ Synergy | Neuronal gene repression | ✓ Multiple |
| **CREB + CBP** | ✓ Synergy | Activity-dependent transcription | ✓ Literature |
| **NeuroD1 + MyoD** | ✗ Interference | Cross-lineage | ✓ C2C12 |

#### E. **Liver/Metabolic (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **HNF4A + HNF1A** | ✓ Synergy | Hepatocyte function | ✓ HepG2 |
| **HNF4A + C/EBPα** | ✓ Synergy | Metabolic genes | ✓ HepG2 |
| **PPARγ + RXR** | ✓ Synergy | Adipogenesis | ✓ Literature |
| **SREBP + ChREBP** | ✓ Synergy | Lipid metabolism | ✓ Literature |
| **HNF4A + MyoD** | ✗ Interference | Cross-lineage | ✓ HepG2, C2C12 |

#### F. **Additional Cross-Context Tests (5 pairs)**
| Pair | Expected Interaction | Rationale | ENCODE Data |
|------|---------------------|-----------|-------------|
| **CTCF + Cohesin** | Structural | Insulation/looping | ✓ Multiple |
| **p53 + p63** | ✓ Synergy | DNA damage response | ✓ Multiple |
| **STAT3 + NF-κB** | ✓ Synergy | Cytokine signaling | ✓ Multiple |
| **AR + FOXA1** | ✓ Synergy | Prostate-specific | ✓ Literature |
| **ER + FOXA1** | ✓ Synergy | Breast-specific | ✓ MCF7 |

---

## 📋 Enhancer Sequences Needed (BED Format)

### What I Need From You:

For each TF listed above, provide **ONE representative enhancer** in BED format with these columns:

```
chr    start    end    name    score    strand
```

**Example format:**
```bed
chr11    5290000    5291001    HS2_enhancer    1000    +
chrX     48649200   48650321   GATA1_enhancer  900     +
chr19    45420000   45420552   KLF1_enhancer   850     +
```

### Enhancer Selection Criteria:
1. **Size:** 500-2000 bp (optimal for AlphaGenome)
2. **Activity:** Strong ENCODE DNase/H3K27ac signal in relevant cell type
3. **Annotation:** Known TF binding from ChIP-seq
4. **Quality:** Avoid repetitive elements (>20% RepeatMasker)
5. **Accessibility:** Public genome assembly (GRCh38/hg38)

### Recommended Sources:

**Option 1: ENCODE Candidate Cis-Regulatory Elements (cCREs)**
- URL: https://screen.encodeproject.org/
- Filter: "Proximal enhancer-like" or "Distal enhancer-like"
- Select: High DNase + H3K27ac signal
- Download: BED file per TF

**Option 2: EPDnew Enhancer Database**
- URL: https://epd.epfl.ch/
- Pre-curated enhancers with TF binding annotations

**Option 3: VISTA Enhancer Browser**
- URL: https://enhancer.lbl.gov/
- Experimentally validated enhancers (gold standard)

**Option 4: Literature Mining**
- PubMed search: "[TF name] enhancer ChIP-seq"
- Extract coordinates from supplementary tables
- Verify with UCSC Genome Browser

---

## 🔧 Part 2: Construct Design Strategy

### A. **Pairwise Cooperativity (Primary Analysis)**

**Design:** 30 new pairs + singles = 91 constructs

```
[Filler] - [Enhancer1] - [5kb spacing] - [Enhancer2] - [95kb filler] - [HBG1 Promoter] - [Filler]
```

**Controls:**
- 30 × Single enhancer (baseline)
- 30 × Paired enhancers (test cooperativity)
- 30 × CTCF-separated pairs (insulation test)
- 1 × Filler only (negative control)

**Total:** 91 new constructs

### B. **Cell-Type Specificity Matrix (Secondary)**

Test key pairs in 3 cell types: K562, HepG2, GM12878

**Example:**
- MyoD + MEF2A in K562 (wrong context)
- MyoD + MEF2A in HepG2 (wrong context)
- PU.1 + C/EBPα in K562 (semi-relevant)
- PU.1 + C/EBPα in GM12878 (correct context)

**Total:** 15-20 constructs (select most interesting pairs)

### C. **Synthetic PWM Constructs (Mechanistic)**

Replace natural enhancers with isolated JASPAR PWM motifs

**Design:**
```
[Filler] - [6x GATA motif] - [5kb] - [6x KLF motif] - [95kb] - [Promoter] - [Filler]
```

**Compare:**
- Natural enhancers (full chromatin context)
- PWM arrays (pure motif cooperativity)
- Scrambled sequences (negative control)

**Total:** 20-30 constructs

### D. **Spacing Variants (Extended)**

Test optimal spacing for top 5 synergistic pairs at:
- 0 bp, 500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 20 kb

**Total:** 35 constructs (5 pairs × 7 distances)

---

## 💻 Part 3: Computational Pipeline

### Step 1: Build Constructs (Week 1)

**Script:** `build_scaled_regulatory_grammar.py`

```python
# Pseudocode structure
1. Load TF enhancer BED files (you provide)
2. Extract sequences from GRCh38 reference
3. Generate all construct combinations:
   - Pairwise cooperativity (91)
   - Cell-type matrix (20)
   - PWM constructs (30)
   - Spacing variants (35)
4. Write FASTA files (176 total)
5. Generate metadata JSON
```

**Output:** 176 FASTA files, 1 manifest JSON

**Time:** 1-2 days (mostly automated)

---

### Step 2: AlphaGenome Predictions (Week 1-2)

**Script:** `run_alphagenome_scaled_predictions.py`

**API calls:** 176 constructs × $0.10 = **$17.60**

**Time:** 6-8 hours (batched)

**Output:** 176 × DNase prediction arrays (.npy files)

---

### Step 3: Enformer Predictions (Week 2-3)

**Two options:**

#### Option A: Enformer via Kipoi (Local)
```bash
conda install -c bioconda kipoi
kipoi predict Enformer/human --source=kipoi
```
- Free but requires GPU (24GB VRAM)
- Runtime: ~8 hours for 176 sequences

#### Option B: Enformer via Google Colab
- Use Deepmind's official notebook
- Free tier: ~20-30 sequences/day
- Total time: 1 week (parallelizable)

**Output format:** Same as AlphaGenome (chromatin tracks)

---

### Step 4: DeepSEA Predictions (Week 2-3)

**Method:** DeepSEA web server or Selene library

```bash
pip install selene-sdk
```

**Runtime:** ~4 hours for 176 sequences

**Note:** DeepSEA has shorter sequence limit (1000 bp) - may need to extract just enhancer regions

---

### Step 5: ChIP-seq Co-Occupancy Analysis (Week 3-4)

**Data sources:**
1. **ENCODE ChIP-seq:** Download processed peaks (BED files)
2. **ReMap 2022:** Pre-computed TF binding atlas
3. **ChIP-Atlas:** Meta-analysis of public ChIP-seq

**Analysis pipeline:**

```python
# Pseudocode
1. Download ChIP-seq peaks for all 40 TFs
2. For each TF pair:
   - Find peaks for TF1 and TF2
   - Calculate co-binding frequency (peaks within 1 kb)
   - Compute enrichment vs random expectation
3. Correlate co-binding with predicted cooperativity
4. Generate scatter plots (observed vs expected)
```

**Metric:** Co-binding enrichment score

**Expected correlation:** If models learn cooperativity correctly, high co-binding should correlate with high additivity scores

---

### Step 6: Statistical Analysis (Week 4)

**Analyses:**

1. **Cooperativity Distribution**
   - Histogram of additivity scores
   - Compare across TF families
   - AlphaGenome vs Enformer vs DeepSEA

2. **Model Agreement**
   - Pairwise correlation between models
   - Concordance on synergy/interference calls
   - Bland-Altman plots

3. **ChIP-seq Validation**
   - Pearson/Spearman correlation
   - ROC curves (predict synergy from co-binding)
   - Confusion matrix

4. **Cell-Type Specificity**
   - ANOVA across cell types
   - Effect size (Cohen's d)
   - Tissue-specific enrichment

5. **PWM Analysis**
   - Natural vs synthetic cooperativity
   - Contribution of chromatin context
   - Motif spacing optimization

**Tools:** Python (scipy, statsmodels), R (ggplot2, ComplexHeatmap)

---

### Step 7: Visualization (Week 4-5)

**Figure 1: Expanded Cooperativity Landscape**
- Heatmap of 40 TF pairs
- AlphaGenome vs Enformer side-by-side
- Color-coded by tissue type
- **Message:** Systematic interference across all models

**Figure 2: Model Comparison**
- 3-way scatter plots (AlphaGenome vs Enformer vs DeepSEA)
- Regression lines + R²
- Outliers highlighted (disagreement cases)
- **Message:** Model consensus or divergence

**Figure 3: ChIP-seq Validation**
- Observed cooperativity (y-axis) vs ChIP-seq co-binding (x-axis)
- Each point = one TF pair
- GATA1+KLF1, HS2+GATA1 labeled
- **Message:** Biological validation of predictions

**Figure 4: PWM vs Natural Enhancers**
- Bar chart comparison for 10 pairs
- Natural, PWM, scrambled side-by-side
- **Message:** Chromatin context contribution

**Figure 5: Cell-Type Specificity**
- Clustered heatmap of 20 pairs × 3 cell types
- Expected vs observed patterns annotated
- **Message:** Limited cell-type differentiation

**Figure 6: Tissue-Specific Insights**
- Box plots by tissue (muscle, immune, liver, neural)
- Compare cooperativity patterns
- **Message:** Tissue-level regulatory logic

---

## 📊 Part 4: ChIP-seq Data Requirements

### Datasets to Download (Public, Free)

#### A. **ENCODE ChIP-seq** (Primary Source)
**URL:** https://www.encodeproject.org/

**For each TF, download:**
```
1. Optimal IDR thresholded peaks (BED narrowPeak)
2. Fold change signal track (bigWig)
3. Cell type: K562 (erythroid)
4. Biological replicates: ≥2
```

**Example search:**
```
Assay: TF ChIP-seq
Target: GATA1-human
Biosample: K562
File format: bed narrowPeak
Status: released
```

**Files needed:** ~40 TF × 1-3 cell types = 120 BED files

**Storage:** ~5-10 GB total

---

#### B. **ReMap 2022** (Alternative/Supplement)
**URL:** https://remap2022.univ-amu.fr/

**Advantage:** Pre-computed TF binding catalog
**Download:** All TFs in one merged BED file
**Filter:** By cell type and TF name

---

#### C. **ChIP-Atlas** (Meta-Analysis)
**URL:** https://chip-atlas.org/

**Advantage:** Quality-filtered peaks across studies
**Use:** Validation and cross-study consensus

---

### Data Processing Pipeline

```python
# pseudocode: analyze_chipseq_cooperativity.py

import pyBigWig
import pybedtools

def calculate_cobinding(tf1_peaks, tf2_peaks, window=1000):
    """
    Calculate co-occupancy frequency
    """
    # Extend peaks by window
    tf1_extended = tf1_peaks.slop(b=window)
    
    # Find overlaps
    overlaps = tf1_extended.intersect(tf2_peaks, c=True)
    
    # Calculate enrichment
    observed = overlaps.count()
    expected = random_expectation(tf1_peaks, tf2_peaks)
    enrichment = observed / expected
    
    return enrichment

# For each TF pair:
for pair in tf_pairs:
    # Load ChIP-seq peaks
    tf1_peaks = load_bed(f"{pair.tf1}_K562_peaks.bed")
    tf2_peaks = load_bed(f"{pair.tf2}_K562_peaks.bed")
    
    # Calculate co-binding
    cobinding = calculate_cobinding(tf1_peaks, tf2_peaks)
    
    # Correlate with AlphaGenome prediction
    alphagenome_score = get_additivity_score(pair)
    
    correlations.append((cobinding, alphagenome_score))

# Statistical test
rho, pval = spearmanr(cobinding_scores, alphagenome_scores)
print(f"Correlation: {rho:.3f}, p={pval:.3e}")
```

---

## 🧬 Part 5: Synthetic PWM Design

### Strategy: Isolate Sequence-Level Grammar

**Hypothesis:** If cooperativity persists with pure PWM motifs, it's sequence-driven. If lost, chromatin context matters.

### JASPAR Motif Selection

**For each TF, download PWM from JASPAR:**
- URL: https://jaspar.genereg.net/
- Format: MEME or PFM
- Select: Highest information content matrix

**Example TFs:**
- GATA1: MA0035.4
- KLF1: MA0493.1
- HNF4A: MA0114.4
- MyoD: MA0499.1
- PU.1: MA0080.5

### Construct Design

**Single PWM construct:**
```
[Filler] - [6× GATA motif, 20bp spacing] - [5kb filler] - [Promoter] - [Filler]
```

**Paired PWM construct:**
```
[Filler] - [6× GATA, 20bp spacing] - [5kb] - [6× KLF, 20bp spacing] - [95kb] - [Promoter] - [Filler]
```

**Controls:**
- Scrambled motif (preserves dinucleotide composition)
- Reverse complement (tests orientation)
- Single motif (tests redundancy)

**Number of copies:** 6× per motif (matches typical enhancer density)

**Spacing within motif array:** 20 bp (allows TF binding without overlap)

### Expected Results

| Condition | Natural Enhancer | PWM Array | Interpretation |
|-----------|-----------------|-----------|----------------|
| GATA1 alone | 1.62 | ? | Baseline motif strength |
| KLF1 alone | 4.78 | ? | Baseline motif strength |
| GATA1+KLF1 | 8.06 (1.26×) | ? | Test if synergy persists |

**If PWM synergy = natural synergy:** Sequence-level cooperativity  
**If PWM synergy < natural synergy:** Chromatin context required  
**If PWM synergy > natural synergy:** Natural enhancers have inhibitory elements

---

## 📦 Deliverables & Timeline

### Week 1: Setup & Construction
- [ ] Receive enhancer BED files from you
- [ ] Extract sequences from GRCh38
- [ ] Build 176 synthetic constructs
- [ ] Validate FASTA files
- [ ] Generate metadata JSON

### Week 2: Predictions (AlphaGenome)
- [ ] Run AlphaGenome API calls (176 constructs)
- [ ] Parse outputs
- [ ] Calculate cooperativity scores
- [ ] Quick QC analysis

### Week 3: Predictions (Enformer/DeepSEA)
- [ ] Set up Enformer pipeline
- [ ] Run Enformer predictions
- [ ] Run DeepSEA predictions
- [ ] Harmonize output formats

### Week 4: ChIP-seq Analysis
- [ ] Download ENCODE ChIP-seq data
- [ ] Process peaks (BEDTools)
- [ ] Calculate co-binding enrichment
- [ ] Correlate with predictions

### Week 5: Statistics & Visualization
- [ ] Comprehensive statistical tests
- [ ] Generate 6 main figures
- [ ] Supplementary tables
- [ ] Draft results section

### Week 6: Manuscript & Documentation
- [ ] Write computational methods
- [ ] Draft discussion
- [ ] Prepare for publication
- [ ] Handoff documentation for wet lab validation

---

## 💰 Budget Breakdown

| Item | Cost | Notes |
|------|------|-------|
| AlphaGenome API (176 constructs) | $17.60 | $0.10 per prediction |
| Enformer (free, local GPU) | $0 | Requires CUDA GPU |
| OR Enformer (Google Cloud) | $20-30 | If GPU unavailable |
| DeepSEA (free, local) | $0 | CPU-based |
| ChIP-seq data | $0 | Public ENCODE data |
| Compute resources | $20-50 | AWS/GCP for heavy analysis |
| **Total** | **$37-97** | Well under $200 budget |

---

## 📊 Expected Outcomes

### Quantitative Metrics
- **40 TF pairs tested** (4× current)
- **176 constructs** (2.7× current)
- **3 AI models compared** (AlphaGenome, Enformer, DeepSEA)
- **40 TF ChIP-seq datasets** analyzed
- **~50 pages** of supplementary data

### Key Questions Answered
1. ✅ Is sub-additive cooperativity general across TF families?
2. ✅ Do all AI models show same patterns? (model comparison)
3. ✅ Does ChIP-seq co-binding predict cooperativity? (validation)
4. ✅ Is cooperativity sequence-driven or chromatin-driven? (PWM analysis)
5. ✅ Can models differentiate cell-type contexts? (specificity)

### Publication Potential
- **Target journals:** Genome Biology, Nucleic Acids Research, Bioinformatics
- **Expected impact:** High (systematic AI model evaluation)
- **Novelty:** First comprehensive cooperativity benchmark for genomic AI
- **Utility:** Community resource for future model development

---

## 🔬 Wet Lab Validation Protocol (For Future)

### Luciferase Reporter Assay Design

**Priority Pairs to Validate:**
1. **GATA1 + KLF1** (predicted synergy 1.26×)
2. **HS2 + GATA1** (predicted interference 0.89×, biology expects synergy)
3. **MyoD + MEF2A** (literature synergy)
4. **PU.1 + GATA1** (cross-lineage interference)

**Constructs:**
```
pGL4-Luc vector:
- Minimal promoter (50 bp)
- Enhancer 1 (natural sequence, 500-1000 bp)
- Enhancer 2 (natural sequence, 500-1000 bp)
- Luciferase reporter
```

**Controls:**
- Empty vector (background)
- Single enhancers (baseline)
- Paired enhancers (test cooperativity)
- Scrambled sequences (specificity)

**Cell lines:**
- K562 (erythroid)
- HepG2 (hepatic)
- C2C12 (myoblast)

**Readout:** Firefly luciferase / Renilla (normalization)

**Statistics:** n=6 biological replicates, t-test, ANOVA

**Timeline:** 3-4 months

**Cost:** ~$5-10K (cloning, reagents, cell culture)

---

## 📋 What I Need From You RIGHT NOW

### Priority 1: TF Enhancer Coordinates (BED files)

**For these 30 TFs, provide one enhancer each:**

#### Must-Have (20 TFs):
1. MyoD (muscle)
2. MEF2A (muscle)
3. GATA4 (cardiac)
4. NKX2-5 (cardiac)
5. PU.1 (myeloid)
6. C/EBPα (myeloid)
7. IRF8 (immune)
8. NF-κB (immune)
9. AP-1 (immune)
10. Oct4 (pluripotency)
11. Sox2 (pluripotency)
12. Nanog (pluripotency)
13. Klf4 (reprogramming)
14. NeuroD1 (neural)
15. Pax6 (neural)
16. HNF1A (liver)
17. PPARγ (adipose)
18. p53 (DNA damage)
19. STAT3 (signaling)
20. REST (repressor)

#### Nice-to-Have (10 TFs):
21. E12 (bHLH)
22. SRF (muscle)
23. TBX5 (cardiac)
24. RUNX1 (hematopoietic)
25. Brn2 (neural)
26. c-Myc (oncogene)
27. Ascl1 (neural)
28. CoREST (repressor)
29. CREB (signaling)
30. FOXA1 (pioneer factor)

**Format:**
```bed
chr    start    end    name              score    strand
chr1   1000000  1001500  MyoD_enhancer_1  900      +
```

### Priority 2: Confirm Cell Types for AlphaGenome

Which cell types are supported by AlphaGenome API?
- K562 ✓ (confirmed)
- HepG2 ✓ (confirmed)
- GM12878 ✓ (confirmed)
- C2C12 (myoblast)? 
- H1-hESC (stem cell)?
- THP-1 (monocyte)?

### Priority 3: Computational Resources

Do you have access to:
- [ ] GPU cluster (for Enformer)
- [ ] AWS/GCP credits
- [ ] High-memory nodes (for large BED files)

---

## 🚀 Next Steps

**Once you provide enhancer coordinates:**

1. I'll generate construct FASTA files (1-2 days)
2. Run AlphaGenome predictions (1 day)
3. Analyze cooperativity scores (1 day)
4. Share preliminary results for feedback

**Total turnaround:** 4-5 days after receiving BED files

---

## 📧 Data Sharing Format

**Please send via:**
1. GitHub repo (preferred)
2. Dropbox/Google Drive link
3. Email attachment (if <10 MB)

**Required files:**
- `enhancer_coordinates.bed` (30 enhancers)
- `enhancer_metadata.csv` (TF name, cell type, source)

**Optional but helpful:**
- ChIP-seq BED files (if already downloaded)
- Literature references for cooperativity
- Known synergy/interference examples

---

## 💡 Questions for You

Before I start building:

1. **What's your rotation end date?** (to prioritize deliverables)
2. **Do you have GPU access** for Enformer, or should I use cloud?
3. **Which 10-15 TF pairs are highest priority** if we need to cut scope?
4. **Do you want PWM analysis included?** (adds 30 constructs)
5. **Cell-type matrix: which pairs to test in 3 cell types?** (adds 15-20 constructs)

---

## 📚 References for TF Pairs

I'll compile these as I build constructs:

- MyoD+MEF2: Molkentin et al. Cell 1995
- PU.1+C/EBP: Zhang et al. Nature 1996  
- Oct4+Sox2: Boyer et al. Cell 2005
- GATA1+KLF1: Drissen et al. EMBO J 2004
- NF-κB+AP-1: Stein et al. Mol Cell 1993

(Full bibliography in final manuscript)

---

**Ready to start as soon as you send the BED files!** 🚀

Let me know if you need help identifying enhancer coordinates or have questions about the design.
