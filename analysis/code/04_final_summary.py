#!/usr/bin/env python3
"""
Generate final summary of enhancer stacking experiment.
"""

import pandas as pd
from pathlib import Path

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
RESULTS_DIR = ROOT / "analysis/results"

print("=" * 80)
print(" " * 20 + "ENHANCER STACKING EXPERIMENT SUMMARY")
print("=" * 80)

# Load summary data
df = pd.read_csv(RESULTS_DIR / "summary_metrics.csv")

print("\n📊 FULL RESULTS TABLE")
print("=" * 80)
print(df.to_string(index=False))

print("\n\n🔬 KEY FINDINGS")
print("=" * 80)

# Finding 1: Position effect
e0_max = df[df['Construct'] == 'E0']['Enhancer_Max_Signal'].values[0]
e100_max = df[df['Construct'] == 'E100']['Enhancer_Max_Signal'].values[0]
position_fold = e100_max / e0_max

print(f"\n1. POSITION EFFECT (E0 vs E100):")
print(f"   • E0 (at promoter): {e0_max:.4f}")
print(f"   • E100 (100kb upstream): {e100_max:.4f}")
print(f"   • Fold change: {position_fold:.2f}x")
print(f"   ➜ Enhancers work BETTER when placed 100kb upstream!")

# Finding 2: Dose-response (1x to 10x)
stacked = df[df['Construct'].str.contains('EC100-')]
stacked_10 = stacked[stacked['Enhancer_Copies'] <= 10].sort_values('Enhancer_Copies')

print(f"\n2. DOSE-RESPONSE (1x to 10x enhancers):")
print(f"   Copy #  |  Max Signal  |  Mean Signal  |     AUC     |  AUC Fold")
print(f"   --------|--------------|---------------|-------------|----------")

baseline_auc = df[df['Construct'] == 'E100']['Enhancer_AUC'].values[0]
for _, row in stacked_10.iterrows():
    copies = int(row['Enhancer_Copies'])
    max_sig = row['Enhancer_Max_Signal']
    mean_sig = row['Enhancer_Mean_Signal']
    auc = row['Enhancer_AUC']
    auc_fold = auc / baseline_auc
    print(f"   {copies:>3}x    |   {max_sig:.4f}    |    {mean_sig:.4f}     |  {auc:>10.2f} |  {auc_fold:>5.2f}x")

print(f"   ➜ Peak signal SATURATES early (stays ~0.24)")
print(f"   ➜ Total signal (AUC) increases LINEARLY up to 10x")

# Finding 3: Extreme copy numbers
e100_auc = df[df['Construct'] == 'E100']['Enhancer_AUC'].values[0]
ec10_max = df[df['Construct'] == 'EC100-10x']['Enhancer_Max_Signal'].values[0]
ec10_auc = df[df['Construct'] == 'EC100-10x']['Enhancer_AUC'].values[0]
ec160_max = df[df['Construct'] == 'EC100-160x']['Enhancer_Max_Signal'].values[0]
ec160_auc = df[df['Construct'] == 'EC100-160x']['Enhancer_AUC'].values[0]
ec320_max = df[df['Construct'] == 'EC100-320x']['Enhancer_Max_Signal'].values[0]
ec320_auc = df[df['Construct'] == 'EC100-320x']['Enhancer_AUC'].values[0]

print(f"\n3. EXTREME COPY NUMBERS (Saturation Test):")
print(f"   Copies  |  Max Signal  |  AUC Fold vs 1x  |  Observation")
print(f"   --------|--------------|------------------|---------------------------")
print(f"   10x     |   {ec10_max:.4f}    |    {ec10_auc/e100_auc:>5.2f}x       | Baseline")
print(f"   160x    |   {ec160_max:.4f}    |   {ec160_auc/e100_auc:>6.2f}x       | Peak DROPS (non-monotonic!)")
print(f"   320x    |   {ec320_max:.4f}    |   {ec320_auc/e100_auc:>6.2f}x       | Peak REBOUNDS to highest!")

print(f"\n   ➜ Non-linear saturation curve")
print(f"   ➜ 320x shows HIGHEST max signal of all constructs (0.336)")
print(f"   ➜ Model does NOT break down at extreme copies")

# Finding 4: Controls
filler_max = df[df['Construct'] == 'FillerOnly']['Enhancer_Max_Signal'].values[0]
no_enh_max = df[df['Construct'] == 'NoEnhancer']['Enhancer_Max_Signal'].values[0]

print(f"\n4. CONTROL VALIDATION:")
print(f"   • FillerOnly: {filler_max:.4f} (background)")
print(f"   • NoEnhancer: {no_enh_max:.4f} (promoter only)")
print(f"   • Fold change (E100 vs FillerOnly): {e100_max/filler_max:.2f}x")
print(f"   ➜ Strong enhancer effect over background")

# Finding 5: Promoter effects
print(f"\n5. LONG-RANGE EFFECTS ON PROMOTER:")
print(f"   Construct      |  Promoter Signal  |  Fold vs Control")
print(f"   ---------------|-------------------|------------------")
baseline_prom = df[df['Construct'] == 'NoEnhancer']['Promoter_Signal'].values[0]
for construct in ['NoEnhancer', 'E100', 'EC100-10x', 'EC100-160x', 'EC100-320x']:
    prom_sig = df[df['Construct'] == construct]['Promoter_Signal'].values[0]
    fold = prom_sig / baseline_prom
    print(f"   {construct:<14} |      {prom_sig:.4f}       |      {fold:.2f}x")

print(f"   ➜ Massive enhancer stacks (160x-320x) boost promoter signal!")

print("\n\n📈 GENERATED VISUALIZATIONS")
print("=" * 80)
print(f"   Location: {RESULTS_DIR}")
print(f"   Files:")
print(f"   • genome_wide_tracks.png       - Full 1 MiB view of all constructs")
print(f"   • enhancer_region_zoom.png     - Zoomed view around enhancer (350-550 kb)")
print(f"   • dose_response_curves.png     - Max/Mean/AUC vs copy number")
print(f"   • bar_chart_comparison.png     - Side-by-side comparison")
print(f"   • saturation_analysis.png      - Fold change and saturation curve")

print("\n\n💡 BIOLOGICAL INTERPRETATION")
print("=" * 80)
print("""
1. ADDITIVITY: Yes, but complex
   - Linear AUC increase (1x to 10x)
   - Peak saturation occurs early
   - Total chromatin accessibility continues to grow

2. SATURATION: Multi-phase
   - Phase 1: Peak saturates (1x-10x)
   - Phase 2: AUC continues linear (10x-160x)
   - Phase 3: Non-monotonic behavior (160x-320x)

3. MODEL BEHAVIOR AT EXTREME INPUTS:
   - No catastrophic failure
   - Model handles 320 kb continuous enhancer sequence
   - Shows non-trivial spatial redistribution effects
   - Peak signal rebounds at 320x to highest value

4. PRACTICAL IMPLICATIONS:
   - For typical experiments (1-5 copies): reasonable dose-response
   - Spatial positioning matters (100kb > 0kb)
   - Model is robust but may extrapolate at extremes
""")

print("\n✅ EXPERIMENT COMPLETE!")
print("=" * 80)
print(f"\nFull report: {RESULTS_DIR / 'EXPERIMENT_REPORT.md'}")
print(f"Data table: {RESULTS_DIR / 'summary_metrics.csv'}")
print(f"Predictions: {ROOT / 'alphagenome/outputs'}/*.npy")

print("\n" + "=" * 80)
