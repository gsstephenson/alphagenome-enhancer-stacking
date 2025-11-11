#!/usr/bin/env python3
"""
Create summary visualization comparing AlphaGenome predictions to biological expectations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Setup
RESULTS_DIR = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/regulatory_grammar/results")
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Load cooperativity data
df = pd.read_csv(RESULTS_DIR / "cooperativity_results.csv")

# Create comprehensive summary figure
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# 1. Additivity scores sorted
ax1 = fig.add_subplot(gs[0, :])
colors = df["interaction"].map({
    "synergy": "#2ecc71",      # Green
    "independent": "#95a5a6",   # Gray
    "interference": "#e74c3c"   # Red
})
bars = ax1.bar(range(len(df)), df["additivity_score"], color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
ax1.axhline(y=1.0, color='black', linestyle='--', linewidth=2, label='Perfect Additivity (1.0)')
ax1.axhline(y=1.1, color='green', linestyle=':', linewidth=1.5, alpha=0.7, label='Synergy Threshold (1.1)')
ax1.axhline(y=0.9, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label='Interference Threshold (0.9)')
ax1.set_xticks(range(len(df)))
ax1.set_xticklabels(df["pair"], rotation=45, ha='right')
ax1.set_ylabel("Additivity Score (Observed/Expected)", fontsize=12, fontweight='bold')
ax1.set_title("Enhancer Pair Cooperativity: AlphaGenome Predictions vs Biological Expectations", 
              fontsize=14, fontweight='bold')
ax1.legend(loc='upper right', fontsize=10)
ax1.grid(axis='y', alpha=0.3)

# Add value labels on bars
for i, (bar, row) in enumerate(zip(bars, df.itertuples())):
    height = bar.get_height()
    label = f"{height:.2f}"
    ax1.text(bar.get_x() + bar.get_width()/2., height,
            label, ha='center', va='bottom', fontsize=9, fontweight='bold')

# 2. Observed vs Expected scatter
ax2 = fig.add_subplot(gs[1, 0])
scatter_colors = df["interaction"].map({
    "synergy": "#2ecc71",
    "independent": "#95a5a6",
    "interference": "#e74c3c"
})
ax2.scatter(df["expected"], df["observed"], c=scatter_colors, s=200, alpha=0.7, edgecolors='black', linewidth=1.5)
max_val = max(df["expected"].max(), df["observed"].max()) * 1.1
ax2.plot([0, max_val], [0, max_val], 'k--', linewidth=2, label='Perfect Additivity')
ax2.set_xlabel("Expected (Sum of Singles)", fontsize=11, fontweight='bold')
ax2.set_ylabel("Observed (Pair)", fontsize=11, fontweight='bold')
ax2.set_title("Cooperativity: Observed vs Expected", fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(alpha=0.3)

# Add labels for key points
for _, row in df.iterrows():
    if row['additivity_score'] > 1.1 or row['additivity_score'] < 0.7:
        ax2.annotate(row['pair'], (row['expected'], row['observed']),
                    fontsize=8, ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

# 3. Interaction type pie chart
ax3 = fig.add_subplot(gs[1, 1])
interaction_counts = df["interaction"].value_counts()
colors_pie = ['#2ecc71', '#95a5a6', '#e74c3c']
explode = (0.1, 0, 0.1)  # Explode synergy and interference
wedges, texts, autotexts = ax3.pie(interaction_counts, labels=interaction_counts.index, autopct='%1.0f%%',
                                    colors=colors_pie, explode=explode, startangle=90,
                                    textprops={'fontsize': 11, 'fontweight': 'bold'})
ax3.set_title("Interaction Type Distribution", fontsize=12, fontweight='bold')

# 4. Biology vs Prediction mismatch
ax4 = fig.add_subplot(gs[2, :])
df_sorted = df.sort_values('additivity_score', ascending=True)
x_pos = range(len(df_sorted))

# Map expected interactions to numeric scores
expected_map = {"synergy": 1.2, "interference": 0.8, "unknown": 1.0}
df_sorted["expected_numeric"] = df_sorted["expected_interaction"].map(expected_map)

# Plot actual vs expected
ax4.scatter(x_pos, df_sorted["additivity_score"], s=200, label='AlphaGenome Predicted',
           color='blue', alpha=0.7, edgecolors='black', linewidth=1.5, marker='o', zorder=3)
ax4.scatter(x_pos, df_sorted["expected_numeric"], s=200, label='Biological Expectation',
           color='orange', alpha=0.7, edgecolors='black', linewidth=1.5, marker='s', zorder=3)

# Connect with lines to show discrepancy
for i, (idx, row) in enumerate(df_sorted.iterrows()):
    ax4.plot([i, i], [row['additivity_score'], row['expected_numeric']], 
            'k:', alpha=0.4, linewidth=1.5, zorder=1)

ax4.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.5, label='Additivity')
ax4.axhline(y=1.1, color='green', linestyle=':', linewidth=1.5, alpha=0.3)
ax4.axhline(y=0.9, color='red', linestyle=':', linewidth=1.5, alpha=0.3)
ax4.set_xticks(x_pos)
ax4.set_xticklabels(df_sorted["pair"], rotation=45, ha='right')
ax4.set_ylabel("Additivity Score", fontsize=11, fontweight='bold')
ax4.set_title("AlphaGenome Predictions vs Biological Expectations (Mismatches Highlighted)", 
             fontsize=12, fontweight='bold')
ax4.legend(fontsize=10, loc='upper left')
ax4.grid(axis='y', alpha=0.3)

# Highlight major mismatches
for i, (idx, row) in enumerate(df_sorted.iterrows()):
    diff = abs(row['additivity_score'] - row['expected_numeric'])
    if diff > 0.3:  # Large mismatch
        ax4.axvspan(i-0.4, i+0.4, alpha=0.15, color='red', zorder=0)

plt.suptitle("Regulatory Grammar Analysis: AlphaGenome's Understanding of Transcriptional Cooperativity",
            fontsize=16, fontweight='bold', y=0.995)

plt.savefig(RESULTS_DIR / "COMPREHENSIVE_SUMMARY.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved comprehensive summary figure: {RESULTS_DIR / 'COMPREHENSIVE_SUMMARY.png'}")

# Create spacing summary
spacing_df = pd.read_csv(RESULTS_DIR / "spacing_results.csv")

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Max DNase
ax1.plot(spacing_df["spacing"], spacing_df["max_dnase"], 'o-', linewidth=3, markersize=10, color='#3498db')
optimal_idx = spacing_df["max_dnase"].idxmax()
optimal_spacing = spacing_df.loc[optimal_idx, "spacing"]
optimal_signal = spacing_df.loc[optimal_idx, "max_dnase"]
ax1.axvline(optimal_spacing, color='red', linestyle='--', linewidth=2, alpha=0.7, 
           label=f'Optimal: {optimal_spacing} bp')
ax1.scatter([optimal_spacing], [optimal_signal], s=300, color='red', marker='*', 
           edgecolors='darkred', linewidth=2, zorder=5, label=f'Max Signal: {optimal_signal:.2f}')
ax1.set_xlabel("Distance Between HS2 and GATA1 (bp)", fontsize=11, fontweight='bold')
ax1.set_ylabel("Max DNase Signal", fontsize=11, fontweight='bold')
ax1.set_title("Spacing Effect on Peak Signal", fontsize=12, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# Relative to optimal
relative = (spacing_df["max_dnase"] / optimal_signal * 100)
colors = ['green' if x >= 90 else 'orange' if x >= 75 else 'red' for x in relative]
ax2.bar(spacing_df["spacing"], relative, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax2.axhline(100, color='red', linestyle='--', linewidth=2, label='100% (Optimal)')
ax2.axhline(90, color='orange', linestyle=':', linewidth=1.5, alpha=0.7, label='90% Threshold')
ax2.set_xlabel("Distance Between HS2 and GATA1 (bp)", fontsize=11, fontweight='bold')
ax2.set_ylabel("Relative Activity (%)", fontsize=11, fontweight='bold')
ax2.set_title("Activity Relative to Optimal Spacing", fontsize=12, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(axis='y', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(ax2.patches, relative)):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{val:.0f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.suptitle("Short-Range Spacing Analysis: HS2+GATA1 Distance Dependence",
            fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(RESULTS_DIR / "SPACING_SUMMARY.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved spacing summary figure: {RESULTS_DIR / 'SPACING_SUMMARY.png'}")

print("\n" + "="*80)
print("SUMMARY VISUALIZATIONS COMPLETE")
print("="*80)
print("\nKey Findings:")
print(f"  • Synergistic pairs: {(df['interaction'] == 'synergy').sum()}/10 (20%)")
print(f"  • Interference pairs: {(df['interaction'] == 'interference').sum()}/10 (50%)")
print(f"  • Mean additivity: {df['additivity_score'].mean():.3f}")
print(f"  • Optimal spacing: {optimal_spacing} bp")
print(f"  • Distance sensitivity: {(1 - spacing_df['max_dnase'].min()/optimal_signal)*100:.0f}% drop over 10kb")
print("\nCritical Mismatch:")
print("  ⚠️  HS2+GATA1: Predicted 0.89× (interference), Expected synergy!")
print("="*80)
