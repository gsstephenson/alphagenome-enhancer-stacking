#!/usr/bin/env python3
"""
Analyze and visualize enhancer stacking experiment results.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path
import seaborn as sns

# Setup
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Define paths
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
DATA_DIR = ROOT / "alphagenome/outputs"
RESULTS_DIR = ROOT / "analysis/results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Construct information
CONSTRUCTS = {
    'FillerOnly': {'copies': 0, 'color': 'gray', 'label': 'Filler Only (Control)'},
    'NoEnhancer': {'copies': 0, 'color': 'lightgray', 'label': 'Promoter Only'},
    'E0': {'copies': 1, 'color': 'blue', 'label': 'E0 (1x, at promoter)'},
    'E100': {'copies': 1, 'color': 'green', 'label': 'E100 (1x, 100kb upstream)'},
    'EC100-2x': {'copies': 2, 'color': 'orange', 'label': 'EC100-2x'},
    'EC100-5x': {'copies': 5, 'color': 'red', 'label': 'EC100-5x'},
    'EC100-10x': {'copies': 10, 'color': 'purple', 'label': 'EC100-10x'},
    'EC100-160x': {'copies': 160, 'color': 'darkred', 'label': 'EC100-160x'},
    'EC100-320x': {'copies': 320, 'color': 'maroon', 'label': 'EC100-320x'},
}

# Genomic positions (adjusted for 1,048,576 bp constructs)
ENHANCER_POS = 400_000  # Enhancer region start
PROMOTER_POS = 500_000  # Promoter position
ENHANCER_SIZE = 1_001   # HS2 enhancer size

def load_predictions(construct_name):
    """Load DNase predictions for a construct."""
    npy_file = DATA_DIR / f"{construct_name}_dnase.npy"
    data = np.load(npy_file)
    # Flatten if needed
    if data.ndim > 1:
        data = data.flatten()
    return data

def extract_region(data, start, end):
    """Extract a specific genomic region from predictions."""
    return data[start:end]

def compute_metrics(data, enhancer_start, enhancer_end, promoter_pos):
    """Compute key metrics for a construct."""
    # Enhancer region metrics
    enhancer_region = data[enhancer_start:enhancer_end]
    enhancer_max = np.max(enhancer_region)
    enhancer_mean = np.mean(enhancer_region)
    enhancer_auc = np.trapz(enhancer_region)  # Area under curve
    
    # Promoter signal (window around promoter)
    promoter_window = 5000
    promoter_start = max(0, promoter_pos - promoter_window)
    promoter_end = min(len(data), promoter_pos + promoter_window)
    promoter_signal = np.mean(data[promoter_start:promoter_end])
    
    # Global metrics
    global_mean = np.mean(data)
    global_max = np.max(data)
    
    return {
        'enhancer_max': enhancer_max,
        'enhancer_mean': enhancer_mean,
        'enhancer_auc': enhancer_auc,
        'promoter_signal': promoter_signal,
        'global_mean': global_mean,
        'global_max': global_max,
    }

def main():
    print("=" * 70)
    print("ENHANCER STACKING ANALYSIS")
    print("=" * 70)
    
    # Load all predictions
    print("\nLoading predictions...")
    predictions = {}
    metrics = {}
    
    for name in CONSTRUCTS.keys():
        print(f"  Loading {name}...")
        predictions[name] = load_predictions(name)
        
        # Compute metrics
        # For EC100-160x and EC100-320x, enhancer region extends further
        if '160x' in name:
            enhancer_end = ENHANCER_POS + (160 * ENHANCER_SIZE)
        elif '320x' in name:
            enhancer_end = ENHANCER_POS + (320 * ENHANCER_SIZE)
        elif name == 'FillerOnly' or name == 'NoEnhancer':
            enhancer_end = ENHANCER_POS + 10_000  # Just a region to compare
        else:
            enhancer_end = ENHANCER_POS + 10_000  # 10kb window around enhancer
        
        metrics[name] = compute_metrics(predictions[name], ENHANCER_POS, enhancer_end, PROMOTER_POS)
    
    # Create comprehensive summary table
    print("\n" + "=" * 70)
    print("SUMMARY METRICS")
    print("=" * 70)
    
    summary_data = []
    for name, info in CONSTRUCTS.items():
        m = metrics[name]
        summary_data.append({
            'Construct': name,
            'Enhancer_Copies': info['copies'],
            'Enhancer_Max_Signal': m['enhancer_max'],
            'Enhancer_Mean_Signal': m['enhancer_mean'],
            'Enhancer_AUC': m['enhancer_auc'],
            'Promoter_Signal': m['promoter_signal'],
            'Global_Mean': m['global_mean'],
            'Global_Max': m['global_max'],
        })
    
    df = pd.DataFrame(summary_data)
    df = df.sort_values('Enhancer_Copies')
    
    print(df.to_string(index=False))
    
    # Save summary table
    csv_path = RESULTS_DIR / "summary_metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Saved summary table: {csv_path}")
    
    # ========================================================================
    # VISUALIZATION 1: Full genome-wide tracks
    # ========================================================================
    print("\n" + "=" * 70)
    print("Generating genome-wide track plot...")
    print("=" * 70)
    
    fig, axes = plt.subplots(len(CONSTRUCTS), 1, figsize=(16, 12), sharex=True, sharey=False)
    
    for idx, (name, info) in enumerate(CONSTRUCTS.items()):
        ax = axes[idx]
        data = predictions[name]
        x = np.arange(len(data)) / 1000  # Convert to kb
        
        ax.plot(x, data, color=info['color'], linewidth=0.5, alpha=0.8)
        ax.fill_between(x, data, color=info['color'], alpha=0.3)
        
        # Mark enhancer and promoter positions
        ax.axvline(ENHANCER_POS / 1000, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(PROMOTER_POS / 1000, color='blue', linestyle='--', alpha=0.5, linewidth=1)
        
        ax.set_ylabel('DNase', fontsize=8)
        ax.set_title(f"{info['label']} | Max: {metrics[name]['enhancer_max']:.4f}", 
                     fontsize=9, loc='left')
        ax.grid(True, alpha=0.3)
    
    axes[-1].set_xlabel('Position (kb)', fontsize=10)
    axes[0].set_title('AlphaGenome DNase Predictions - Full Genome View', 
                       fontsize=12, fontweight='bold', loc='center', pad=20)
    
    plt.tight_layout()
    fig_path = RESULTS_DIR / "genome_wide_tracks.png"
    plt.savefig(fig_path, bbox_inches='tight')
    print(f"✓ Saved: {fig_path}")
    plt.close()
    
    # ========================================================================
    # VISUALIZATION 2: Zoomed view around enhancer region
    # ========================================================================
    print("Generating zoomed enhancer region plot...")
    
    zoom_start = ENHANCER_POS - 50_000
    zoom_end = ENHANCER_POS + 150_000
    
    fig, axes = plt.subplots(len(CONSTRUCTS), 1, figsize=(14, 12), sharex=True)
    
    for idx, (name, info) in enumerate(CONSTRUCTS.items()):
        ax = axes[idx]
        data = predictions[name][zoom_start:zoom_end]
        x = np.arange(zoom_start, zoom_end) / 1000
        
        ax.plot(x, data, color=info['color'], linewidth=1.5, label=info['label'])
        ax.fill_between(x, data, color=info['color'], alpha=0.3)
        
        # Mark positions
        ax.axvline(ENHANCER_POS / 1000, color='red', linestyle='--', alpha=0.7, linewidth=1.5, label='Enhancer')
        ax.axvline(PROMOTER_POS / 1000, color='blue', linestyle='--', alpha=0.7, linewidth=1.5, label='Promoter')
        
        ax.set_ylabel('DNase Signal', fontsize=8)
        ax.set_title(f"{name} ({info['copies']} copies)", fontsize=9, loc='left')
        ax.legend(loc='upper right', fontsize=6)
        ax.grid(True, alpha=0.3)
    
    axes[-1].set_xlabel('Position (kb)', fontsize=10)
    axes[0].set_title('Zoomed View: Enhancer Region (350-550 kb)', 
                       fontsize=12, fontweight='bold', loc='center', pad=20)
    
    plt.tight_layout()
    fig_path = RESULTS_DIR / "enhancer_region_zoom.png"
    plt.savefig(fig_path, bbox_inches='tight')
    print(f"✓ Saved: {fig_path}")
    plt.close()
    
    # ========================================================================
    # VISUALIZATION 3: Enhancer copy number vs signal (dose-response)
    # ========================================================================
    print("Generating dose-response curves...")
    
    # Filter to only stacked enhancers (exclude controls and E0)
    stacked = df[df['Enhancer_Copies'].isin([1, 2, 5, 10, 160, 320])]
    stacked = stacked[~stacked['Construct'].isin(['E0'])]  # Exclude E0 (different position)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Max signal vs copy number
    ax = axes[0, 0]
    ax.plot(stacked['Enhancer_Copies'], stacked['Enhancer_Max_Signal'], 
            marker='o', markersize=8, linewidth=2, color='darkblue')
    ax.set_xlabel('Enhancer Copy Number', fontsize=11)
    ax.set_ylabel('Max DNase Signal', fontsize=11)
    ax.set_title('Max Signal vs Copy Number', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Mean signal vs copy number
    ax = axes[0, 1]
    ax.plot(stacked['Enhancer_Copies'], stacked['Enhancer_Mean_Signal'], 
            marker='s', markersize=8, linewidth=2, color='darkgreen')
    ax.set_xlabel('Enhancer Copy Number', fontsize=11)
    ax.set_ylabel('Mean DNase Signal', fontsize=11)
    ax.set_title('Mean Signal vs Copy Number', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # AUC vs copy number
    ax = axes[1, 0]
    ax.plot(stacked['Enhancer_Copies'], stacked['Enhancer_AUC'], 
            marker='^', markersize=8, linewidth=2, color='darkred')
    ax.set_xlabel('Enhancer Copy Number', fontsize=11)
    ax.set_ylabel('Area Under Curve (AUC)', fontsize=11)
    ax.set_title('Total Signal (AUC) vs Copy Number', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Promoter signal vs copy number
    ax = axes[1, 1]
    ax.plot(stacked['Enhancer_Copies'], stacked['Promoter_Signal'], 
            marker='D', markersize=8, linewidth=2, color='purple')
    ax.set_xlabel('Enhancer Copy Number', fontsize=11)
    ax.set_ylabel('Promoter Region Signal', fontsize=11)
    ax.set_title('Promoter Signal vs Copy Number', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    plt.suptitle('Dose-Response: Enhancer Stacking Effects', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    fig_path = RESULTS_DIR / "dose_response_curves.png"
    plt.savefig(fig_path, bbox_inches='tight')
    print(f"✓ Saved: {fig_path}")
    plt.close()
    
    # ========================================================================
    # VISUALIZATION 4: Bar chart comparison
    # ========================================================================
    print("Generating bar chart comparison...")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Max signal
    ax = axes[0]
    bars = ax.bar(range(len(df)), df['Enhancer_Max_Signal'], 
                   color=[CONSTRUCTS[name]['color'] for name in df['Construct']])
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Construct'], rotation=45, ha='right')
    ax.set_ylabel('Max DNase Signal', fontsize=11)
    ax.set_title('Max Signal by Construct', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Mean signal
    ax = axes[1]
    bars = ax.bar(range(len(df)), df['Enhancer_Mean_Signal'], 
                   color=[CONSTRUCTS[name]['color'] for name in df['Construct']])
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Construct'], rotation=45, ha='right')
    ax.set_ylabel('Mean DNase Signal', fontsize=11)
    ax.set_title('Mean Signal by Construct', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # AUC
    ax = axes[2]
    bars = ax.bar(range(len(df)), df['Enhancer_AUC'], 
                   color=[CONSTRUCTS[name]['color'] for name in df['Construct']])
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Construct'], rotation=45, ha='right')
    ax.set_ylabel('Area Under Curve', fontsize=11)
    ax.set_title('Total Signal (AUC) by Construct', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Enhancer Stacking: Signal Comparison', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig_path = RESULTS_DIR / "bar_chart_comparison.png"
    plt.savefig(fig_path, bbox_inches='tight')
    print(f"✓ Saved: {fig_path}")
    plt.close()
    
    # ========================================================================
    # VISUALIZATION 5: Saturation analysis
    # ========================================================================
    print("Generating saturation analysis...")
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot with error bars showing saturation
    stacked_sorted = stacked.sort_values('Enhancer_Copies')
    
    ax.plot(stacked_sorted['Enhancer_Copies'], stacked_sorted['Enhancer_Max_Signal'], 
            marker='o', markersize=10, linewidth=2.5, color='darkblue', label='Max Signal')
    
    # Add reference line for FillerOnly (background)
    filler_max = df[df['Construct'] == 'FillerOnly']['Enhancer_Max_Signal'].values[0]
    ax.axhline(filler_max, color='gray', linestyle=':', linewidth=2, alpha=0.7, label='Filler Background')
    
    # Calculate fold-change relative to 1x
    baseline = stacked_sorted[stacked_sorted['Enhancer_Copies'] == 1]['Enhancer_Max_Signal'].values[0]
    fold_changes = stacked_sorted['Enhancer_Max_Signal'] / baseline
    
    ax2 = ax.twinx()
    ax2.plot(stacked_sorted['Enhancer_Copies'], fold_changes, 
             marker='s', markersize=8, linewidth=2, color='red', alpha=0.7, label='Fold Change')
    ax2.set_ylabel('Fold Change vs 1x', fontsize=12, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    ax.set_xlabel('Enhancer Copy Number', fontsize=12)
    ax.set_ylabel('Max DNase Signal', fontsize=12)
    ax.set_title('Saturation Analysis: Max Signal vs Copy Number', fontsize=13, fontweight='bold')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)
    ax2.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    fig_path = RESULTS_DIR / "saturation_analysis.png"
    plt.savefig(fig_path, bbox_inches='tight')
    print(f"✓ Saved: {fig_path}")
    plt.close()
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE!")
    print("=" * 70)
    print(f"\nAll results saved to: {RESULTS_DIR}")
    print("\nGenerated files:")
    print("  • summary_metrics.csv")
    print("  • genome_wide_tracks.png")
    print("  • enhancer_region_zoom.png")
    print("  • dose_response_curves.png")
    print("  • bar_chart_comparison.png")
    print("  • saturation_analysis.png")

if __name__ == "__main__":
    main()
