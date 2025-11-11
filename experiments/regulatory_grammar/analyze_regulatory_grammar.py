#!/usr/bin/env python3
"""
Analyze regulatory grammar experiment results.

Generates:
1. Cell-type specificity heatmap
2. Cooperativity additivity scores
3. Spacing response curves
4. Orientation comparison plots
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Paths
EXPERIMENT_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/regulatory_grammar")
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
RESULTS_DIR = EXPERIMENT_ROOT / "results"

# Setup plotting
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def load_predictions(construct_name: str, cell_type: str = "K562") -> Dict:
    """Load predictions and stats for a construct."""
    base_name = f"{construct_name}_{cell_type}"
    
    # Load numpy array
    npy_path = OUTPUT_DIR / f"{base_name}_dnase.npy"
    if not npy_path.exists():
        return None
    
    predictions = np.load(npy_path)
    
    # Load stats
    stats_path = OUTPUT_DIR / f"{base_name}_stats.txt"
    stats = {}
    if stats_path.exists():
        with open(stats_path, 'r') as f:
            for line in f:
                if ':' in line and not line.startswith('#') and '=' not in line:
                    key, value = line.strip().split(':', 1)
                    try:
                        stats[key.strip()] = float(value.strip())
                    except:
                        pass
    
    return {
        "predictions": predictions,
        "stats": stats,
        "cell_type": cell_type,
    }


def analyze_celltype_specificity(manifest: List[Dict]) -> pd.DataFrame:
    """Analyze cell-type specificity experiments."""
    print("\n" + "=" * 80)
    print("ANALYSIS 1: Cell-Type Specificity")
    print("=" * 80)
    
    # Filter to cell-type experiments
    celltype_constructs = [c for c in manifest if c["experiment"] == "cell_type_specificity"]
    print(f"Found {len(celltype_constructs)} cell-type constructs")
    
    # Collect data
    data = []
    for construct in celltype_constructs:
        name = construct["construct"]
        result = load_predictions(name, construct["cell_type"])
        
        if result:
            data.append({
                "promoter": construct["promoter"],
                "enhancer": construct["enhancer"],
                "cell_type": construct["cell_type"],
                "expected": construct["expected"],
                "max_dnase": result["stats"].get("max", 0),
                "mean_dnase": result["stats"].get("mean", 0),
            })
    
    df = pd.DataFrame(data)
    
    # Create pivot table for heatmap
    for metric in ["max_dnase", "mean_dnase"]:
        pivot = df.pivot_table(
            values=metric,
            index=["promoter", "enhancer"],
            columns="cell_type",
            aggfunc="mean"
        )
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(pivot, annot=True, fmt=".4f", cmap="YlOrRd", ax=ax)
        ax.set_title(f"Cell-Type Specificity: {metric.replace('_', ' ').title()}")
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f"celltype_heatmap_{metric}.png")
        plt.close()
    
    # Save table
    df.to_csv(RESULTS_DIR / "celltype_specificity_results.csv", index=False)
    print(f"✓ Saved results to {RESULTS_DIR}")
    
    return df


def analyze_cooperativity(manifest: List[Dict]) -> pd.DataFrame:
    """Analyze pairwise cooperativity."""
    print("\n" + "=" * 80)
    print("ANALYSIS 2: Motif Cooperativity")
    print("=" * 80)
    
    # Get singles and pairs
    singles = {c["enhancer1"]: c for c in manifest 
               if c["experiment"] == "pairwise_cooperativity" and c["enhancer2"] is None}
    pairs = [c for c in manifest 
             if c["experiment"] == "pairwise_cooperativity" and c["enhancer2"] is not None]
    
    print(f"Found {len(singles)} single enhancers, {len(pairs)} pairs")
    
    # Load single enhancer signals
    single_signals = {}
    for enh_name, construct in singles.items():
        result = load_predictions(construct["construct"])
        if result:
            single_signals[enh_name] = result["stats"].get("max", 0)
    
    print("Single enhancer signals:")
    for enh, signal in single_signals.items():
        print(f"  {enh}: {signal:.4f}")
    
    # Calculate additivity scores
    data = []
    for construct in pairs:
        enh1 = construct["enhancer1"]
        enh2 = construct["enhancer2"]
        
        result = load_predictions(construct["construct"])
        if result and enh1 in single_signals and enh2 in single_signals:
            observed = result["stats"].get("max", 0)
            expected = single_signals[enh1] + single_signals[enh2]
            additivity = observed / expected if expected > 0 else 0
            
            data.append({
                "pair": f"{enh1}+{enh2}",
                "enhancer1": enh1,
                "enhancer2": enh2,
                "observed": observed,
                "expected": expected,
                "additivity_score": additivity,
                "interaction": "synergy" if additivity > 1.1 else ("interference" if additivity < 0.9 else "independent"),
                "expected_interaction": construct.get("expected", "unknown"),
            })
    
    df = pd.DataFrame(data)
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = df["interaction"].map({"synergy": "green", "interference": "red", "independent": "gray"})
    ax.barh(df["pair"], df["additivity_score"], color=colors)
    ax.axvline(1.0, color='black', linestyle='--', label='Additive (1.0)')
    ax.axvline(1.1, color='green', linestyle=':', alpha=0.5, label='Synergy threshold')
    ax.axvline(0.9, color='red', linestyle=':', alpha=0.5, label='Interference threshold')
    ax.set_xlabel("Additivity Score")
    ax.set_title("Motif Cooperativity: Observed / Expected Signal")
    ax.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "cooperativity_additivity_scores.png")
    plt.close()
    
    # Save table
    df.to_csv(RESULTS_DIR / "cooperativity_results.csv", index=False)
    print(f"✓ Saved results to {RESULTS_DIR}")
    
    return df


def analyze_spacing(manifest: List[Dict]) -> pd.DataFrame:
    """Analyze short-range spacing effects."""
    print("\n" + "=" * 80)
    print("ANALYSIS 3: Short-Range Spacing")
    print("=" * 80)
    
    # Filter to spacing experiments
    spacing_constructs = [c for c in manifest if c["experiment"] == "short_range_spacing"]
    spacing_constructs.sort(key=lambda x: x["spacing"])
    
    print(f"Found {len(spacing_constructs)} spacing constructs")
    
    # Collect data
    data = []
    for construct in spacing_constructs:
        result = load_predictions(construct["construct"])
        if result:
            data.append({
                "spacing": construct["spacing"],
                "max_dnase": result["stats"].get("max", 0),
                "mean_dnase": result["stats"].get("mean", 0),
            })
    
    df = pd.DataFrame(data)
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for ax, metric in zip(axes, ["max_dnase", "mean_dnase"]):
        ax.plot(df["spacing"], df[metric], 'o-', linewidth=2, markersize=8)
        ax.set_xlabel("Spacing between HS2 and GATA1 (bp)")
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.set_title(f"Spacing Effect on {metric.replace('_', ' ').title()}")
        ax.grid(True, alpha=0.3)
        
        # Mark potential optimal spacing
        optimal_idx = df[metric].idxmax()
        optimal_spacing = df.loc[optimal_idx, "spacing"]
        ax.axvline(optimal_spacing, color='red', linestyle='--', alpha=0.5, 
                   label=f'Optimal: {optimal_spacing} bp')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "spacing_response_curves.png")
    plt.close()
    
    # Save table
    df.to_csv(RESULTS_DIR / "spacing_results.csv", index=False)
    print(f"✓ Saved results to {RESULTS_DIR}")
    
    return df


def analyze_orientation(manifest: List[Dict]) -> pd.DataFrame:
    """Analyze orientation effects."""
    print("\n" + "=" * 80)
    print("ANALYSIS 4: Orientation Effects")
    print("=" * 80)
    
    # Filter to orientation experiments
    orient_constructs = [c for c in manifest if c["experiment"] == "orientation_effects"]
    
    print(f"Found {len(orient_constructs)} orientation constructs")
    
    # Collect data
    data = []
    for construct in orient_constructs:
        result = load_predictions(construct["construct"])
        if result:
            data.append({
                "pair": f"{construct['enhancer1']}+{construct['enhancer2']}",
                "orientation": f"({construct['orientation1']}, {construct['orientation2']})",
                "max_dnase": result["stats"].get("max", 0),
                "mean_dnase": result["stats"].get("mean", 0),
                "expected": construct.get("expected", "unknown"),
            })
    
    df = pd.DataFrame(data)
    
    # Plot grouped by pair
    pairs = df["pair"].unique()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, pair in enumerate(pairs[:4]):
        ax = axes[idx]
        pair_data = df[df["pair"] == pair]
        
        x = range(len(pair_data))
        ax.bar(x, pair_data["max_dnase"])
        ax.set_xticks(x)
        ax.set_xticklabels(pair_data["orientation"], rotation=45)
        ax.set_ylabel("Max DNase")
        ax.set_title(f"Orientation Effects: {pair}")
        ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "orientation_effects.png")
    plt.close()
    
    # Save table
    df.to_csv(RESULTS_DIR / "orientation_results.csv", index=False)
    print(f"✓ Saved results to {RESULTS_DIR}")
    
    return df


def main():
    print("=" * 80)
    print("Regulatory Grammar Analysis")
    print("=" * 80)
    
    # Setup
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load manifest
    with open(MANIFEST_PATH, 'r') as f:
        manifest = json.load(f)
    
    # Run analyses
    celltype_df = analyze_celltype_specificity(manifest)
    coop_df = analyze_cooperativity(manifest)
    spacing_df = analyze_spacing(manifest)
    orient_df = analyze_orientation(manifest)
    
    # Summary report
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Results saved to: {RESULTS_DIR}")
    print("\nGenerated files:")
    print("  - celltype_heatmap_max_dnase.png")
    print("  - celltype_heatmap_mean_dnase.png")
    print("  - cooperativity_additivity_scores.png")
    print("  - spacing_response_curves.png")
    print("  - orientation_effects.png")
    print("  - *.csv (data tables)")
    print("=" * 80)


if __name__ == "__main__":
    main()
