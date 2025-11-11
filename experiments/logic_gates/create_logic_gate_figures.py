#!/usr/bin/env python3
"""
Generate visualization figures for logic gate analysis.

Creates:
1. Truth table heatmaps for each TF pair
2. Logic score comparison plots
3. Confusion matrix (expected vs predicted gates)
4. Synergy scatter plots
5. Summary dashboard
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

EXPERIMENT_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/logic_gates")
ANALYSIS_PATH = EXPERIMENT_ROOT / "logic_gate_analysis.json"
RESULTS_DIR = EXPERIMENT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# Ideal logic gate patterns
IDEAL_GATES = {
    "AND": np.array([0, 0, 0, 1]),
    "OR": np.array([0, 1, 1, 1]),
    "NOT": np.array([1, 0, 1, 0]),
    "XOR": np.array([0, 1, 1, 0]),
}

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def load_results() -> List[Dict]:
    """Load analysis results."""
    with open(ANALYSIS_PATH, 'r') as f:
        return json.load(f)


def plot_truth_table_heatmap(result: Dict, ax: plt.Axes):
    """Plot truth table as 2x2 heatmap."""
    signals = result["signals_normalized"]
    
    # Reshape to 2x2 (TF_B × TF_A)
    truth_table = np.array(signals).reshape(2, 2)
    
    # Plot heatmap
    sns.heatmap(
        truth_table,
        annot=True,
        fmt=".2f",
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        cbar=False,
        ax=ax,
        square=True,
    )
    
    # Labels
    ax.set_xlabel(f"{result['tf_a']}", fontsize=10)
    ax.set_ylabel(f"{result['tf_b']}", fontsize=10)
    ax.set_xticklabels(["0", "1"])
    ax.set_yticklabels(["0", "1"])
    
    # Title with classification
    correct = "✓" if result["correct_classification"] else "✗"
    ax.set_title(
        f"{result['tf_a']}×{result['tf_b']}\n"
        f"Expected: {result['expected_gate']}, Got: {result['best_gate']} {correct}",
        fontsize=9
    )


def plot_logic_scores_bar(result: Dict, ax: plt.Axes):
    """Plot logic scores as bar chart."""
    scores = result["logic_scores"]
    expected = result["expected_gate"]
    best = result["best_gate"]
    
    gates = list(scores.keys())
    values = [scores[g] for g in gates]
    colors = ["#2ecc71" if g == expected else "#3498db" for g in gates]
    
    bars = ax.bar(gates, values, color=colors, alpha=0.7, edgecolor="black", linewidth=1)
    
    # Highlight best fit
    for i, gate in enumerate(gates):
        if gate == best:
            bars[i].set_edgecolor("red")
            bars[i].set_linewidth(2)
    
    ax.set_ylabel("Logic Score (R²)", fontsize=10)
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5, linewidth=1)
    ax.grid(axis="y", alpha=0.3)
    
    # Legend
    from matplotlib.patches import Rectangle
    legend_elements = [
        Rectangle((0,0),1,1, facecolor="#2ecc71", alpha=0.7, label="Expected"),
        Rectangle((0,0),1,1, facecolor="white", edgecolor="red", linewidth=2, label="Best fit"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="upper right")


def create_individual_pair_figure(result: Dict):
    """Create comprehensive figure for a single TF pair."""
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.3)
    
    # Truth table heatmap
    ax1 = fig.add_subplot(gs[0, 0])
    plot_truth_table_heatmap(result, ax1)
    
    # Logic scores
    ax2 = fig.add_subplot(gs[0, 1])
    plot_logic_scores_bar(result, ax2)
    
    # Ideal vs observed comparison
    ax3 = fig.add_subplot(gs[0, 2])
    plot_ideal_vs_observed(result, ax3)
    
    # Overall title
    fig.suptitle(
        f"Logic Gate Analysis: {result['tf_a']} × {result['tf_b']} ({result['cell_type']})",
        fontsize=14,
        fontweight="bold"
    )
    
    # Save
    filename = f"logic_gate_{result['tf_a']}_{result['tf_b']}_{result['cell_type']}.png"
    fig.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    return filename


def plot_ideal_vs_observed(result: Dict, ax: plt.Axes):
    """Plot ideal vs observed patterns."""
    observed = np.array(result["signals_normalized"])
    expected_gate = result["expected_gate"]
    ideal = IDEAL_GATES[expected_gate]
    
    x = np.arange(4)
    width = 0.35
    
    ax.bar(x - width/2, ideal, width, label="Ideal", color="#e74c3c", alpha=0.7)
    ax.bar(x + width/2, observed, width, label="Observed", color="#3498db", alpha=0.7)
    
    ax.set_ylabel("Normalized Signal", fontsize=10)
    ax.set_xlabel("Truth Table Condition", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(["00", "01", "10", "11"])
    ax.set_ylim(0, 1.2)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    
    # R² score
    score = result["logic_scores"][expected_gate]
    ax.text(
        0.98, 0.98, f"R² = {score:.3f}",
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        fontsize=10
    )


def create_all_truth_tables_grid(results: List[Dict]):
    """Create grid of all truth table heatmaps."""
    n = len(results)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    axes = axes.flatten() if n > 1 else [axes]
    
    for i, result in enumerate(results):
        plot_truth_table_heatmap(result, axes[i])
    
    # Hide empty subplots
    for i in range(n, len(axes)):
        axes[i].axis("off")
    
    fig.suptitle("Logic Gate Truth Tables - All TF Pairs", fontsize=16, fontweight="bold", y=0.995)
    fig.savefig(FIGURES_DIR / "all_truth_tables.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"  ✓ all_truth_tables.png")


def create_confusion_matrix(results: List[Dict]):
    """Create confusion matrix of expected vs predicted gate types."""
    expected = [r["expected_gate"] for r in results]
    predicted = [r["best_gate"] for r in results]
    
    gate_types = sorted(IDEAL_GATES.keys())
    confusion = np.zeros((len(gate_types), len(gate_types)), dtype=int)
    
    for exp, pred in zip(expected, predicted):
        i = gate_types.index(exp)
        j = gate_types.index(pred)
        confusion[i, j] += 1
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sns.heatmap(
        confusion,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=gate_types,
        yticklabels=gate_types,
        ax=ax,
        cbar_kws={"label": "Count"},
    )
    
    ax.set_xlabel("Predicted Gate Type", fontsize=12)
    ax.set_ylabel("Expected Gate Type", fontsize=12)
    ax.set_title("Logic Gate Classification Confusion Matrix", fontsize=14, fontweight="bold")
    
    # Compute accuracy
    accuracy = np.trace(confusion) / np.sum(confusion)
    ax.text(
        0.5, -0.15, f"Overall Accuracy: {accuracy:.1%}",
        transform=ax.transAxes,
        ha="center", fontsize=12, fontweight="bold"
    )
    
    fig.savefig(FIGURES_DIR / "confusion_matrix.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"  ✓ confusion_matrix.png")


def create_synergy_scatter(results: List[Dict]):
    """Create scatter plot of synergy metrics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Extract data
    expected_gates = [r["expected_gate"] for r in results]
    additivity = [r["synergy_metrics"]["additivity_ratio"] for r in results]
    eom = [r["synergy_metrics"]["excess_over_max"] for r in results]
    best_scores = [r["best_score"] for r in results]
    
    # Color by expected gate
    gate_colors = {"AND": "#2ecc71", "OR": "#3498db", "NOT": "#e74c3c", "XOR": "#9b59b6"}
    colors = [gate_colors[g] for g in expected_gates]
    
    # Plot 1: Additivity vs Logic Score
    ax1 = axes[0]
    scatter1 = ax1.scatter(additivity, best_scores, c=colors, s=100, alpha=0.7, edgecolors="black", linewidth=1)
    ax1.axvline(1.0, color="gray", linestyle="--", alpha=0.5, label="Additive (1.0)")
    ax1.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax1.set_xlabel("Additivity Ratio", fontsize=12)
    ax1.set_ylabel("Best Logic Score (R²)", fontsize=12)
    ax1.set_title("Synergy vs Logic Gate Performance", fontsize=13, fontweight="bold")
    ax1.grid(alpha=0.3)
    
    # Plot 2: Excess over Max vs Logic Score
    ax2 = axes[1]
    scatter2 = ax2.scatter(eom, best_scores, c=colors, s=100, alpha=0.7, edgecolors="black", linewidth=1)
    ax2.axvline(0, color="gray", linestyle="--", alpha=0.5, label="No excess")
    ax2.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Excess Over Maximum", fontsize=12)
    ax2.set_ylabel("Best Logic Score (R²)", fontsize=12)
    ax2.set_title("Cooperativity vs Logic Gate Performance", fontsize=13, fontweight="bold")
    ax2.grid(alpha=0.3)
    
    # Shared legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=gate, edgecolor="black") 
                      for gate, color in gate_colors.items()]
    fig.legend(handles=legend_elements, loc="upper center", ncol=4, 
              bbox_to_anchor=(0.5, 0.02), fontsize=11, frameon=True)
    
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(FIGURES_DIR / "synergy_analysis.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"  ✓ synergy_analysis.png")


def create_logic_scores_comparison(results: List[Dict]):
    """Create grouped bar chart comparing logic scores across all pairs."""
    # Organize data
    pairs = [f"{r['tf_a']}×{r['tf_b']}" for r in results]
    gate_types = list(IDEAL_GATES.keys())
    
    fig, ax = plt.subplots(figsize=(16, 6))
    
    x = np.arange(len(pairs))
    width = 0.2
    
    for i, gate in enumerate(gate_types):
        scores = [r["logic_scores"][gate] for r in results]
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, scores, width, label=gate, alpha=0.8)
        
        # Highlight expected gates
        for j, (score, result) in enumerate(zip(scores, results)):
            if result["expected_gate"] == gate:
                bars[j].set_edgecolor("red")
                bars[j].set_linewidth(2)
    
    ax.set_ylabel("Logic Score (R²)", fontsize=12)
    ax.set_xlabel("TF Pair", fontsize=12)
    ax.set_title("Logic Gate Scores Across All TF Pairs", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=45, ha="right", fontsize=9)
    ax.legend(title="Gate Type", fontsize=10)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5, linewidth=1)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, 1)
    
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "logic_scores_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"  ✓ logic_scores_comparison.png")


def create_summary_dashboard(results: List[Dict]):
    """Create comprehensive summary dashboard."""
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. Classification accuracy by gate type
    ax1 = fig.add_subplot(gs[0, 0])
    gate_types = sorted(IDEAL_GATES.keys())
    accuracies = []
    for gate in gate_types:
        gate_results = [r for r in results if r["expected_gate"] == gate]
        if gate_results:
            acc = sum(1 for r in gate_results if r["correct_classification"]) / len(gate_results)
            accuracies.append(acc * 100)
        else:
            accuracies.append(0)
    
    bars = ax1.bar(gate_types, accuracies, color=["#2ecc71", "#3498db", "#e74c3c", "#9b59b6"], alpha=0.7)
    ax1.axhline(50, color="gray", linestyle="--", alpha=0.5, label="Random (50%)")
    ax1.set_ylabel("Classification Accuracy (%)", fontsize=11)
    ax1.set_title("Accuracy by Expected Gate Type", fontsize=12, fontweight="bold")
    ax1.set_ylim(0, 100)
    ax1.legend(fontsize=9)
    
    # 2. Mean logic scores by gate type
    ax2 = fig.add_subplot(gs[0, 1])
    mean_scores = []
    for gate in gate_types:
        scores = [r["logic_scores"][gate] for r in results if r["expected_gate"] == gate]
        mean_scores.append(np.mean(scores) if scores else 0)
    
    ax2.bar(gate_types, mean_scores, color=["#2ecc71", "#3498db", "#e74c3c", "#9b59b6"], alpha=0.7)
    ax2.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax2.set_ylabel("Mean Logic Score (R²)", fontsize=11)
    ax2.set_title("Mean Logic Score for Expected Gate", fontsize=12, fontweight="bold")
    ax2.set_ylim(0, 1)
    
    # 3. Synergy distribution
    ax3 = fig.add_subplot(gs[0, 2])
    synergy_classes = ["synergy", "additive", "interference"]
    synergy_counts = [sum(1 for r in results if r["synergy_metrics"]["synergy_class"] == sc) 
                     for sc in synergy_classes]
    colors_syn = ["#2ecc71", "#f39c12", "#e74c3c"]
    ax3.pie(synergy_counts, labels=synergy_classes, autopct="%1.1f%%", 
           colors=colors_syn, startangle=90)
    ax3.set_title("Synergy Classification", fontsize=12, fontweight="bold")
    
    # 4-6. Example truth tables (best AND, OR, XOR)
    example_gates = ["AND", "OR", "XOR"]
    for i, gate in enumerate(example_gates):
        ax = fig.add_subplot(gs[1, i])
        gate_results = [r for r in results if r["expected_gate"] == gate]
        if gate_results:
            best_result = max(gate_results, key=lambda r: r["best_score"])
            plot_truth_table_heatmap(best_result, ax)
        else:
            ax.text(0.5, 0.5, f"No {gate} examples", ha="center", va="center")
            ax.axis("off")
    
    fig.suptitle("Logic Gate Experiment - Summary Dashboard", fontsize=16, fontweight="bold")
    fig.savefig(FIGURES_DIR / "summary_dashboard.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"  ✓ summary_dashboard.png")


def main():
    print("=" * 80)
    print("LOGIC GATE FIGURE GENERATION")
    print("=" * 80)
    
    # Create output directory
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load results
    results = load_results()
    print(f"\n✓ Loaded {len(results)} TF pair results\n")
    
    print("Generating figures...")
    
    # 1. Individual pair figures
    print("\n1. Individual TF pair figures:")
    for result in results:
        filename = create_individual_pair_figure(result)
        print(f"  ✓ {filename}")
    
    # 2. Grid of all truth tables
    print("\n2. Overview figures:")
    create_all_truth_tables_grid(results)
    
    # 3. Confusion matrix
    create_confusion_matrix(results)
    
    # 4. Synergy analysis
    create_synergy_scatter(results)
    
    # 5. Logic scores comparison
    create_logic_scores_comparison(results)
    
    # 6. Summary dashboard
    create_summary_dashboard(results)
    
    print(f"\n{'='*80}")
    print(f"✓ All figures saved to: {FIGURES_DIR}")
    print(f"{'='*80}\n")
    print(f"🎉 Complete! Review figures and prepare manuscript.")


if __name__ == "__main__":
    main()
