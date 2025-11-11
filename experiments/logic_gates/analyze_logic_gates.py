#!/usr/bin/env python3
"""
Analyze logic gate experiments - compute truth tables and logic scores.

For each TF pair, computes:
1. Truth table (4 conditions: 00, 01, 10, 11)
2. Logic score for each gate type (AND, OR, NOT, XOR)
3. Best-fit gate classification
4. Synergy/additivity metrics
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error

EXPERIMENT_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/logic_gates")
MANIFEST_PATH = EXPERIMENT_ROOT / "logic_gate_manifest.json"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
RESULTS_DIR = EXPERIMENT_ROOT / "results"
ANALYSIS_PATH = EXPERIMENT_ROOT / "logic_gate_analysis.json"

# Ideal logic gate patterns
IDEAL_GATES = {
    "AND": np.array([0, 0, 0, 1]),
    "OR": np.array([0, 1, 1, 1]),
    "NOT": np.array([1, 0, 1, 0]),  # B input inverts
    "XOR": np.array([0, 1, 1, 0]),
}

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def load_manifest() -> List[Dict]:
    """Load construct manifest."""
    with open(MANIFEST_PATH, 'r') as f:
        return json.load(f)


def load_prediction_stats(construct_name: str, cell_type: str) -> Dict:
    """Load prediction statistics for a construct."""
    stats_path = OUTPUT_DIR / f"{construct_name}_{cell_type}_stats.txt"
    
    if not stats_path.exists():
        return None
    
    stats = {}
    with open(stats_path, 'r') as f:
        for line in f:
            if ':' in line and not line.startswith('#'):
                key, value = line.strip().split(':', 1)
                try:
                    stats[key.strip()] = float(value.strip())
                except:
                    pass
    return stats


def normalize_signals(signals: np.ndarray) -> np.ndarray:
    """
    Normalize signals to [0, 1] range.
    Sets baseline (00) to 0, max to 1.
    """
    min_val = signals[0]  # 00 condition is baseline
    max_val = np.max(signals)
    
    if max_val - min_val < 1e-6:  # Avoid division by zero
        return np.zeros_like(signals)
    
    return (signals - min_val) / (max_val - min_val)


def compute_logic_score(observed: np.ndarray, gate_type: str) -> float:
    """
    Compute how well observed pattern matches ideal gate.
    Returns R² score (1.0 = perfect match, 0.0 = no correlation).
    """
    ideal = IDEAL_GATES[gate_type]
    
    # Normalize observed to [0, 1]
    observed_norm = normalize_signals(observed)
    
    # Compute R² score
    ss_res = np.sum((observed_norm - ideal) ** 2)
    ss_tot = np.sum((observed_norm - np.mean(observed_norm)) ** 2)
    
    if ss_tot < 1e-6:
        return 0.0
    
    r2 = 1 - (ss_res / ss_tot)
    return max(0.0, r2)  # Clip to [0, 1]


def compute_synergy_metrics(signals: Dict[str, float]) -> Dict:
    """
    Compute various synergy/cooperativity metrics.
    
    signals: {"00": baseline, "01": B_only, "10": A_only, "11": both}
    """
    s00, s01, s10, s11 = signals["00"], signals["01"], signals["10"], signals["11"]
    
    # Additivity ratio
    expected_additive = s01 + s10
    if expected_additive > 0:
        additivity = s11 / expected_additive
    else:
        additivity = np.nan
    
    # Excess over maximum
    max_single = max(s01, s10)
    eom = s11 - max_single
    
    # Bliss independence (multiplicative model)
    if s01 > 0 and s10 > 0:
        expected_bliss = s01 * s10 / s00 if s00 > 0 else s01 * s10
        bliss_excess = s11 - expected_bliss
    else:
        bliss_excess = np.nan
    
    # Synergy classification
    if additivity > 1.1:
        synergy_class = "synergy"
    elif additivity > 0.9:
        synergy_class = "additive"
    else:
        synergy_class = "interference"
    
    return {
        "additivity_ratio": additivity,
        "excess_over_max": eom,
        "bliss_excess": bliss_excess,
        "synergy_class": synergy_class,
    }


def analyze_tf_pair(pair_data: List[Dict]) -> Dict:
    """
    Analyze a single TF pair across all 4 truth table conditions.
    
    pair_data: List of 4 construct entries (00, 01, 10, 11)
    """
    # Sort by binary code to ensure correct order
    pair_data = sorted(pair_data, key=lambda x: x["binary_code"])
    
    # Extract signals (using max DNase as proxy for activity)
    signals = {}
    signal_array = np.zeros(4)
    
    for i, entry in enumerate(pair_data):
        binary_code = entry["binary_code"]
        stats = load_prediction_stats(entry["construct"], entry["cell_type"])
        
        if stats is None:
            return None  # Missing data
        
        signal = stats["max"]
        signals[binary_code] = signal
        signal_array[i] = signal
    
    # Compute logic scores for all gate types
    logic_scores = {}
    for gate_type in IDEAL_GATES.keys():
        score = compute_logic_score(signal_array, gate_type)
        logic_scores[gate_type] = score
    
    # Best-fit gate
    best_gate = max(logic_scores, key=logic_scores.get)
    best_score = logic_scores[best_gate]
    
    # Second-best for confidence
    scores_sorted = sorted(logic_scores.items(), key=lambda x: x[1], reverse=True)
    confidence = scores_sorted[0][1] - scores_sorted[1][1] if len(scores_sorted) > 1 else 1.0
    
    # Expected gate from experimental design
    expected_gate = pair_data[0]["gate_type"]
    correct_classification = (best_gate == expected_gate)
    
    # Synergy metrics
    synergy = compute_synergy_metrics(signals)
    
    # Normalized signals
    normalized = normalize_signals(signal_array)
    
    return {
        "tf_a": pair_data[0]["tf_a"],
        "tf_b": pair_data[0]["tf_b"],
        "cell_type": pair_data[0]["cell_type"],
        "expected_gate": expected_gate,
        "signals_raw": signals,
        "signals_array": signal_array.tolist(),
        "signals_normalized": normalized.tolist(),
        "logic_scores": logic_scores,
        "best_gate": best_gate,
        "best_score": best_score,
        "confidence": confidence,
        "correct_classification": correct_classification,
        "synergy_metrics": synergy,
        "biological_rationale": pair_data[0]["biological_rationale"],
    }


def analyze_all_pairs(manifest: List[Dict]) -> List[Dict]:
    """Analyze all TF pairs in manifest."""
    
    # Group by (tf_a, tf_b, cell_type)
    pairs = {}
    for entry in manifest:
        key = (entry["tf_a"], entry["tf_b"], entry["cell_type"])
        if key not in pairs:
            pairs[key] = []
        pairs[key].append(entry)
    
    print(f"\n{'='*80}")
    print(f"Analyzing {len(pairs)} TF pairs")
    print(f"{'='*80}\n")
    
    results = []
    for key, pair_data in pairs.items():
        tf_a, tf_b, cell_type = key
        
        if len(pair_data) != 4:
            print(f"⚠️  {tf_a}×{tf_b} ({cell_type}): Missing conditions (found {len(pair_data)}/4)")
            continue
        
        print(f"Analyzing: {tf_a} × {tf_b} ({cell_type})")
        
        result = analyze_tf_pair(pair_data)
        
        if result is None:
            print(f"  ✗ Missing prediction data")
            continue
        
        results.append(result)
        
        print(f"  Expected: {result['expected_gate']:8s} | Best fit: {result['best_gate']:8s} (score={result['best_score']:.3f})")
        print(f"  Signals: 00={result['signals_raw']['00']:.3f}, 01={result['signals_raw']['01']:.3f}, "
              f"10={result['signals_raw']['10']:.3f}, 11={result['signals_raw']['11']:.3f}")
        print(f"  Synergy: {result['synergy_metrics']['synergy_class']} (additivity={result['synergy_metrics']['additivity_ratio']:.2f})")
        
        if result['correct_classification']:
            print(f"  ✓ Correct classification!")
        else:
            print(f"  ✗ Misclassified (expected {result['expected_gate']})")
        print()
    
    return results


def generate_summary_statistics(results: List[Dict]) -> pd.DataFrame:
    """Generate summary statistics across all pairs."""
    
    # Overall accuracy
    total = len(results)
    correct = sum(1 for r in results if r["correct_classification"])
    accuracy = correct / total if total > 0 else 0
    
    print(f"\n{'='*80}")
    print(f"SUMMARY STATISTICS")
    print(f"{'='*80}\n")
    print(f"Overall classification accuracy: {correct}/{total} = {accuracy:.1%}")
    
    # Accuracy by expected gate type
    gate_stats = {}
    for gate_type in IDEAL_GATES.keys():
        gate_results = [r for r in results if r["expected_gate"] == gate_type]
        if gate_results:
            gate_correct = sum(1 for r in gate_results if r["correct_classification"])
            gate_stats[gate_type] = {
                "total": len(gate_results),
                "correct": gate_correct,
                "accuracy": gate_correct / len(gate_results),
                "mean_score": np.mean([r["best_score"] for r in gate_results]),
            }
    
    print(f"\nAccuracy by gate type:")
    for gate_type, stats in gate_stats.items():
        print(f"  {gate_type:8s}: {stats['correct']}/{stats['total']} = {stats['accuracy']:.1%} "
              f"(mean score: {stats['mean_score']:.3f})")
    
    # Synergy distribution
    synergy_counts = {}
    for r in results:
        sclass = r["synergy_metrics"]["synergy_class"]
        synergy_counts[sclass] = synergy_counts.get(sclass, 0) + 1
    
    print(f"\nSynergy classification:")
    for sclass, count in sorted(synergy_counts.items()):
        print(f"  {sclass:15s}: {count} pairs ({count/total:.1%})")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    return df


def save_results(results: List[Dict], df: pd.DataFrame):
    """Save analysis results."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save detailed JSON
    with open(ANALYSIS_PATH, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save CSV summary
    csv_path = RESULTS_DIR / "logic_gate_summary.csv"
    summary_df = df[[
        "tf_a", "tf_b", "cell_type", "expected_gate", "best_gate", "best_score",
        "confidence", "correct_classification"
    ]].copy()
    
    # Add synergy columns
    summary_df["additivity"] = df["synergy_metrics"].apply(lambda x: x["additivity_ratio"])
    summary_df["synergy_class"] = df["synergy_metrics"].apply(lambda x: x["synergy_class"])
    
    summary_df.to_csv(csv_path, index=False)
    
    print(f"\n✓ Results saved:")
    print(f"  Detailed JSON: {ANALYSIS_PATH}")
    print(f"  CSV summary: {csv_path}")


def main():
    print("=" * 80)
    print("LOGIC GATE ANALYSIS")
    print("=" * 80)
    
    # Load manifest
    manifest = load_manifest()
    print(f"\n✓ Loaded manifest: {len(manifest)} constructs")
    
    # Analyze all pairs
    results = analyze_all_pairs(manifest)
    
    if not results:
        print("\n✗ No results to analyze. Run predictions first!")
        return
    
    # Generate summary
    df = generate_summary_statistics(results)
    
    # Save results
    save_results(results, df)
    
    print(f"\n{'='*80}")
    print(f"🎉 Analysis complete!")
    print(f"{'='*80}\n")
    print(f"Next step: python create_logic_gate_figures.py")


if __name__ == "__main__":
    main()
