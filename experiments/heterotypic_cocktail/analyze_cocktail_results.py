#!/usr/bin/env python3
"""Analyze heterotypic enhancer cocktail predictions."""

import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "heterotypic_cocktail"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
RESULTS_DIR = EXPERIMENT_ROOT / "results"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 10,
})

MODULE_WINDOW = 1_000
PROMOTER_WINDOW = 5_000


def load_manifest() -> List[Dict]:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Manifest not found at {MANIFEST_PATH}")
    with open(MANIFEST_PATH, "r") as handle:
        return json.load(handle)


def load_predictions(construct: str) -> np.ndarray:
    npy_path = OUTPUT_DIR / f"{construct}_dnase.npy"
    arr = np.load(npy_path)
    if arr.ndim > 1:
        arr = arr.squeeze()
    return arr


def slice_window(arr: np.ndarray, start: int, end: int) -> np.ndarray:
    start = max(start, 0)
    end = min(end, len(arr))
    return arr[start:end]


def compute_metrics(preds: np.ndarray, features: List[Dict]) -> Dict:
    metrics: Dict[str, float] = {}
    module_means: Dict[str, List[float]] = {}

    promoter_feat = next((f for f in features if f.get("label") == "promoter"), None)

    for feat in features:
        if feat.get("label") == "enhancer_module":
            module_name = feat.get("module")
            module_means.setdefault(module_name, [])
            region = preds[feat["start"]:feat["end"]]
            module_means[module_name].append(float(region.mean()))

    for module_name, values in module_means.items():
        if values:
            metrics[f"module_mean_{module_name}"] = float(np.mean(values))
            metrics[f"module_max_{module_name}"] = float(np.max(values))

    if promoter_feat:
        center = (promoter_feat["start"] + promoter_feat["end"]) // 2
        window = slice_window(preds, center - PROMOTER_WINDOW, center + PROMOTER_WINDOW)
        metrics["promoter_mean"] = float(window.mean())
        metrics["promoter_max"] = float(window.max())
    else:
        metrics["promoter_mean"] = np.nan
        metrics["promoter_max"] = np.nan

    metrics["global_mean"] = float(preds.mean())
    metrics["global_max"] = float(preds.max())
    metrics["global_std"] = float(preds.std())

    return metrics


def plot_tracks(predictions: Dict[str, np.ndarray], manifest: List[Dict]) -> Path:
    constructs = [entry["construct"] for entry in manifest]
    fig, axes = plt.subplots(len(constructs), 1, figsize=(14, 12), sharex=True)

    for ax, entry in zip(axes, manifest):
        preds = predictions[entry["construct"]]
        x = np.arange(len(preds)) / 1_000
        ax.plot(x, preds, color="black", linewidth=0.6)
        for feat in entry["features"]:
            if feat.get("label") == "enhancer_module":
                ax.axvspan(feat["start"] / 1_000, feat["end"] / 1_000, alpha=0.1, color="red")
            elif feat.get("label") == "promoter":
                ax.axvspan(feat["start"] / 1_000, feat["end"] / 1_000, alpha=0.15, color="blue")
            elif feat.get("label") == "ctcf_bracket":
                ax.axvspan(feat["start"] / 1_000, feat["end"] / 1_000, alpha=0.1, color="orange")
        ax.set_ylabel("DNase")
        ax.set_title(entry["construct"], loc="left")
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Position (kb)")
    plt.tight_layout()
    out_path = RESULTS_DIR / "cocktail_tracks.png"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_module_heatmap(df: pd.DataFrame) -> Path:
    module_cols = [col for col in df.columns if col.startswith("module_mean_")]
    heatmap_df = df.set_index("construct")[module_cols]
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.heatmap(heatmap_df, cmap="viridis", annot=True, fmt=".4f", ax=ax)
    ax.set_title("Module mean accessibility")
    plt.tight_layout()
    out_path = RESULTS_DIR / "cocktail_module_heatmap.png"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    manifest = load_manifest()
    predictions = {entry["construct"]: load_predictions(entry["construct"]) for entry in manifest}

    records = []
    for entry in manifest:
        metrics = compute_metrics(predictions[entry["construct"]], entry.get("features", []))
        metrics.update({"construct": entry["construct"], "description": entry.get("description", "")})
        records.append(metrics)

    df = pd.DataFrame(records)
    metrics_path = RESULTS_DIR / "cocktail_metrics.csv"
    df.to_csv(metrics_path, index=False)

    track_path = plot_tracks(predictions, manifest)
    heatmap_path = plot_module_heatmap(df)

    print(f"Metrics saved to {metrics_path}")
    print(f"Tracks saved to {track_path}")
    print(f"Module heatmap saved to {heatmap_path}")


if __name__ == "__main__":
    main()
