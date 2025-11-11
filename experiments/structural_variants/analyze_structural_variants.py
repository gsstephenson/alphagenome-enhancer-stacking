#!/usr/bin/env python3
"""Analyze AlphaGenome structural-variant enhancer constructs."""

import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "structural_variants"
PREDICTION_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
RESULTS_DIR = EXPERIMENT_ROOT / "results"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 10,
})

PROMOTER_WINDOW = 5_000
ANCHOR_WINDOW = 1_000


def load_manifest() -> List[Dict]:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Manifest not found at {MANIFEST_PATH}")
    with open(MANIFEST_PATH, "r") as handle:
        manifest = json.load(handle)
    return manifest


def load_predictions(construct: str) -> np.ndarray:
    npy_path = PREDICTION_DIR / f"{construct}_dnase.npy"
    data = np.load(npy_path)
    if data.ndim > 1:
        data = data.squeeze()
    return data


def clamp_interval(start: int, end: int, length: int) -> slice:
    start = max(start, 0)
    end = min(end, length)
    return slice(start, end)


def compute_metrics(preds: np.ndarray, features: List[Dict]) -> Dict:
    length = len(preds)
    enhancer_feats = [f for f in features if f.get("label") == "hs2_block"]
    promoter_feat = next((f for f in features if f.get("label") == "promoter"), None)
    anchor_feats = sorted(
        (f for f in features if f.get("label") == "ctcf_anchor"),
        key=lambda item: item["start"],
    )

    metrics: Dict[str, float] = {}

    if enhancer_feats:
        enh = enhancer_feats[0]
        region = preds[enh["start"]:enh["end"]]
        metrics["enhancer_start"] = enh["start"]
        metrics["enhancer_end"] = enh["end"]
        metrics["enhancer_copies"] = enh.get("copies", np.nan)
        metrics["enhancer_max"] = float(region.max())
        metrics["enhancer_mean"] = float(region.mean())
        metrics["enhancer_auc"] = float(np.trapz(region))
    else:
        metrics.update({
            "enhancer_start": np.nan,
            "enhancer_end": np.nan,
            "enhancer_copies": 0,
            "enhancer_max": np.nan,
            "enhancer_mean": np.nan,
            "enhancer_auc": np.nan,
        })

    if promoter_feat:
        center = (promoter_feat["start"] + promoter_feat["end"]) // 2
        window = clamp_interval(center - PROMOTER_WINDOW, center + PROMOTER_WINDOW, length)
        metrics["promoter_start"] = promoter_feat["start"]
        metrics["promoter_end"] = promoter_feat["end"]
        metrics["promoter_mean"] = float(preds[window].mean())
        metrics["promoter_max"] = float(preds[window].max())
    else:
        metrics.update({
            "promoter_start": np.nan,
            "promoter_end": np.nan,
            "promoter_mean": np.nan,
            "promoter_max": np.nan,
        })

    if len(anchor_feats) >= 2:
        left, right = anchor_feats[0], anchor_feats[-1]
        metrics["loop_span_start"] = left["start"]
        metrics["loop_span_end"] = right["end"]
        span = preds[left["start"]:right["end"]]
        metrics["loop_span_mean"] = float(span.mean())
    else:
        metrics.update({
            "loop_span_start": np.nan,
            "loop_span_end": np.nan,
            "loop_span_mean": np.nan,
        })

    for idx, anchor in enumerate(anchor_feats):
        center = (anchor["start"] + anchor["end"]) // 2
        window = clamp_interval(center - ANCHOR_WINDOW, center + ANCHOR_WINDOW, length)
        metrics[f"anchor{idx+1}_mean"] = float(preds[window].mean())
        metrics[f"anchor{idx+1}_max"] = float(preds[window].max())
        metrics[f"anchor{idx+1}_orientation"] = anchor.get("orientation", "")

    metrics["global_mean"] = float(preds.mean())
    metrics["global_max"] = float(preds.max())
    metrics["global_std"] = float(preds.std())

    return metrics


def add_feature_spans(ax, features: List[Dict]) -> None:
    for feat in features:
        start = feat["start"] / 1_000
        end = feat["end"] / 1_000
        label = feat.get("label")
        if label == "hs2_block":
            ax.axvspan(start, end, color="red", alpha=0.15, label="Enhancer block")
        elif label == "promoter":
            ax.axvspan(start, end, color="blue", alpha=0.15, label="Promoter")
        elif label == "ctcf_anchor":
            ax.axvspan(start, end, color="orange", alpha=0.1, label=f"CTCF ({feat.get('orientation','')})")


def plot_tracks(predictions: Dict[str, np.ndarray], manifest: List[Dict]) -> Path:
    constructs = [entry["construct"] for entry in manifest]
    fig, axes = plt.subplots(len(constructs), 1, figsize=(14, 10), sharex=True)

    for ax, entry in zip(axes, manifest):
        name = entry["construct"]
        preds = predictions[name]
        x = np.arange(len(preds)) / 1_000  # kb
        ax.plot(x, preds, color="black", linewidth=0.7)
        add_feature_spans(ax, entry.get("features", []))
        ax.set_ylabel("DNase")
        ax.set_title(name, loc="left", fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Position (kb)")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels, fontsize=7, loc="upper right")
    plt.tight_layout()
    out_path = RESULTS_DIR / "structural_variant_tracks.png"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_metric_bars(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    sns.barplot(data=df, x="construct", y="enhancer_max", ax=axes[0], color="#b2182b")
    axes[0].set_title("Enhancer Max")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=30)

    sns.barplot(data=df, x="construct", y="promoter_mean", ax=axes[1], color="#2166ac")
    axes[1].set_title("Promoter Mean (±5 kb)")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=30)

    sns.barplot(data=df, x="construct", y="loop_span_mean", ax=axes[2], color="#4daf4a")
    axes[2].set_title("Loop Span Mean")
    axes[2].set_xlabel("")
    axes[2].tick_params(axis="x", rotation=30)

    for ax in axes:
        ax.set_ylabel("DNase")
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    out_path = RESULTS_DIR / "structural_variant_metrics.png"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    manifest = load_manifest()
    predictions = {entry["construct"]: load_predictions(entry["construct"]) for entry in manifest}

    records = []
    for entry in manifest:
        name = entry["construct"]
        metrics = compute_metrics(predictions[name], entry.get("features", []))
        metrics.update({
            "construct": name,
            "description": entry.get("description", ""),
        })
        records.append(metrics)

    df = pd.DataFrame(records)
    df_path = RESULTS_DIR / "structural_variant_metrics.csv"
    df.sort_values("construct").to_csv(df_path, index=False)

    track_path = plot_tracks(predictions, manifest)
    metric_path = plot_metric_bars(df)

    print(f"Metrics table saved to: {df_path}")
    print(f"Genome-wide tracks saved to: {track_path}")
    print(f"Metric comparison saved to: {metric_path}")


if __name__ == "__main__":
    main()
