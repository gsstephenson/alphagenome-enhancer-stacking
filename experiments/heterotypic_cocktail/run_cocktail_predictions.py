#!/usr/bin/env python3
"""Run AlphaGenome predictions for heterotypic enhancer cocktails."""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from dotenv import load_dotenv

try:
    from alphagenome.models import dna_client
except ImportError:
    print("AlphaGenome package not available. Activate alphagenome-env before running.")
    sys.exit(1)

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "heterotypic_cocktail"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
LOG_DIR = EXPERIMENT_ROOT / "logs"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

ONTOLOGY_TERM = "EFO:0002067"  # K562


def load_manifest() -> list:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Manifest not found at {MANIFEST_PATH}")
    with open(MANIFEST_PATH, "r") as handle:
        return json.load(handle)


def load_sequence(fasta_path: Path) -> str:
    with open(fasta_path, "r") as handle:
        lines = [line.strip() for line in handle if not line.startswith(">")]
    return "".join(lines).upper()


def save_outputs(predictions: np.ndarray, construct: str) -> None:
    flat = predictions.reshape(-1)
    np.save(OUTPUT_DIR / f"{construct}_dnase.npy", predictions)
    np.savetxt(OUTPUT_DIR / f"{construct}_dnase.txt", flat, fmt="%.6f")
    stats_path = OUTPUT_DIR / f"{construct}_stats.txt"
    with open(stats_path, "w") as handle:
        handle.write(f"Construct: {construct}\n")
        handle.write(f"Timestamp: {datetime.now().isoformat()}\n")
        handle.write(f"Shape: {predictions.shape}\n")
        handle.write(f"Min: {flat.min():.6f}\n")
        handle.write(f"Max: {flat.max():.6f}\n")
        handle.write(f"Mean: {flat.mean():.6f}\n")
        handle.write(f"Std: {flat.std():.6f}\n")


def main() -> None:
    load_dotenv()
    api_key = os.getenv("ALPHA_GENOME_API_KEY") or os.getenv("ALPHA_GENOME_KEY")
    if not api_key:
        print("AlphaGenome API key not configured.")
        sys.exit(1)

    manifest = load_manifest()
    client = dna_client.create(api_key)
    log_path = LOG_DIR / f"predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    for entry in manifest:
        name = entry["construct"]
        fasta_rel = entry["fasta"]
        fasta_path = EXPERIMENT_ROOT / fasta_rel
        print(f"Predicting {name}...")
        sequence = load_sequence(fasta_path)
        try:
            output = client.predict_sequence(
                sequence=sequence,
                requested_outputs=[dna_client.OutputType.DNASE],
                ontology_terms=[ONTOLOGY_TERM],
            )
        except Exception as exc:
            print(f"Prediction failed for {name}: {exc}")
            with open(log_path, "a") as handle:
                handle.write(f"{datetime.now().isoformat()} | FAILED | {name} | {exc}\n")
            continue

        save_outputs(output.dnase.values, name)
        with open(log_path, "a") as handle:
            handle.write(f"{datetime.now().isoformat()} | SUCCESS | {name}\n")

    print("Prediction sweep finished.")


if __name__ == "__main__":
    main()
