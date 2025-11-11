#!/usr/bin/env python3
"""Run AlphaGenome predictions for structural-variant constructs."""

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
    print("AlphaGenome package not available. Activate the correct environment before running.")
    sys.exit(1)

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "structural_variants"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
LOG_DIR = EXPERIMENT_ROOT / "logs"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

ONTOLOGY_TERM = "EFO:0002067"  # K562


def load_manifest() -> list:
    if MANIFEST_PATH.exists():
        with open(MANIFEST_PATH, "r") as handle:
            data = json.load(handle)
        return data
    constructs = []
    for fasta_path in sorted(CONSTRUCT_DIR.glob("*_construct.fa")):
        name = fasta_path.stem.replace("_construct", "")
        constructs.append({
            "construct": name,
            "description": "",
            "length": None,
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
        })
    return constructs


def load_fasta_sequence(path: Path) -> str:
    with open(path, "r") as handle:
        lines = [line.strip() for line in handle if not line.startswith(">")]
    return "".join(lines).upper()


def save_outputs(predictions: np.ndarray, construct_name: str) -> None:
    flat = predictions.reshape(-1)
    np.save(OUTPUT_DIR / f"{construct_name}_dnase.npy", predictions)
    np.savetxt(OUTPUT_DIR / f"{construct_name}_dnase.txt", flat, fmt="%.6f")
    stats_path = OUTPUT_DIR / f"{construct_name}_stats.txt"
    with open(stats_path, "w") as handle:
        handle.write(f"Construct: {construct_name}\n")
        handle.write(f"Timestamp: {datetime.now().isoformat()}\n")
        handle.write(f"Array shape: {predictions.shape}\n")
        handle.write(f"Min: {np.min(flat):.6f}\n")
        handle.write(f"Max: {np.max(flat):.6f}\n")
        handle.write(f"Mean: {np.mean(flat):.6f}\n")
        handle.write(f"Std: {np.std(flat):.6f}\n")


def main() -> None:
    load_dotenv()
    api_key = os.getenv("ALPHA_GENOME_API_KEY") or os.getenv("ALPHA_GENOME_KEY")
    if not api_key:
        print("AlphaGenome API key not found in environment. Set ALPHA_GENOME_KEY before running.")
        sys.exit(1)

    manifest = load_manifest()
    if not manifest:
        print("No constructs found. Run build_structural_variant_constructs.py first.")
        sys.exit(1)

    client = dna_client.create(api_key)
    log_file = LOG_DIR / f"predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    for entry in manifest:
        name = entry["construct"]
        fasta_rel = entry.get("fasta")
        fasta_path = EXPERIMENT_ROOT / fasta_rel if fasta_rel else CONSTRUCT_DIR / f"{name}_construct.fa"
        if not fasta_path.exists():
            print(f"Skipping {name}: FASTA not found ({fasta_path})")
            continue

        print(f"Predicting {name} ({fasta_path.name})")
        sequence = load_fasta_sequence(fasta_path)
        try:
            output = client.predict_sequence(
                sequence=sequence,
                requested_outputs=[dna_client.OutputType.DNASE],
                ontology_terms=[ONTOLOGY_TERM],
            )
        except Exception as exc:
            print(f"Prediction failed for {name}: {exc}")
            with open(log_file, "a") as handle:
                handle.write(f"{datetime.now().isoformat()} | FAILED | {name} | {exc}\n")
            continue

        dnase = output.dnase.values
        save_outputs(dnase, name)
        with open(log_file, "a") as handle:
            handle.write(f"{datetime.now().isoformat()} | SUCCESS | {name}\n")

    print("Prediction sweep finished.")


if __name__ == "__main__":
    main()
