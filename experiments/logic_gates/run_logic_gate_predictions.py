#!/usr/bin/env python3
"""
Run AlphaGenome predictions on all logic gate constructs using the Python API.
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np
from dotenv import load_dotenv

try:
    from alphagenome.models import dna_client
except ImportError:
    print("AlphaGenome package not available. Activate alphagenome-env before running.")
    import sys
    sys.exit(1)

# Load environment variables
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
load_dotenv(ROOT / ".env")

EXPERIMENT_ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking/experiments/logic_gates")
MANIFEST_PATH = EXPERIMENT_ROOT / "logic_gate_manifest.json"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
LOGS_DIR = EXPERIMENT_ROOT / "logs"
RESULTS_PATH = EXPERIMENT_ROOT / "prediction_results.json"

# Cell type mapping (EFO IDs)
CELL_TYPES = {
    "K562": "EFO:0002067",  # Erythroleukemia
    "HepG2": "EFO:0001187",  # Hepatocellular carcinoma
    "GM12878": "EFO:0002784",  # B-lymphocyte
}

# API setup
API_KEY = os.getenv("ALPHA_GENOME_KEY") or os.getenv("ALPHA_GENOME_API_KEY")
if not API_KEY:
    print("=" * 80)
    print("ERROR: AlphaGenome API key not found!")
    print("=" * 80)
    print(f"Looked for .env file at: {ROOT / '.env'}")
    print("Expected variable: ALPHA_GENOME_KEY or ALPHA_GENOME_API_KEY")
    print("=" * 80)
    raise ValueError("API key not configured")


def load_manifest() -> List[Dict]:
    """Load construct manifest."""
    with open(MANIFEST_PATH, 'r') as f:
        return json.load(f)


def load_fasta(path: Path) -> str:
    """Load sequence from FASTA file."""
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('>')]
    return "".join(lines).upper()


def run_prediction(
    sequence: str,
    construct_name: str,
    cell_type: str,
    client
) -> Dict:
    """Run AlphaGenome prediction for a single construct."""
    print(f"  Predicting {construct_name} in {cell_type} (length: {len(sequence)} bp)...")
    
    try:
        # Call AlphaGenome API
        cell_type_id = CELL_TYPES[cell_type]
        result = client.predict_sequence(
            sequence=sequence,
            requested_outputs=[dna_client.OutputType.DNASE],
            ontology_terms=[cell_type_id]
        )
        
        # Extract DNase predictions
        dnase_predictions = result.dnase.values
        
        # Get predictions for specific cell type
        cell_predictions = dnase_predictions[:, 0] if len(dnase_predictions.shape) > 1 else dnase_predictions
        
        # Also compute mean across all cell types for comparison
        mean_predictions = np.mean(dnase_predictions, axis=1) if len(dnase_predictions.shape) > 1 else dnase_predictions
        
        # Basic stats
        stats = {
            "max": float(np.max(cell_predictions)),
            "mean": float(np.mean(cell_predictions)),
            "std": float(np.std(cell_predictions)),
            "median": float(np.median(cell_predictions)),
            "global_max": float(np.max(mean_predictions)),
            "global_mean": float(np.mean(mean_predictions)),
        }
        
        return {
            "success": True,
            "predictions": cell_predictions,
            "mean_predictions": mean_predictions,
            "stats": stats,
            "cell_type": cell_type,
            "cell_type_id": cell_type_id,
        }
        
    except Exception as e:
        print(f"    ❌ Error: {e}")
        return {
            "success": False,
            "error": str(e),
            "cell_type": cell_type,
        }


def save_predictions(
    construct_name: str,
    cell_type: str,
    predictions: np.ndarray,
    mean_predictions: np.ndarray,
    stats: Dict
) -> None:
    """Save predictions to disk."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Filename with cell type
    base_name = f"{construct_name}_{cell_type}"
    
    # Save numpy array
    npy_path = OUTPUT_DIR / f"{base_name}_dnase.npy"
    np.save(npy_path, predictions)
    
    # Save mean across all cell types
    npy_mean_path = OUTPUT_DIR / f"{base_name}_dnase_mean.npy"
    np.save(npy_mean_path, mean_predictions)
    
    # Save text format (first 100 values for inspection)
    txt_path = OUTPUT_DIR / f"{base_name}_dnase.txt"
    with open(txt_path, 'w') as f:
        f.write(f"# DNase predictions for {construct_name} in {cell_type}\n")
        f.write(f"# Total bins: {len(predictions)}\n")
        f.write(f"# First 100 values:\n")
        for i, val in enumerate(predictions[:100]):
            f.write(f"{i}\t{val:.6f}\n")
    
    # Save stats
    stats_path = OUTPUT_DIR / f"{base_name}_stats.txt"
    with open(stats_path, 'w') as f:
        f.write(f"# AlphaGenome prediction statistics\n")
        f.write(f"# Construct: {construct_name}\n")
        f.write(f"# Cell type: {cell_type}\n")
        f.write(f"#\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")


def run_all_predictions(manifest: List[Dict], client):
    """Run predictions for all constructs in manifest."""
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    
    total = len(manifest)
    results = []
    
    print(f"\n{'='*80}")
    print(f"Running predictions for {total} constructs")
    print(f"{'='*80}\n")
    
    start_time = datetime.now()
    
    for i, entry in enumerate(manifest, 1):
        construct_name = entry["construct"]
        fasta_path = Path(entry["fasta_path"])
        cell_type = entry["cell_type"]
        
        print(f"[{i}/{total}] {construct_name}")
        print(f"  Gate: {entry['gate_type']}, Condition: {entry['binary_code']}")
        
        # Check if already exists
        output_npy = OUTPUT_DIR / f"{construct_name}_{cell_type}_dnase.npy"
        if output_npy.exists():
            print(f"  ⏭  Already exists, skipping")
            # Load existing stats for summary
            try:
                stats_path = OUTPUT_DIR / f"{construct_name}_{cell_type}_stats.txt"
                stats = {}
                with open(stats_path, 'r') as f:
                    for line in f:
                        if ':' in line and not line.startswith('#'):
                            key, val = line.strip().split(':', 1)
                            stats[key.strip()] = float(val.strip())
                results.append({
                    "construct": construct_name,
                    "cell_type": cell_type,
                    "success": True,
                    "skipped": True,
                    "stats": stats
                })
            except:
                results.append({
                    "construct": construct_name,
                    "cell_type": cell_type,
                    "success": True,
                    "skipped": True
                })
            continue
        
        # Load sequence
        try:
            sequence = load_fasta(fasta_path)
        except Exception as e:
            print(f"  ❌ Error loading FASTA: {e}")
            results.append({
                "construct": construct_name,
                "cell_type": cell_type,
                "success": False,
                "error": f"FASTA load error: {str(e)}"
            })
            continue
        
        # Run prediction
        pred_start = time.time()
        result = run_prediction(sequence, construct_name, cell_type, client)
        elapsed = time.time() - pred_start
        
        if result["success"]:
            # Save predictions
            save_predictions(
                construct_name,
                cell_type,
                result["predictions"],
                result["mean_predictions"],
                result["stats"]
            )
            
            print(f"  ✓ Success ({elapsed:.1f}s)")
            print(f"    Stats: max={result['stats']['max']:.3f}, mean={result['stats']['mean']:.6f}")
            
            results.append({
                "construct": construct_name,
                "cell_type": cell_type,
                "success": True,
                "stats": result["stats"],
                "elapsed_seconds": elapsed,
                "timestamp": datetime.now().isoformat()
            })
        else:
            print(f"  ✗ Failed: {result.get('error', 'Unknown error')}")
            results.append({
                "construct": construct_name,
                "cell_type": cell_type,
                "success": False,
                "error": result.get('error', 'Unknown error'),
                "elapsed_seconds": elapsed,
                "timestamp": datetime.now().isoformat()
            })
        
        # Brief pause to avoid overwhelming API
        time.sleep(0.5)
    
    # Save results summary
    with open(RESULTS_PATH, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Summary
    print(f"\n{'='*80}")
    success_count = sum(1 for r in results if r["success"])
    failed_count = total - success_count
    
    print(f"✓ Completed: {success_count}/{total} successful")
    if failed_count > 0:
        print(f"✗ Failed: {failed_count}/{total}")
    
    total_time = (datetime.now() - start_time).total_seconds()
    print(f"⏱  Total time: {total_time/60:.1f} minutes")
    print(f"✓ Results saved to: {RESULTS_PATH}")
    print(f"{'='*80}\n")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Run AlphaGenome predictions on logic gate constructs")
    parser.add_argument("--limit", type=int, help="Limit to first N constructs (for testing)")
    args = parser.parse_args()
    
    print("=" * 80)
    print("LOGIC GATE PREDICTION RUNNER")
    print("=" * 80)
    
    # Load manifest
    manifest = load_manifest()
    print(f"\n✓ Loaded manifest: {len(manifest)} constructs")
    
    if args.limit:
        manifest = manifest[:args.limit]
        print(f"  → Limited to first {args.limit} constructs")
    
    # Initialize AlphaGenome client
    print(f"\n✓ Initializing AlphaGenome client...")
    client = dna_client.create(API_KEY)
    
    # Run predictions
    run_all_predictions(manifest, client)
    
    print("\n🎉 Done! Next step:")
    print("   python analyze_logic_gates.py")


if __name__ == "__main__":
    main()
