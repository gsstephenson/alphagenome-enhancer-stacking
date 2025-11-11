#!/usr/bin/env python3
"""
Run AlphaGenome predictions on all logic gate constructs.

Reads manifest, submits predictions, and organizes outputs.
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


def run_alphagenome_prediction(
    fasta_path: Path,
    cell_type: str,
    output_prefix: str,
    track: str = "dnase"
) -> Dict:
    """
    Run AlphaGenome prediction on a single construct.
    
    NOTE: Adjust this function based on your AlphaGenome API!
    This is a template - you'll need to modify for your specific setup.
    """
    
    # Example command structure (MODIFY FOR YOUR SETUP):
    # python predict.py --fasta input.fa --cell-type K562 --track dnase --output out.npy
    
    output_npy = OUTPUT_DIR / f"{output_prefix}_{cell_type}_{track}.npy"
    output_txt = OUTPUT_DIR / f"{output_prefix}_{cell_type}_{track}.txt"
    output_stats = OUTPUT_DIR / f"{output_prefix}_{cell_type}_stats.txt"
    
    cmd = [
        "python", str(ALPHAGENOME_SCRIPT),
        "--fasta", str(fasta_path),
        "--cell-type", cell_type,
        "--track", track,
        "--output", str(output_npy),
    ]
    
    try:
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 min timeout
        )
        
        if result.returncode != 0:
            print(f"  ✗ AlphaGenome failed: {result.stderr}")
            return {"success": False, "error": result.stderr}
        
        # Load predictions and compute stats
        predictions = np.load(output_npy)
        
        stats = {
            "max": float(np.max(predictions)),
            "mean": float(np.mean(predictions)),
            "std": float(np.std(predictions)),
            "median": float(np.median(predictions)),
        }
        
        # Save stats
        with open(output_stats, 'w') as f:
            f.write(f"# AlphaGenome prediction statistics\n")
            f.write(f"# Construct: {output_prefix}\n")
            f.write(f"# Cell type: {cell_type}\n")
            f.write(f"# Track: {track}\n")
            f.write(f"#\n")
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        
        # Save human-readable predictions (first 100 values)
        with open(output_txt, 'w') as f:
            f.write(f"# First 100 prediction values\n")
            for i, val in enumerate(predictions.flatten()[:100]):
                f.write(f"{i}\t{val}\n")
        
        return {"success": True, "stats": stats, "output": str(output_npy)}
        
    except subprocess.TimeoutExpired:
        print(f"  ✗ Prediction timed out")
        return {"success": False, "error": "Timeout"}
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return {"success": False, "error": str(e)}


def run_all_predictions(manifest: List[Dict], dry_run: bool = False):
    """Run predictions for all constructs in manifest."""
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    
    total = len(manifest)
    results = []
    
    print(f"\n{'='*80}")
    print(f"Running predictions for {total} constructs")
    print(f"{'='*80}\n")
    
    for i, entry in enumerate(manifest, 1):
        construct_name = entry["construct"]
        fasta_path = Path(entry["fasta_path"])
        cell_type = entry["cell_type"]
        
        print(f"[{i}/{total}] {construct_name}")
        print(f"  Gate: {entry['gate_type']}, Condition: {entry['binary_code']}")
        
        if dry_run:
            print(f"  [DRY RUN] Would run: {fasta_path} → {cell_type}")
            results.append({"construct": construct_name, "success": True, "dry_run": True})
            continue
        
        # Check if already exists
        output_npy = OUTPUT_DIR / f"{construct_name}_{cell_type}_dnase.npy"
        if output_npy.exists():
            print(f"  ⏭  Already exists, skipping")
            results.append({"construct": construct_name, "success": True, "skipped": True})
            continue
        
        # Run prediction
        start_time = time.time()
        result = run_alphagenome_prediction(fasta_path, cell_type, construct_name)
        elapsed = time.time() - start_time
        
        result["construct"] = construct_name
        result["elapsed_seconds"] = elapsed
        results.append(result)
        
        if result["success"]:
            print(f"  ✓ Success ({elapsed:.1f}s)")
            print(f"    Stats: max={result['stats']['max']:.3f}, mean={result['stats']['mean']:.6f}")
        else:
            print(f"  ✗ Failed: {result.get('error', 'Unknown error')}")
        
        # Brief pause to avoid overwhelming API
        time.sleep(1)
    
    # Save results summary
    results_path = EXPERIMENT_ROOT / "prediction_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Summary
    print(f"\n{'='*80}")
    success_count = sum(1 for r in results if r["success"])
    print(f"✓ Completed: {success_count}/{total} successful")
    print(f"✓ Results saved to: {results_path}")
    print(f"{'='*80}\n")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Run AlphaGenome predictions on logic gate constructs")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be run without executing")
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
    
    # Check if AlphaGenome script exists
    if not args.dry_run and not ALPHAGENOME_SCRIPT.exists():
        print(f"\n⚠️  WARNING: AlphaGenome script not found at {ALPHAGENOME_SCRIPT}")
        print(f"   Please update ALPHAGENOME_SCRIPT path in this file!")
        print(f"   Running in dry-run mode...\n")
        args.dry_run = True
    
    # Run predictions
    run_all_predictions(manifest, dry_run=args.dry_run)
    
    print("\n🎉 Done! Next step:")
    print("   python analyze_logic_gates.py")


if __name__ == "__main__":
    main()
