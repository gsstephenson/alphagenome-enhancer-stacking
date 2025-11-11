#!/usr/bin/env python3
"""
Run AlphaGenome predictions for all regulatory grammar constructs.

Handles multiple cell types for cell-type specificity experiments.
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

# Load environment variables from .env file
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
load_dotenv(ROOT / ".env")

# Paths
EXPERIMENT_ROOT = ROOT / "experiments" / "regulatory_grammar"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"
SEQUENCE_DIR = EXPERIMENT_ROOT / "sequences"
OUTPUT_DIR = EXPERIMENT_ROOT / "alphagenome_outputs"
LOG_DIR = EXPERIMENT_ROOT / "logs"

# Cell type mapping (EFO IDs)
CELL_TYPES = {
    "K562": "EFO:0002067",  # Erythroleukemia
    "HepG2": "EFO:0001187",  # Hepatocellular carcinoma
    "GM12878": "EFO:0002784",  # B-lymphocyte
}

# API setup - loads from .env file
API_KEY = os.getenv("ALPHA_GENOME_KEY") or os.getenv("ALPHA_GENOME_API_KEY")
if not API_KEY:
    print("=" * 80)
    print("ERROR: AlphaGenome API key not found!")
    print("=" * 80)
    print(f"Looked for .env file at: {ROOT / '.env'}")
    print("Expected variable: ALPHA_GENOME_KEY or ALPHA_GENOME_API_KEY")
    print("\nMake sure the .env file exists and contains:")
    print("  ALPHA_GENOME_KEY=your_api_key_here")
    print("=" * 80)
    raise ValueError("API key not configured")


def load_fasta(path: Path) -> str:
    """Load sequence from FASTA file."""
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('>')]
    return "".join(lines).upper()


def get_cell_types_for_construct(construct: Dict) -> List[str]:
    """Determine which cell types to run for a construct."""
    experiment = construct["experiment"]
    
    if experiment == "cell_type_specificity":
        # Only run for the specific cell type indicated
        return [construct["cell_type"]]
    else:
        # All other experiments: just run K562 (erythroid context)
        return ["K562"]


def run_prediction(sequence: str, construct_name: str, cell_type: str, client) -> Dict:
    """Run AlphaGenome prediction for a single construct and cell type."""
    print(f"  Predicting {construct_name} in {cell_type}...")
    
    try:
        # Call AlphaGenome API
        cell_type_id = CELL_TYPES[cell_type]
        result = client.predict_sequence(
            sequence=sequence,
            requested_outputs=[dna_client.OutputType.DNASE],
            ontology_terms=[cell_type_id]
        )
        
        # Extract DNase predictions - use .values to get numpy array
        dnase_predictions = result.dnase.values
        
        # Get predictions for specific cell type (should be first since we requested it)
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


def save_predictions(construct_name: str, cell_type: str, predictions: np.ndarray, 
                     mean_predictions: np.ndarray, stats: Dict) -> None:
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
        f.write(f"Statistics for {construct_name} in {cell_type}\n")
        f.write(f"{'=' * 60}\n")
        for key, value in stats.items():
            f.write(f"{key:20s}: {value:.6f}\n")


def main():
    print("=" * 80)
    print("AlphaGenome Regulatory Grammar Predictions")
    print("=" * 80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Setup logging
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / f"predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Load manifest
    print(f"\nLoading manifest from {MANIFEST_PATH}...")
    with open(MANIFEST_PATH, 'r') as f:
        manifest = json.load(f)
    print(f"Found {len(manifest)} constructs")
    
    # Initialize client
    print("\nInitializing AlphaGenome client...")
    client = dna_client.create(API_KEY)
    
    # Track results
    results_summary = []
    total_predictions = sum(len(get_cell_types_for_construct(c)) for c in manifest)
    completed = 0
    failed = 0
    
    print(f"\nTotal predictions to run: {total_predictions}")
    print(f"Estimated time: ~{total_predictions * 1.5:.0f} minutes")
    print("=" * 80)
    
    # Process each construct
    for i, construct in enumerate(manifest, 1):
        construct_name = construct["construct"]
        experiment = construct["experiment"]
        fasta_path = EXPERIMENT_ROOT / construct["fasta"]
        
        print(f"\n[{i}/{len(manifest)}] {construct_name}")
        print(f"  Experiment: {experiment}")
        
        # Load sequence
        try:
            sequence = load_fasta(fasta_path)
            print(f"  Sequence length: {len(sequence):,} bp")
        except Exception as e:
            print(f"  ❌ Failed to load sequence: {e}")
            failed += 1
            continue
        
        # Get cell types to test
        cell_types_to_run = get_cell_types_for_construct(construct)
        print(f"  Cell types: {', '.join(cell_types_to_run)}")
        
        # Run predictions for each cell type
        for cell_type in cell_types_to_run:
            # Check if already exists
            output_file = OUTPUT_DIR / f"{construct_name}_{cell_type}_dnase.npy"
            if output_file.exists():
                print(f"  ⚠️  {cell_type}: Already exists, skipping...")
                completed += 1
                continue
            
            # Run prediction
            result = run_prediction(sequence, construct_name, cell_type, client)
            
            if result["success"]:
                # Save predictions
                save_predictions(
                    construct_name,
                    cell_type,
                    result["predictions"],
                    result["mean_predictions"],
                    result["stats"]
                )
                
                # Record summary
                results_summary.append({
                    "construct": construct_name,
                    "cell_type": cell_type,
                    "experiment": experiment,
                    "stats": result["stats"],
                    "timestamp": datetime.now().isoformat(),
                })
                
                completed += 1
                print(f"    ✓ Max DNase: {result['stats']['max']:.4f}")
                print(f"    ✓ Mean DNase: {result['stats']['mean']:.6f}")
            else:
                failed += 1
            
            # Rate limiting
            time.sleep(1)
    
    # Save results summary
    summary_path = EXPERIMENT_ROOT / "prediction_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    # Final report
    print("\n" + "=" * 80)
    print("PREDICTION SUMMARY")
    print("=" * 80)
    print(f"Total constructs: {len(manifest)}")
    print(f"Total predictions: {total_predictions}")
    print(f"Completed: {completed}")
    print(f"Failed: {failed}")
    print(f"Success rate: {100 * completed / total_predictions:.1f}%")
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Results saved to: {EXPERIMENT_ROOT}")
    print("=" * 80)


if __name__ == "__main__":
    main()
