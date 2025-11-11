#!/usr/bin/env python3
"""
Distance Decay Experiment - AlphaGenome Predictions

Runs AlphaGenome DNase predictions for distance series constructs.

Author: Grant Stephenson
Date: November 11, 2025
Lab: Layer Laboratory, CU Boulder
"""

import json
import os
import sys
import time
from pathlib import Path
import numpy as np
from dotenv import load_dotenv

# Load API key
env_path = Path(__file__).parent.parent.parent / ".env"
load_dotenv(env_path)

api_key = os.getenv("ALPHA_GENOME_KEY")
if not api_key:
    print("ERROR: ALPHA_GENOME_KEY not found in environment")
    print(f"Looked in: {env_path}")
    sys.exit(1)

# Import AlphaGenome
try:
    from alphagenome.models.dna_client import create
    from alphagenome.models.dna_output import OutputType
except ImportError:
    print("ERROR: alphagenome package not found")
    print("Install with: pip install alphagenome")
    sys.exit(1)

# Paths
BASE_DIR = Path(__file__).parent
SEQUENCES_DIR = BASE_DIR / "sequences"
OUTPUT_DIR = BASE_DIR / "alphagenome_outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

LOGS_DIR = BASE_DIR / "logs"
LOGS_DIR.mkdir(exist_ok=True)

# Cell type for K562
CELL_TYPE = "EFO:0002067"  # K562 ontology term


def load_fasta(fasta_path):
    """Load sequence from FASTA file."""
    with open(fasta_path) as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    name = header[1:] if header.startswith('>') else header
    seq = ''.join(line.strip() for line in lines[1:])
    
    return name, seq


def run_prediction(client, name, sequence):
    """Run AlphaGenome prediction on a sequence."""
    print()
    print("=" * 80)
    print(f"Predicting: {name}")
    print("=" * 80)
    print(f"Sequence length: {len(sequence):,} bp")
    
    try:
        print("Calling AlphaGenome API...")
        start_time = time.time()
        
        # Run prediction
        predictions = client.predict_sequence(
            sequence=sequence,
            requested_outputs=[OutputType.DNASE],
            ontology_terms=[CELL_TYPE]
        )
        
        elapsed = time.time() - start_time
        print(f"✓ Prediction completed in {elapsed:.1f} seconds")
        
        # Extract DNase predictions (from heterotypic_cocktail pattern)
        dnase_array = predictions.dnase.values  # Get the DNase values
        k562_dnase = dnase_array.mean(axis=1)  # Average across cell types
        
        print(f"Output shape: {dnase_array.shape}")
        print(f"K562 DNase stats:")
        print(f"  Min: {k562_dnase.min():.6f}")
        print(f"  Max: {k562_dnase.max():.6f}")
        print(f"  Mean: {k562_dnase.mean():.6f}")
        print(f"  Std: {k562_dnase.std():.6f}")
        
        # Save outputs
        npy_path = OUTPUT_DIR / f"{name}_dnase.npy"
        np.save(npy_path, k562_dnase)
        print(f"Saved: {npy_path}")
        
        txt_path = OUTPUT_DIR / f"{name}_dnase.txt"
        np.savetxt(txt_path, k562_dnase, fmt='%.6f')
        print(f"Saved: {txt_path}")
        
        stats_path = OUTPUT_DIR / f"{name}_stats.txt"
        with open(stats_path, 'w') as f:
            f.write(f"Construct: {name}\n")
            f.write(f"Sequence length: {len(sequence)} bp\n")
            f.write(f"Prediction time: {elapsed:.1f} seconds\n")
            f.write(f"\nDNase Statistics:\n")
            f.write(f"  Min: {k562_dnase.min():.6f}\n")
            f.write(f"  Max: {k562_dnase.max():.6f}\n")
            f.write(f"  Mean: {k562_dnase.mean():.6f}\n")
            f.write(f"  Std: {k562_dnase.std():.6f}\n")
        print(f"Saved: {stats_path}")
        
        return True
        
    except Exception as e:
        print(f"✗ Prediction failed: {e}")
        return False


def main():
    print("=" * 80)
    print("Distance Decay Experiment - AlphaGenome Predictions")
    print("=" * 80)
    print()
    
    # Load manifest
    manifest_path = BASE_DIR / "construct_manifest.json"
    if not manifest_path.exists():
        print(f"ERROR: Manifest not found: {manifest_path}")
        print("Run build_distance_constructs.py first")
        sys.exit(1)
    
    with open(manifest_path) as f:
        manifest = json.load(f)
    
    constructs = list(manifest['constructs'].keys())
    print(f"Found {len(constructs)} constructs to predict")
    print()
    
    # Check for existing predictions
    to_predict = []
    for name in constructs:
        npy_path = OUTPUT_DIR / f"{name}_dnase.npy"
        if npy_path.exists():
            print(f"  {name}: Already predicted (skipping)")
        else:
            print(f"  {name}: Needs prediction")
            to_predict.append(name)
    
    if not to_predict:
        print("\n✓ All predictions already complete!")
        return
    
    print(f"\nWill predict {len(to_predict)} constructs")
    print()
    
    # Create AlphaGenome client
    print("Creating AlphaGenome client...")
    client = create(api_key=api_key)
    print("✓ Client created")
    
    # Run predictions
    success_count = 0
    fail_count = 0
    
    for i, name in enumerate(to_predict, 1):
        print(f"\n[{i}/{len(to_predict)}] Processing {name}...")
        
        # Load sequence
        fasta_path = SEQUENCES_DIR / f"{name}.fa"
        _, sequence = load_fasta(fasta_path)
        
        # Run prediction
        success = run_prediction(client, name, sequence)
        
        if success:
            success_count += 1
        else:
            fail_count += 1
        
        # Brief pause between predictions
        if i < len(to_predict):
            print("\nWaiting 2 seconds before next prediction...")
            time.sleep(2)
    
    # Summary
    print("\n" + "=" * 80)
    print("PREDICTION SUMMARY")
    print("=" * 80)
    print(f"  Successful: {success_count}/{len(to_predict)}")
    print(f"  Failed: {fail_count}/{len(to_predict)}")
    print("=" * 80)
    
    if fail_count == 0:
        print("\n✓ All predictions completed successfully!")
    else:
        print(f"\n⚠ {fail_count} predictions failed")


if __name__ == "__main__":
    main()
