#!/usr/bin/env python3
"""
Run AlphaGenome predictions for all enhancer stacking constructs.
Predicts DNase accessibility in K562 cell type.
"""

import os
import sys
import numpy as np
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# AlphaGenome imports
try:
    from alphagenome.models import dna_client
except ImportError:
    print("Error: AlphaGenome not found. Make sure you're in the alphagenome-env conda environment.")
    sys.exit(1)

# Define paths
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
CONSTRUCT_DIR = ROOT / "sequences/constructs"
OUTPUT_DIR = ROOT / "alphagenome/outputs"
LOG_DIR = ROOT / "logs"

# Ensure output directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Define all constructs
CONSTRUCTS = [
    'FillerOnly',
    'NoEnhancer',
    'E0',
    'E100',
    'EC100-2x',
    'EC100-5x',
    'EC100-10x',
    'EC100-160x',
    'EC100-320x'
]

CELL_TYPE = 'K562'

def load_fasta_sequence(fasta_path):
    """Load DNA sequence from FASTA file."""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    # Skip header line (starting with >)
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()

def run_alphagenome_prediction(dna_model, sequence, construct_name, ontology_term='EFO:0002067'):
    """
    Run AlphaGenome prediction on a DNA sequence.
    
    Args:
        dna_model: AlphaGenome DNA client model
        sequence: DNA sequence string
        construct_name: Name of the construct
        ontology_term: Cell type ontology (default: EFO:0002067 for K562)
    
    Returns:
        predictions: numpy array of DNase accessibility predictions
    """
    print(f"\n{'='*70}")
    print(f"Predicting {construct_name}")
    print(f"{'='*70}")
    print(f"Sequence length: {len(sequence):,} bp")
    print(f"Cell type ontology: {ontology_term} (K562)")
    
    try:
        # Run prediction
        print(f"Running AlphaGenome prediction...")
        output = dna_model.predict_sequence(
            sequence=sequence,
            requested_outputs=[dna_client.OutputType.DNASE],
            ontology_terms=[ontology_term],
        )
        
        # Extract DNase predictions
        dnase_predictions = output.dnase.values
        print(f"✓ Prediction complete: {dnase_predictions.shape}")
        
        return dnase_predictions
        
    except Exception as e:
        print(f"❌ Error during prediction: {e}")
        import traceback
        traceback.print_exc()
        return None

def save_predictions(predictions, construct_name, output_dir):
    """Save predictions in multiple formats."""
    
    # Flatten if multi-dimensional (predictions may have shape (length, 1) or (1, length))
    if predictions.ndim > 1:
        predictions_flat = predictions.flatten()
    else:
        predictions_flat = predictions
    
    # Save as numpy array (.npy)
    npy_path = output_dir / f"{construct_name}_dnase.npy"
    np.save(npy_path, predictions)
    print(f"✓ Saved numpy array: {npy_path.name}")
    
    # Save as text file for easy inspection
    txt_path = output_dir / f"{construct_name}_dnase.txt"
    np.savetxt(txt_path, predictions_flat, fmt='%.6f')
    print(f"✓ Saved text file: {txt_path.name}")
    
    # Save summary statistics
    stats_path = output_dir / f"{construct_name}_stats.txt"
    with open(stats_path, 'w') as f:
        f.write(f"Construct: {construct_name}\n")
        f.write(f"Cell type: K562 (EFO:0002067)\n")
        f.write(f"Prediction timestamp: {datetime.now().isoformat()}\n")
        f.write(f"\nPrediction array shape: {predictions.shape}\n")
        f.write(f"Flattened shape: {predictions_flat.shape}\n")
        f.write(f"Data type: {predictions.dtype}\n")
        f.write(f"\nStatistics:\n")
        f.write(f"  Min:    {np.min(predictions_flat):.6f}\n")
        f.write(f"  Max:    {np.max(predictions_flat):.6f}\n")
        f.write(f"  Mean:   {np.mean(predictions_flat):.6f}\n")
        f.write(f"  Median: {np.median(predictions_flat):.6f}\n")
        f.write(f"  Std:    {np.std(predictions_flat):.6f}\n")
    print(f"✓ Saved statistics: {stats_path.name}")
    
    return npy_path

def main():
    """Main prediction pipeline."""
    
    print("=" * 70)
    print("AlphaGenome Enhancer Stacking Predictions")
    print("=" * 70)
    print(f"Cell type: {CELL_TYPE}")
    print(f"Number of constructs: {len(CONSTRUCTS)}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Load API key from environment
    load_dotenv()
    api_key = os.getenv('ALPHA_GENOME_API_KEY') or os.getenv('ALPHA_GENOME_KEY')
    if not api_key:
        print("❌ Error: ALPHA_GENOME_API_KEY or ALPHA_GENOME_KEY not found in environment")
        print("Please set it in a .env file or export it as an environment variable")
        sys.exit(1)
    print(f"✓ API key loaded (length: {len(api_key)})")
    
    # Initialize AlphaGenome DNA client
    print("\n" + "=" * 70)
    print("Initializing AlphaGenome DNA client...")
    print("=" * 70)
    try:
        dna_model = dna_client.create(api_key)
        print("✓ DNA client initialized successfully")
    except Exception as e:
        print(f"❌ Error initializing DNA client: {e}")
        sys.exit(1)
    
    # K562 ontology term
    K562_ONTOLOGY = 'EFO:0002067'
    
    # Create log file
    log_file = LOG_DIR / f"predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    results = {}
    
    for i, construct_name in enumerate(CONSTRUCTS, 1):
        print(f"\n\n{'#'*70}")
        print(f"# [{i}/{len(CONSTRUCTS)}] Processing {construct_name}")
        print(f"{'#'*70}")
        
        # Load construct sequence
        fasta_path = CONSTRUCT_DIR / f"{construct_name}_construct.fa"
        
        if not fasta_path.exists():
            print(f"❌ Error: FASTA file not found: {fasta_path}")
            continue
        
        print(f"Loading sequence from {fasta_path.name}...")
        sequence = load_fasta_sequence(fasta_path)
        print(f"✓ Loaded {len(sequence):,} bp")
        
        # Run prediction
        predictions = run_alphagenome_prediction(dna_model, sequence, construct_name, K562_ONTOLOGY)
        
        if predictions is not None:
            # Save predictions
            output_path = save_predictions(predictions, construct_name, OUTPUT_DIR)
            results[construct_name] = {
                'sequence_length': len(sequence),
                'predictions': predictions,
                'output_path': output_path
            }
            
            # Log success
            with open(log_file, 'a') as f:
                f.write(f"{datetime.now().isoformat()} | SUCCESS | {construct_name}\n")
        else:
            # Log failure
            with open(log_file, 'a') as f:
                f.write(f"{datetime.now().isoformat()} | FAILED  | {construct_name}\n")
    
    # Final summary
    print("\n\n" + "=" * 70)
    print("PREDICTION SUMMARY")
    print("=" * 70)
    
    successful = len(results)
    failed = len(CONSTRUCTS) - successful
    
    print(f"Successful: {successful}/{len(CONSTRUCTS)}")
    print(f"Failed:     {failed}/{len(CONSTRUCTS)}")
    
    if results:
        print("\nResults saved to:")
        for name, info in results.items():
            print(f"  • {name}: {info['output_path'].name}")
    
    print(f"\nLog file: {log_file}")
    print("\n✅ Prediction pipeline complete!")
    
    return results

if __name__ == "__main__":
    results = main()
