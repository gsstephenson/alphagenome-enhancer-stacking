#!/usr/bin/env python3
"""
Distance Decay Experiment - Construct Builder (with Technical Replicates)

Tests how AlphaGenome predictions change with enhancer-promoter distance.
Creates 8 distances × 3 technical replicates = 24 constructs.
Each replicate uses different random filler sequences but identical enhancer/promoter positions.

Author: Grant Stephenson
Date: November 11, 2025
Lab: Layer Laboratory, CU Boulder
"""

import json
import random
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent
SEQUENCES_DIR = BASE_DIR / "sequences"
SEQUENCES_DIR.mkdir(exist_ok=True)

# Paths
BASE_DIR = Path(__file__).parent
SEQUENCES_DIR = BASE_DIR / "sequences"
CLEAN_SEQ_DIR = BASE_DIR / "clean_sequences"
ENHANCER_FILE = CLEAN_SEQ_DIR / "HS2_enhancer.fa"
PROMOTER_FILE = CLEAN_SEQ_DIR / "HBG1_promoter.fa"
FILLER_FILE = CLEAN_SEQ_DIR / "filler.fa"

CONSTRUCT_LENGTH = 1048576  # 1 MiB (2^20)
PROMOTER_CENTER = CONSTRUCT_LENGTH // 2  # 524,288 bp

# Distance configurations (in bp, upstream of promoter)
DISTANCES = [1000, 5000, 10000, 25000, 50000, 100000, 200000, 500000]

# Technical replicates with different filler seeds
NUM_REPLICATES = 3
REPLICATE_SEEDS = [42, 123, 987]  # Random seeds for reproducibility


def load_sequence(filepath):
    """Load DNA sequence from FASTA or text file."""
    with open(filepath) as f:
        lines = f.readlines()
    
    # Handle FASTA format
    if lines[0].startswith('>'):
        seq = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    else:
        seq = ''.join(line.strip() for line in lines)
    
    return seq.upper()


def shuffle_filler(filler_seq, seed):
    """
    Shuffle filler sequence with given seed for technical replicates.
    Maintains nucleotide composition but randomizes order.
    """
    random.seed(seed)
    filler_list = list(filler_seq)
    random.shuffle(filler_list)
    return ''.join(filler_list)


def build_construct(enhancer_seq, promoter_seq, filler, distance_bp, replicate_seed=None):
    """
    Build a construct with enhancer at specified distance upstream of promoter.
    
    Layout:
        [filler][enhancer][spacer][promoter][filler]
    
    The promoter is centered at PROMOTER_CENTER.
    The enhancer is placed distance_bp upstream.
    
    Args:
        replicate_seed: If provided, shuffle filler with this seed for technical replicate
    """
    # Shuffle filler for technical replicates
    if replicate_seed is not None:
        filler = shuffle_filler(filler, replicate_seed)
    enhancer_len = len(enhancer_seq)
    promoter_len = len(promoter_seq)
    
    # Ensure filler is at least construct length
    if len(filler) < CONSTRUCT_LENGTH:
        # Repeat filler if needed
        repeats = (CONSTRUCT_LENGTH // len(filler)) + 1
        filler = filler * repeats
    
    # Promoter position (centered)
    promoter_start = PROMOTER_CENTER - promoter_len // 2
    promoter_end = promoter_start + promoter_len
    
    # Enhancer position (distance_bp upstream of promoter)
    enhancer_end = promoter_start - distance_bp
    enhancer_start = enhancer_end - enhancer_len
    
    if enhancer_start < 0:
        raise ValueError(f"Enhancer would start at negative position: {enhancer_start}")
    
    if promoter_end > CONSTRUCT_LENGTH:
        raise ValueError(f"Promoter would extend beyond construct: {promoter_end} > {CONSTRUCT_LENGTH}")
    
    # Build construct by filling in regions
    # Start with full filler, then replace enhancer and promoter regions
    construct_array = bytearray(filler[:CONSTRUCT_LENGTH].encode())
    
    # Insert enhancer
    construct_array[enhancer_start:enhancer_end] = enhancer_seq.encode()
    
    # Insert promoter
    construct_array[promoter_start:promoter_end] = promoter_seq.encode()
    
    final_seq = construct_array.decode()
    
    if len(final_seq) != CONSTRUCT_LENGTH:
        raise ValueError(f"Construct length {len(final_seq)} != {CONSTRUCT_LENGTH}")
    
    return final_seq, {
        'enhancer_start': enhancer_start,
        'enhancer_end': enhancer_end,
        'promoter_start': promoter_start,
        'promoter_end': promoter_end,
        'distance_bp': distance_bp,
        'enhancer_length': enhancer_len,
        'promoter_length': promoter_len
    }


def main():
    print("=" * 80)
    print("Distance Decay Experiment - Building Constructs")
    print("=" * 80)
    print()
    
    # Load components
    print("Loading components...")
    enhancer = load_sequence(ENHANCER_FILE)
    promoter = load_sequence(PROMOTER_FILE)
    filler = load_sequence(FILLER_FILE)
    
    print(f"  HS2 enhancer: {len(enhancer)} bp")
    print(f"  HBG1 promoter: {len(promoter)} bp")
    print(f"  Filler: {len(filler)} bp")
    print()
    
    # Build constructs with replicates
    manifest = {
        'experiment': 'distance_decay_replicates',
        'date': '2025-11-11',
        'construct_length': CONSTRUCT_LENGTH,
        'promoter_center': PROMOTER_CENTER,
        'num_replicates': NUM_REPLICATES,
        'replicate_seeds': REPLICATE_SEEDS,
        'enhancer': str(ENHANCER_FILE),
        'promoter': str(PROMOTER_FILE),
        'constructs': {}
    }
    
    total_constructs = len(DISTANCES) * NUM_REPLICATES
    print(f"Building {len(DISTANCES)} distances × {NUM_REPLICATES} replicates = {total_constructs} constructs...")
    print()
    
    for distance in DISTANCES:
        print(f"Distance {distance//1000}kb:")
        
        for rep_idx, seed in enumerate(REPLICATE_SEEDS, 1):
            name = f"Distance_{distance//1000}kb_rep{rep_idx}"
            print(f"  {name} (seed={seed})...", end=" ")
            
            try:
                seq, metadata = build_construct(enhancer, promoter, filler, distance, 
                                               replicate_seed=seed)
                
                # Save FASTA
                fasta_path = SEQUENCES_DIR / f"{name}.fa"
                with open(fasta_path, 'w') as f:
                    f.write(f">{name}\n")
                    # Write in 80 bp lines
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + '\n')
                
                # Store metadata
                metadata['replicate'] = rep_idx
                metadata['seed'] = seed
                manifest['constructs'][name] = metadata
                
                print(f"✓ Saved to {fasta_path.name}")
                
            except Exception as e:
                print(f"✗ ERROR: {e}")
        
        print()
    
    # Save manifest
    manifest_path = BASE_DIR / "construct_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"Manifest saved: {manifest_path}")
    print()
    print("=" * 80)
    print(f"✓ Successfully built {len(manifest['constructs'])} constructs")
    print(f"  {len(DISTANCES)} distances × {NUM_REPLICATES} replicates = {total_constructs} total")
    print("=" * 80)


if __name__ == "__main__":
    main()
