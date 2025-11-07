#!/usr/bin/env python3
"""
Parse XML FASTA sequences and build all enhancer stacking constructs.
Each construct is exactly 1,000,000 bp.
"""

import re
import os
from pathlib import Path

# Define paths
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
ENHANCER_FILE = ROOT / "sequences/enhancers/HS2_enhancer.fa"
PROMOTER_FILE = ROOT / "sequences/promoters/HBG1_promoter.fa"
FILLER_FILE = ROOT / "filler/1M_filler.txt"
OUTPUT_DIR = ROOT / "sequences/constructs"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_xml_fasta(filepath):
    """Extract clean DNA sequence from DAS XML format."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract DNA content between <DNA> tags
    match = re.search(r'<DNA[^>]*>(.*?)</DNA>', content, re.DOTALL)
    if not match:
        raise ValueError(f"No DNA sequence found in {filepath}")
    
    dna_content = match.group(1)
    # Remove all whitespace and newlines
    sequence = ''.join(dna_content.split()).upper()
    
    return sequence

def load_filler():
    """Load the 1M bp filler sequence and extend if needed."""
    with open(FILLER_FILE, 'r') as f:
        filler = f.read().strip().upper()
    
    # If we need more filler, repeat it cyclically
    # Need at least 1,048,576 bp for a single construct
    min_needed = 1_100_000  # A bit more than we need
    if len(filler) < min_needed:
        repeats = (min_needed // len(filler)) + 1
        filler = (filler * repeats)[:min_needed]
    
    return filler

def save_construct(name, sequence, output_dir):
    """Save a construct as a single-line FASTA file."""
    filepath = output_dir / f"{name}_construct.fa"
    with open(filepath, 'w') as f:
        f.write(f">{name}_construct\n")
        f.write(sequence + "\n")
    print(f"✓ Saved {name}: {len(sequence):,} bp → {filepath.name}")
    return filepath

def build_constructs():
    """Build all enhancer stacking constructs."""
    
    print("=" * 70)
    print("STEP 1: Parsing sequences from XML FASTA files")
    print("=" * 70)
    
    # Parse enhancer and promoter
    enhancer = parse_xml_fasta(ENHANCER_FILE)
    promoter = parse_xml_fasta(PROMOTER_FILE)
    filler = load_filler()
    
    print(f"✓ HS2 Enhancer: {len(enhancer)} bp")
    print(f"✓ HBG1 Promoter: {len(promoter)} bp")
    print(f"✓ Filler: {len(filler):,} bp")
    
    # Target construct length (must match AlphaGenome's supported length)
    # 1048576 = 2^20 = 1 MiB (not 1 MB)
    CONSTRUCT_LENGTH = 1_048_576
    
    # Position parameters
    PROMOTER_POS = 500_000  # Promoter centered at 500kb
    ENHANCER_POS = 400_000  # Enhancers at 400kb (100kb upstream)
    
    print("\n" + "=" * 70)
    print("STEP 2: Building constructs")
    print("=" * 70)
    
    constructs = {}
    
    # 1. FillerOnly: Just 1Mb of filler
    print("\n[1/10] Building FillerOnly...")
    constructs['FillerOnly'] = filler[:CONSTRUCT_LENGTH]
    save_construct('FillerOnly', constructs['FillerOnly'], OUTPUT_DIR)
    
    # 2. NoEnhancer: Promoter at center, rest is filler
    print("\n[2/10] Building NoEnhancer...")
    upstream_filler = filler[:PROMOTER_POS]
    downstream_length = CONSTRUCT_LENGTH - PROMOTER_POS - len(promoter)
    downstream_filler = filler[:downstream_length]
    constructs['NoEnhancer'] = upstream_filler + promoter + downstream_filler
    save_construct('NoEnhancer', constructs['NoEnhancer'], OUTPUT_DIR)
    
    # 3. E0: Enhancer immediately upstream of promoter (at PROMOTER_POS - enhancer length)
    print("\n[3/10] Building E0...")
    e0_enhancer_pos = PROMOTER_POS - len(enhancer)
    upstream_filler = filler[:e0_enhancer_pos]
    downstream_length = CONSTRUCT_LENGTH - e0_enhancer_pos - len(enhancer) - len(promoter)
    downstream_filler = filler[:downstream_length]
    constructs['E0'] = upstream_filler + enhancer + promoter + downstream_filler
    save_construct('E0', constructs['E0'], OUTPUT_DIR)
    
    # 4. E100: Single enhancer at 100kb upstream (ENHANCER_POS)
    print("\n[4/10] Building E100...")
    upstream_filler = filler[:ENHANCER_POS]
    middle_filler_length = PROMOTER_POS - ENHANCER_POS - len(enhancer)
    middle_filler = filler[:middle_filler_length]
    downstream_length = CONSTRUCT_LENGTH - PROMOTER_POS - len(promoter)
    downstream_filler = filler[:downstream_length]
    constructs['E100'] = upstream_filler + enhancer + middle_filler + promoter + downstream_filler
    save_construct('E100', constructs['E100'], OUTPUT_DIR)
    
    # 5-9. EC100-Nx: Multiple stacked enhancers at 100kb upstream
    copy_counts = [2, 5, 10, 160, 320]
    for i, n_copies in enumerate(copy_counts, start=5):
        name = f"EC100-{n_copies}x"
        print(f"\n[{i}/10] Building {name}...")
        
        # Stack n copies of the enhancer
        stacked_enhancers = enhancer * n_copies
        
        upstream_filler = filler[:ENHANCER_POS]
        middle_filler_length = PROMOTER_POS - ENHANCER_POS - len(stacked_enhancers)
        
        if middle_filler_length < 0:
            print(f"  ⚠️  Warning: {n_copies} enhancers ({len(stacked_enhancers):,} bp) overflow into promoter region!")
            print(f"     Adjusting to place enhancers just before promoter...")
            # Place enhancers immediately before promoter
            middle_filler_length = 0
            upstream_length = PROMOTER_POS - len(stacked_enhancers)
            upstream_filler = filler[:upstream_length]
        
        middle_filler = filler[:middle_filler_length]
        downstream_length = CONSTRUCT_LENGTH - len(upstream_filler) - len(stacked_enhancers) - len(middle_filler) - len(promoter)
        downstream_filler = filler[:downstream_length]
        
        constructs[name] = upstream_filler + stacked_enhancers + middle_filler + promoter + downstream_filler
        
        # Verify length
        actual_length = len(constructs[name])
        if actual_length != CONSTRUCT_LENGTH:
            print(f"  ⚠️  Length mismatch: {actual_length:,} bp (expected {CONSTRUCT_LENGTH:,} bp)")
            # Adjust by trimming/padding downstream filler
            if actual_length < CONSTRUCT_LENGTH:
                extra_filler = filler[:CONSTRUCT_LENGTH - actual_length]
                constructs[name] += extra_filler
            else:
                constructs[name] = constructs[name][:CONSTRUCT_LENGTH]
        
        save_construct(name, constructs[name], OUTPUT_DIR)
    
    print("\n" + "=" * 70)
    print("SUMMARY: All constructs built successfully!")
    print("=" * 70)
    
    # Print summary table
    print("\n{:<20} {:>15} {:>15}".format("Construct", "Length (bp)", "File"))
    print("-" * 70)
    for name, seq in constructs.items():
        print("{:<20} {:>15,} {:>15}".format(name, len(seq), f"{name}_construct.fa"))
    
    # Verify all sequences
    print("\n" + "=" * 70)
    print("VERIFICATION")
    print("=" * 70)
    all_valid = True
    for name, seq in constructs.items():
        if len(seq) != CONSTRUCT_LENGTH:
            print(f"❌ {name}: {len(seq):,} bp (expected {CONSTRUCT_LENGTH:,})")
            all_valid = False
        else:
            # Check for invalid characters
            invalid_chars = set(seq) - set('ATGC')
            if invalid_chars:
                print(f"❌ {name}: Contains invalid characters: {invalid_chars}")
                all_valid = False
    
    if all_valid:
        print(f"✅ All constructs are exactly {CONSTRUCT_LENGTH:,} bp with valid DNA sequences!")
    
    return constructs

if __name__ == "__main__":
    constructs = build_constructs()
    print(f"\n✅ Done! {len(constructs)} constructs saved to {OUTPUT_DIR}")
