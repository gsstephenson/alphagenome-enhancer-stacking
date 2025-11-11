#!/usr/bin/env python3
"""
Convert all sequences to clean FASTA format.

Handles:
- XML format (HS2 enhancer from DAS)
- Plain text (filler)
- Already-FASTA (promoter)
"""

from pathlib import Path
import re

# Paths
BASE_DIR = Path(__file__).parent
PARENT_DIR = BASE_DIR.parent.parent
ENHANCERS_DIR = PARENT_DIR / "sequences" / "enhancers"
PROMOTERS_DIR = PARENT_DIR / "sequences" / "promoters"
FILLER_DIR = PARENT_DIR / "filler"
OUTPUT_DIR = BASE_DIR / "clean_sequences"
OUTPUT_DIR.mkdir(exist_ok=True)


def extract_sequence_from_xml(filepath):
    """Extract DNA sequence from XML/DAS format."""
    with open(filepath) as f:
        content = f.read()
    
    # Find DNA content between <DNA> tags
    dna_match = re.search(r'<DNA[^>]*>(.*?)</DNA>', content, re.DOTALL)
    if dna_match:
        dna_content = dna_match.group(1)
        # Remove all whitespace and newlines
        sequence = ''.join(dna_content.split())
        return sequence.upper()
    
    return None


def load_fasta(filepath):
    """Load sequence from FASTA file."""
    with open(filepath) as f:
        lines = f.readlines()
    
    if lines[0].startswith('>'):
        seq = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    else:
        seq = ''.join(line.strip() for line in lines)
    
    return seq.upper()


def save_fasta(sequence, name, output_path, description=""):
    """Save sequence as FASTA with 80 bp per line."""
    with open(output_path, 'w') as f:
        if description:
            f.write(f">{name} {description}\n")
        else:
            f.write(f">{name}\n")
        
        # Write 80 bp per line
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')


def main():
    print("=" * 80)
    print("Converting Sequences to Clean FASTA Format")
    print("=" * 80)
    print()
    
    # 1. Convert HS2 enhancer (XML format)
    hs2_input = ENHANCERS_DIR / "HS2_enhancer.fa"
    hs2_output = OUTPUT_DIR / "HS2_enhancer.fa"
    
    print(f"Converting HS2 enhancer from XML...")
    hs2_seq = extract_sequence_from_xml(hs2_input)
    if hs2_seq:
        save_fasta(hs2_seq, "HS2_enhancer", hs2_output, "chr11:5290000-5291000")
        print(f"  ✓ {len(hs2_seq)} bp → {hs2_output}")
    else:
        print(f"  ✗ Failed to extract sequence from {hs2_input}")
    
    # 2. Convert HBG1 promoter (check if already FASTA)
    hbg1_input = PROMOTERS_DIR / "HBG1_promoter.fa"
    hbg1_output = OUTPUT_DIR / "HBG1_promoter.fa"
    
    print(f"Converting HBG1 promoter...")
    # Check if it's XML or FASTA
    with open(hbg1_input) as f:
        first_line = f.readline()
    
    if first_line.startswith('<?xml'):
        hbg1_seq = extract_sequence_from_xml(hbg1_input)
    else:
        hbg1_seq = load_fasta(hbg1_input)
    
    if hbg1_seq:
        save_fasta(hbg1_seq, "HBG1_promoter", hbg1_output, "chr11:5273600-5273900")
        print(f"  ✓ {len(hbg1_seq)} bp → {hbg1_output}")
    else:
        print(f"  ✗ Failed to load promoter sequence")
    
    # 3. Convert filler (plain text, no line breaks)
    filler_input = FILLER_DIR / "1M_filler.txt"
    filler_output = OUTPUT_DIR / "filler.fa"
    
    print(f"Converting filler sequence...")
    filler_seq = load_fasta(filler_input)
    save_fasta(filler_seq, "filler", filler_output, "A/T-rich background (1 MB)")
    print(f"  ✓ {len(filler_seq)} bp → {filler_output}")
    
    print()
    print("=" * 80)
    print("✓ All sequences converted to clean FASTA format")
    print("=" * 80)
    print()
    print("Verification:")
    print(f"  HS2 enhancer: {len(hs2_seq)} bp")
    print(f"  HBG1 promoter: {len(hbg1_seq)} bp")
    print(f"  Filler: {len(filler_seq)} bp")
    print()
    print("Next step: Run build_distance_constructs.py with updated paths")


if __name__ == "__main__":
    main()
