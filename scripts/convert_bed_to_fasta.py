#!/usr/bin/env python3
"""
Convert BED files to FASTA format by fetching sequences from UCSC Genome Browser.

This script reads BED files containing genomic coordinates and retrieves the
corresponding DNA sequences using the UCSC REST API, then writes them to FASTA format.

Usage:
    python convert_bed_to_fasta.py

The script will process all .BED files in the enhancers directory and create
corresponding .fa files with the sequences.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List, Tuple, Optional
import urllib.request
import urllib.error
import json

# Paths
ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
ENHANCER_DIR = ROOT / "sequences" / "enhancers"
OUTPUT_DIR = ROOT / "sequences" / "enhancers"

# UCSC API endpoint
UCSC_API = "https://api.genome.ucsc.edu/getData/sequence"

# RefSeq to UCSC chromosome mapping
REFSEQ_TO_UCSC = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}


def parse_bed_line(line: str) -> Optional[Tuple[str, int, int, str]]:
    """Parse a single BED line and return (chrom, start, end, name)."""
    parts = line.strip().split('\t')
    if len(parts) < 3:
        return None
    
    chrom = parts[0]
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    
    # Get feature name if available, otherwise use coordinates
    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
    
    return (chrom, start, end, name)


def refseq_to_ucsc_chrom(refseq: str) -> Optional[str]:
    """Convert RefSeq accession to UCSC chromosome name."""
    return REFSEQ_TO_UCSC.get(refseq)


def fetch_sequence_from_ucsc(chrom: str, start: int, end: int, genome: str = "hg38") -> Optional[str]:
    """
    Fetch DNA sequence from UCSC Genome Browser API.
    
    Args:
        chrom: Chromosome name (e.g., 'chr16')
        start: Start position (0-based)
        end: End position (exclusive)
        genome: Genome assembly (default: hg38)
    
    Returns:
        DNA sequence as string, or None if request fails
    """
    url = f"{UCSC_API}?genome={genome};chrom={chrom};start={start};end={end}"
    
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.loads(response.read().decode())
            if 'dna' in data:
                return data['dna'].upper()
            else:
                print(f"  ⚠️  No sequence data in response for {chrom}:{start}-{end}")
                return None
    except urllib.error.URLError as e:
        print(f"  ❌ Error fetching {chrom}:{start}-{end}: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"  ❌ Error parsing JSON for {chrom}:{start}-{end}: {e}")
        return None


def process_bed_file(bed_path: Path, merge_strategy: str = "first", feature_filter: str = "enhancer") -> List[Tuple[str, str]]:
    """
    Process a BED file and fetch sequences for all regions.
    
    Args:
        bed_path: Path to BED file
        merge_strategy: How to handle multiple regions:
            - "first": Use only the first region
            - "longest": Use the longest region
            - "all": Concatenate all regions
            - "separate": Create separate FASTA entries
        feature_filter: Only include regions with this feature type (default: "enhancer")
                       Set to None to include all features
    
    Returns:
        List of (name, sequence) tuples
    """
    print(f"\nProcessing {bed_path.name}...")
    
    regions = []
    filtered_regions = []
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parsed = parse_bed_line(line)
            if parsed:
                regions.append(parsed)
                # Filter by feature type
                chrom, start, end, name = parsed
                if feature_filter is None or feature_filter.lower() in name.lower():
                    filtered_regions.append(parsed)
    
    if not regions:
        print(f"  ⚠️  No valid regions found in {bed_path.name}")
        return []
    
    print(f"  Found {len(regions)} region(s) total")
    
    # Report filtering results
    if feature_filter:
        if len(filtered_regions) > 0 and len(filtered_regions) < len(regions):
            print(f"  Filtered to {len(filtered_regions)} region(s) matching '{feature_filter}'")
            excluded = len(regions) - len(filtered_regions)
            print(f"  Excluded {excluded} non-enhancer region(s) (silencers, etc.)")
            regions = filtered_regions
        elif len(filtered_regions) == 0:
            # No exact matches for "enhancer", try broader match
            print(f"  ⚠️  No regions with '{feature_filter}' label found")
            print(f"  Using 'transcriptional_cis_regulatory_region' as fallback (excluding silencers)")
            for parsed in regions:
                chrom, start, end, name = parsed
                # Accept transcriptional regulatory regions but NOT silencers
                if "transcriptional_cis_regulatory_region" in name.lower() or ("regulatory" in name.lower() and "silencer" not in name.lower()):
                    filtered_regions.append(parsed)
            if filtered_regions:
                print(f"  Found {len(filtered_regions)} regulatory region(s)")
                regions = filtered_regions
            else:
                print(f"  ⚠️  No regulatory regions found")
                return []
        else:
            regions = filtered_regions
    
    if not regions:
        print(f"  ⚠️  No valid regions after filtering")
        return []
    
    # Convert RefSeq to UCSC if needed
    processed_regions = []
    for chrom, start, end, name in regions:
        if chrom.startswith("NC_"):
            ucsc_chrom = refseq_to_ucsc_chrom(chrom)
            if ucsc_chrom:
                print(f"  Converting {chrom} → {ucsc_chrom}")
                chrom = ucsc_chrom
            else:
                print(f"  ⚠️  Unknown RefSeq accession: {chrom}")
                continue
        processed_regions.append((chrom, start, end, name))
    
    if not processed_regions:
        print(f"  ⚠️  No valid regions after conversion")
        return []
    
    # Apply merge strategy
    if merge_strategy == "first":
        processed_regions = [processed_regions[0]]
        print(f"  Using first region only: {processed_regions[0][3]}")
    elif merge_strategy == "longest":
        longest = max(processed_regions, key=lambda r: r[2] - r[1])
        processed_regions = [longest]
        print(f"  Using longest region ({longest[2] - longest[1]} bp): {longest[3]}")
    
    # Fetch sequences
    sequences = []
    for chrom, start, end, name in processed_regions:
        print(f"  Fetching {chrom}:{start}-{end} ({end - start} bp)...")
        seq = fetch_sequence_from_ucsc(chrom, start, end)
        
        if seq:
            sequences.append((name, seq))
            print(f"    ✓ Retrieved {len(seq)} bp")
        else:
            print(f"    ✗ Failed to retrieve sequence")
        
        # Rate limiting - be nice to UCSC servers
        time.sleep(0.5)
    
    return sequences


def write_fasta(sequences: List[Tuple[str, str]], output_path: Path, base_name: str) -> None:
    """
    Write sequences to FASTA format.
    
    Args:
        sequences: List of (name, sequence) tuples
        output_path: Output file path
        base_name: Base name for the enhancer (used in header)
    """
    if not sequences:
        print(f"  ⚠️  No sequences to write")
        return
    
    with open(output_path, 'w') as f:
        if len(sequences) == 1:
            # Single sequence
            name, seq = sequences[0]
            f.write(f">{base_name}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
        else:
            # Multiple sequences
            for idx, (name, seq) in enumerate(sequences, 1):
                f.write(f">{base_name}_region{idx} {name}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
    
    print(f"  ✓ Wrote {output_path.name} ({sum(len(s[1]) for s in sequences)} bp total)")


def main():
    """Main entry point."""
    print("=" * 70)
    print("BED to FASTA Converter")
    print("=" * 70)
    
    # Find all BED files
    bed_files = list(ENHANCER_DIR.glob("*.BED")) + list(ENHANCER_DIR.glob("*.bed"))
    
    if not bed_files:
        print("\n❌ No BED files found in", ENHANCER_DIR)
        return
    
    print(f"\nFound {len(bed_files)} BED file(s):")
    for bf in bed_files:
        print(f"  - {bf.name}")
    
    # Process each BED file
    for bed_path in bed_files:
        # Determine output name
        base_name = bed_path.stem.replace("_enhancer", "")
        output_path = OUTPUT_DIR / f"{base_name}_module.fa"
        
        # Check if output already exists
        if output_path.exists():
            print(f"\n⚠️  {output_path.name} already exists. Skipping...")
            continue
        
        # Process BED file
        # For enhancers with multiple regions, use the longest one
        # Filter to only include enhancers (not silencers)
        sequences = process_bed_file(bed_path, merge_strategy="longest", feature_filter="enhancer")
        
        # Write FASTA
        if sequences:
            write_fasta(sequences, output_path, base_name)
        else:
            print(f"  ❌ No sequences retrieved for {bed_path.name}")
    
    print("\n" + "=" * 70)
    print("Conversion complete!")
    print("=" * 70)
    
    # Summary
    fa_files = list(OUTPUT_DIR.glob("*_module.fa"))
    print(f"\nTotal FASTA files in {OUTPUT_DIR}:")
    for fa in sorted(fa_files):
        print(f"  ✓ {fa.name}")


if __name__ == "__main__":
    main()
