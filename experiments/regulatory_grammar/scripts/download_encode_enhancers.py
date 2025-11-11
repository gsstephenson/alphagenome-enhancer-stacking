#!/usr/bin/env python3
"""
Download TF ChIP-seq enhancers from ENCODE

This script helps identify high-quality enhancer regions for each TF
by querying ENCODE cCREs (candidate cis-regulatory elements).

Usage:
    python download_encode_enhancers.py --tf GATA1 --cell-type K562 --output enhancers/
"""

import argparse
import requests
import json
import pandas as pd
from pathlib import Path


def query_encode_screen(tf_name, cell_type="K562", element_type="dELS"):
    """
    Query ENCODE SCREEN database for enhancers
    
    Args:
        tf_name: Transcription factor name (e.g., "GATA1")
        cell_type: Cell type (e.g., "K562")
        element_type: cCRE type
            - "dELS" = distal enhancer-like signature
            - "pELS" = proximal enhancer-like signature
            - "PLS" = promoter-like signature
    
    Returns:
        List of candidate enhancer coordinates
    """
    
    # ENCODE SCREEN API endpoint
    base_url = "https://api.wenglab.org/screen_v13/fftab"
    
    # Build query
    params = {
        "assembly": "GRCh38",
        "accession": "",
        "chromosome": "",
        "start": "",
        "end": "",
        "element_type": element_type,
        "rank_method": "top",
        "rank_start": 0,
        "rank_end": 100,
    }
    
    print(f"Querying ENCODE SCREEN for {tf_name} in {cell_type}...")
    print(f"Element type: {element_type}")
    print("\nNote: This is a template. You'll need to:")
    print("1. Go to https://screen.encodeproject.org/")
    print("2. Search for your TF and cell type")
    print("3. Download BED file of top enhancers")
    print("4. Or use the example template below\n")
    
    return []


def download_chipseq_peaks(tf_name, cell_type="K562"):
    """
    Download ChIP-seq peaks from ENCODE
    
    Returns URL to download BED file
    """
    
    # ENCODE portal API
    base_url = "https://www.encodeproject.org"
    search_url = f"{base_url}/search/"
    
    params = {
        "type": "Experiment",
        "assay_title": "TF ChIP-seq",
        "target.label": tf_name,
        "biosample_ontology.term_name": cell_type,
        "status": "released",
        "format": "json"
    }
    
    print(f"\nSearching ENCODE for {tf_name} ChIP-seq in {cell_type}...")
    
    try:
        response = requests.get(search_url, params=params)
        data = response.json()
        
        if '@graph' in data and len(data['@graph']) > 0:
            experiments = data['@graph']
            print(f"Found {len(experiments)} experiments")
            
            # Get first experiment
            exp = experiments[0]
            exp_id = exp['accession']
            
            print(f"\nExperiment: {exp_id}")
            print(f"URL: {base_url}{exp['@id']}")
            print(f"\nTo download peak files:")
            print(f"1. Visit: {base_url}/experiments/{exp_id}/")
            print(f"2. Click 'Files' tab")
            print(f"3. Download 'optimal IDR thresholded peaks' BED file")
            
            return exp_id
        else:
            print("No experiments found. Try manual search:")
            print(f"https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&target.label={tf_name}&biosample_ontology.term_name={cell_type}")
            
    except Exception as e:
        print(f"Error querying ENCODE: {e}")
        print("\nFalling back to manual search...")
    
    return None


def create_enhancer_template():
    """
    Create template BED file with example coordinates
    """
    
    template_data = {
        # TF enhancers from your existing experiments
        "HS2": {
            "chr": "chr11",
            "start": 5290000,
            "end": 5291001,
            "name": "HS2_beta_globin_LCR",
            "score": 1000,
            "strand": "+",
            "cell_type": "K562",
            "source": "Known β-globin LCR"
        },
        "GATA1": {
            "chr": "chrX",
            "start": 48649200,
            "end": 48650321,
            "name": "GATA1_enhancer",
            "score": 950,
            "strand": "+",
            "cell_type": "K562",
            "source": "ENCODE:ENCSR000EIL"
        },
        "KLF1": {
            "chr": "chr19",
            "start": 12990000,
            "end": 12990552,
            "name": "KLF1_enhancer",
            "score": 900,
            "strand": "+",
            "cell_type": "K562",
            "source": "ENCODE ChIP-seq"
        },
        "TAL1": {
            "chr": "chr1",
            "start": 47680000,
            "end": 47680974,
            "name": "TAL1_enhancer",
            "score": 880,
            "strand": "+",
            "cell_type": "K562",
            "source": "ENCODE ChIP-seq"
        },
        "HNF4A": {
            "chr": "chr20",
            "start": 44360000,
            "end": 44360502,
            "name": "HNF4A_enhancer",
            "score": 920,
            "strand": "+",
            "cell_type": "HepG2",
            "source": "ENCODE:ENCSR000EAF"
        },
    }
    
    return template_data


def generate_bed_template(output_file="enhancer_coordinates_template.bed"):
    """
    Generate BED file template with example enhancers
    """
    
    template = create_enhancer_template()
    
    # Convert to DataFrame
    rows = []
    for tf, coords in template.items():
        rows.append([
            coords['chr'],
            coords['start'],
            coords['end'],
            coords['name'],
            coords['score'],
            coords['strand']
        ])
    
    df = pd.DataFrame(rows, columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
    
    # Save as BED
    df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"\n✅ Created template BED file: {output_file}")
    print("\nBED format (tab-separated):")
    print(df.to_string(index=False))
    
    # Also save metadata
    metadata_file = output_file.replace('.bed', '_metadata.csv')
    metadata_rows = []
    for tf, coords in template.items():
        metadata_rows.append({
            'TF': tf,
            'name': coords['name'],
            'cell_type': coords['cell_type'],
            'source': coords['source'],
            'coordinates': f"{coords['chr']}:{coords['start']}-{coords['end']}"
        })
    
    metadata_df = pd.DataFrame(metadata_rows)
    metadata_df.to_csv(metadata_file, index=False)
    
    print(f"\n✅ Created metadata file: {metadata_file}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(description="Download TF enhancer coordinates")
    parser.add_argument("--tf", help="Transcription factor name")
    parser.add_argument("--cell-type", default="K562", help="Cell type")
    parser.add_argument("--output", default="enhancers/", help="Output directory")
    parser.add_argument("--template", action="store_true", help="Generate template BED file")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    if args.template:
        # Generate template
        output_file = output_dir / "enhancer_coordinates_template.bed"
        generate_bed_template(str(output_file))
        
        print("\n" + "="*60)
        print("INSTRUCTIONS FOR FINDING MORE ENHANCERS")
        print("="*60)
        print("\nMethod 1: ENCODE SCREEN (Recommended)")
        print("--------------------------------------")
        print("1. Go to: https://screen.encodeproject.org/")
        print("2. Search box: Enter TF name (e.g., 'MyoD')")
        print("3. Select cell type from dropdown")
        print("4. Filter: 'Distal enhancer-like signature' (dELS)")
        print("5. Sort by DNase signal (high to low)")
        print("6. Click on top enhancer → view coordinates")
        print("7. Download BED file or copy coordinates")
        
        print("\nMethod 2: ENCODE ChIP-seq")
        print("-------------------------")
        print("1. Go to: https://www.encodeproject.org/")
        print("2. Search: 'TF ChIP-seq [TF name] [cell type]'")
        print("3. Select experiment with highest quality")
        print("4. Click 'Files' tab")
        print("5. Download 'optimal IDR thresholded peaks' (BED narrowPeak)")
        print("6. Extract top peaks: sort -k7 -rn peaks.bed | head -10")
        
        print("\nMethod 3: Literature Mining")
        print("---------------------------")
        print("1. PubMed: '[TF] enhancer coordinates'")
        print("2. Check supplementary tables")
        print("3. Convert from hg19 to GRCh38 if needed")
        print("4. Use UCSC LiftOver: https://genome.ucsc.edu/cgi-bin/hgLiftOver")
        
        print("\nMethod 4: VISTA Enhancer Browser")
        print("---------------------------------")
        print("1. Go to: https://enhancer.lbl.gov/")
        print("2. Search by TF or tissue type")
        print("3. Use experimentally validated enhancers (gold standard)")
        
        print("\n" + "="*60)
        print("RECOMMENDED TF ENHANCERS TO FIND")
        print("="*60)
        print("\nPriority 1 (Muscle/Cardiac):")
        print("  - MyoD: Core regulatory region (chr11:17,700,000-17,800,000)")
        print("  - MEF2A: Muscle-specific enhancers")
        print("  - GATA4: Cardiac enhancers")
        print("  - NKX2-5: Heart development enhancers")
        
        print("\nPriority 2 (Immune/Myeloid):")
        print("  - PU.1 (SPI1): Myeloid enhancers")
        print("  - C/EBPα (CEBPA): Myeloid/adipose enhancers")
        print("  - NF-κB (RELA): Inflammatory enhancers")
        
        print("\nPriority 3 (Pluripotency):")
        print("  - Oct4 (POU5F1): Stem cell enhancers")
        print("  - Sox2: Neural/pluripotency enhancers")
        print("  - Nanog: Self-renewal enhancers")
        
        print("\n" + "="*60)
        
    else:
        # Query for specific TF
        if args.tf:
            download_chipseq_peaks(args.tf, args.cell_type)
            query_encode_screen(args.tf, args.cell_type)
        else:
            print("Error: Provide --tf or use --template")
            print("Example: python download_encode_enhancers.py --template")


if __name__ == "__main__":
    main()
