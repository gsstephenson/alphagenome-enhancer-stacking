#!/usr/bin/env python3
"""
Build comprehensive regulatory grammar test constructs for AlphaGenome.

Tests:
1. Cell-type specificity (3 promoters × 3 enhancers × 3 cell types)
2. Pairwise motif cooperativity (10 TF pairs + 5 singles + 2 CTCF controls)
3. Short-range spacing (10 distances for HS2+GATA1)
4. Orientation effects (4 pairs × 4 orientations)

Total: 78 constructs
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "regulatory_grammar"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

# Input sequences
ENHANCER_DIR = ROOT / "sequences" / "enhancers"
PROMOTER_DIR = ROOT / "sequences" / "promoters"
FILLER_PATH = ROOT / "filler" / "1M_filler.txt"

# Construct parameters
CONSTRUCT_LENGTH = 1_048_576
PROMOTER_POS = 500_000
ENHANCER_DOMAIN_START = 250_000


@dataclass
class Module:
    name: str
    sequence: str

    def oriented(self, orientation: str) -> str:
        if orientation == "+":
            return self.sequence
        if orientation == "-":
            return reverse_complement(self.sequence)
        raise ValueError(f"Unknown orientation '{orientation}'")


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(table)[::-1]


def load_fasta_sequence(path: Path) -> str:
    """Load sequence from FASTA file (handles both plain and XML-style)."""
    text = path.read_text()
    
    # Try XML format first (like HS2)
    match = re.search(r"<DNA[^>]*>(.*?)</DNA>", text, re.DOTALL)
    if match:
        return "".join(match.group(1).split()).upper()
    
    # Plain FASTA format
    lines = [line.strip() for line in text.split('\n') if not line.startswith('>')]
    return "".join(lines).upper()


def load_filler(path: Path) -> str:
    return path.read_text().strip().upper()


class SequenceBuilder:
    def __init__(self, filler: str):
        self._filler = filler
        self._filler_idx = 0
        self.parts: List[str] = []
        self.cursor = 0
        self.features: List[Dict] = []

    def _take_filler(self, length: int) -> str:
        if length <= 0:
            return ""
        chunks = []
        remaining = length
        while remaining > 0:
            chunk_len = min(len(self._filler) - self._filler_idx, remaining)
            chunks.append(self._filler[self._filler_idx:self._filler_idx + chunk_len])
            self._filler_idx = (self._filler_idx + chunk_len) % len(self._filler)
            remaining -= chunk_len
        return "".join(chunks)

    def append_filler(self, length: int, label: Optional[str] = None) -> None:
        if length <= 0:
            return
        seq = self._take_filler(length)
        self._append(seq, label)

    def append_module(self, module: Module, orientation: str, label: str) -> None:
        seq = module.oriented(orientation)
        meta = {"module": module.name, "orientation": orientation}
        self._append(seq, label, meta)

    def _append(self, seq: str, label: Optional[str], metadata: Optional[Dict] = None) -> None:
        if not seq:
            return
        start = self.cursor
        self.parts.append(seq)
        self.cursor += len(seq)
        if label:
            feature = {"label": label, "start": start, "end": self.cursor}
            if metadata:
                feature.update(metadata)
            self.features.append(feature)

    def finish(self, total_length: int) -> str:
        if self.cursor > total_length:
            raise ValueError(f"Construct exceeds target ({self.cursor} > {total_length})")
        self.append_filler(total_length - self.cursor)
        return "".join(self.parts)


def build_celltype_construct(
    promoter: Module,
    enhancer: Module,
    filler: str,
    spacing: int = 100_000
) -> Dict:
    """Build promoter + single enhancer construct."""
    builder = SequenceBuilder(filler)
    builder.append_filler(ENHANCER_DOMAIN_START, "upstream_filler")
    builder.append_module(enhancer, "+", "enhancer")
    builder.append_filler(spacing - len(enhancer.sequence), "spacer")
    builder.append_module(promoter, "+", "promoter")
    sequence = builder.finish(CONSTRUCT_LENGTH)
    return {"sequence": sequence, "features": builder.features}


def build_pairwise_construct(
    enhancer1: Module,
    enhancer2: Optional[Module],
    promoter: Module,
    filler: str,
    spacing: int = 5_000,
    separator: Optional[Module] = None
) -> Dict:
    """Build construct with two enhancers (or one if enhancer2 is None)."""
    builder = SequenceBuilder(filler)
    builder.append_filler(ENHANCER_DOMAIN_START, "upstream_filler")
    
    builder.append_module(enhancer1, "+", "enhancer1")
    
    if enhancer2:
        if separator:
            builder.append_module(separator, "+", "separator")
        builder.append_filler(spacing, "inter_enhancer_spacing")
        builder.append_module(enhancer2, "+", "enhancer2")
        
        # Calculate remaining space to promoter
        used = builder.cursor - ENHANCER_DOMAIN_START
        remaining = (PROMOTER_POS - ENHANCER_DOMAIN_START) - used
        builder.append_filler(remaining, "spacer_to_promoter")
    else:
        # Single enhancer
        builder.append_filler(PROMOTER_POS - builder.cursor, "spacer_to_promoter")
    
    builder.append_module(promoter, "+", "promoter")
    sequence = builder.finish(CONSTRUCT_LENGTH)
    return {"sequence": sequence, "features": builder.features}


def build_spacing_construct(
    enhancer1: Module,
    enhancer2: Module,
    promoter: Module,
    filler: str,
    spacing: int
) -> Dict:
    """Build construct with specific spacing between enhancers."""
    return build_pairwise_construct(enhancer1, enhancer2, promoter, filler, spacing, None)


def build_orientation_construct(
    enhancer1: Module,
    enhancer2: Module,
    promoter: Module,
    filler: str,
    orient1: str,
    orient2: str,
    spacing: int = 5_000
) -> Dict:
    """Build construct with specified orientations."""
    builder = SequenceBuilder(filler)
    builder.append_filler(ENHANCER_DOMAIN_START, "upstream_filler")
    builder.append_module(enhancer1, orient1, "enhancer1")
    builder.append_filler(spacing, "inter_enhancer_spacing")
    builder.append_module(enhancer2, orient2, "enhancer2")
    
    used = builder.cursor - ENHANCER_DOMAIN_START
    remaining = (PROMOTER_POS - ENHANCER_DOMAIN_START) - used
    builder.append_filler(remaining, "spacer_to_promoter")
    builder.append_module(promoter, "+", "promoter")
    sequence = builder.finish(CONSTRUCT_LENGTH)
    return {"sequence": sequence, "features": builder.features}


def save_construct(name: str, sequence: str) -> Path:
    CONSTRUCT_DIR.mkdir(parents=True, exist_ok=True)
    path = CONSTRUCT_DIR / f"{name}.fa"
    with open(path, 'w') as f:
        f.write(f">{name}\n")
        f.write(sequence + "\n")
    return path


def main():
    print("=" * 80)
    print("Building Regulatory Grammar Test Constructs")
    print("=" * 80)
    
    # Load modules
    print("\nLoading modules...")
    modules = {
        "HS2": Module("HS2", load_fasta_sequence(ENHANCER_DIR / "HS2_enhancer.fa")),
        "GATA1": Module("GATA1", load_fasta_sequence(ENHANCER_DIR / "GATA1_module.fa")),
        "KLF1": Module("KLF1", load_fasta_sequence(ENHANCER_DIR / "KLF1_module.fa")),
        "TAL1": Module("TAL1", load_fasta_sequence(ENHANCER_DIR / "TAL1_module.fa")),
        "HNF4A": Module("HNF4A", load_fasta_sequence(ENHANCER_DIR / "HNF4A_module.fa")),
        "CTCF": Module("CTCF", load_fasta_sequence(ENHANCER_DIR / "CTCF_module.fa")),
    }
    
    promoters = {
        "HBG1": Module("HBG1", load_fasta_sequence(PROMOTER_DIR / "HBG1_promoter.fa")),
        "ALB": Module("ALB", load_fasta_sequence(PROMOTER_DIR / "ALB_premoter.fa")),
        "CD19": Module("CD19", load_fasta_sequence(PROMOTER_DIR / "CD19_premoter.fa")),
    }
    
    filler = load_filler(FILLER_PATH)
    
    for name, mod in modules.items():
        print(f"  ✓ {name}: {len(mod.sequence)} bp")
    for name, prom in promoters.items():
        print(f"  ✓ {name} promoter: {len(prom.sequence)} bp")
    
    manifest = []
    
    # ========== PART 1: CELL-TYPE SPECIFICITY (27 constructs) ==========
    print("\n" + "=" * 80)
    print("PART 1: Cell-Type Specificity Matrix")
    print("=" * 80)
    
    celltype_combos = [
        # Erythroid promoter (HBG1)
        ("HBG1", "HS2", "K562", "correct"),
        ("HBG1", "GATA1", "K562", "correct"),
        ("HBG1", "KLF1", "K562", "correct"),
        ("HBG1", "TAL1", "K562", "correct"),
        ("HBG1", "HNF4A", "K562", "wrong_enhancer"),
        ("HBG1", "HS2", "HepG2", "wrong_cell"),
        ("HBG1", "HNF4A", "HepG2", "mismatch"),
        ("HBG1", "HS2", "GM12878", "wrong_cell"),
        ("HBG1", "HNF4A", "GM12878", "mismatch"),
        
        # Hepatic promoter (ALB)
        ("ALB", "HNF4A", "HepG2", "correct"),
        ("ALB", "HS2", "HepG2", "wrong_enhancer"),
        ("ALB", "GATA1", "HepG2", "wrong_enhancer"),
        ("ALB", "HNF4A", "K562", "wrong_cell"),
        ("ALB", "HS2", "K562", "mismatch"),
        ("ALB", "HNF4A", "GM12878", "wrong_cell"),
        ("ALB", "HS2", "GM12878", "mismatch"),
        
        # B-cell promoter (CD19)
        ("CD19", "HS2", "GM12878", "wrong_enhancer"),
        ("CD19", "HNF4A", "GM12878", "wrong_enhancer"),
        ("CD19", "GATA1", "GM12878", "wrong_enhancer"),
        ("CD19", "HS2", "K562", "mismatch"),
        ("CD19", "HNF4A", "K562", "mismatch"),
        ("CD19", "HS2", "HepG2", "mismatch"),
        ("CD19", "HNF4A", "HepG2", "mismatch"),
    ]
    
    for prom_name, enh_name, cell_type, label in celltype_combos:
        name = f"CellType_{prom_name}_{enh_name}_{cell_type}"
        print(f"Building {name} ({label})...")
        result = build_celltype_construct(
            promoters[prom_name],
            modules[enh_name],
            filler
        )
        fasta_path = save_construct(name, result["sequence"])
        manifest.append({
            "construct": name,
            "experiment": "cell_type_specificity",
            "promoter": prom_name,
            "enhancer": enh_name,
            "cell_type": cell_type,
            "expected": label,
            "length": len(result["sequence"]),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": result["features"],
        })
    
    print(f"✓ Built {len(celltype_combos)} cell-type constructs")
    
    # ========== PART 2: PAIRWISE COOPERATIVITY (17 constructs) ==========
    print("\n" + "=" * 80)
    print("PART 2: Pairwise Motif Cooperativity")
    print("=" * 80)
    
    # Singles (controls)
    single_enhancers = ["HS2", "GATA1", "KLF1", "TAL1", "HNF4A"]
    for enh_name in single_enhancers:
        name = f"Single_{enh_name}"
        print(f"Building {name}...")
        result = build_pairwise_construct(
            modules[enh_name],
            None,
            promoters["HBG1"],
            filler
        )
        fasta_path = save_construct(name, result["sequence"])
        manifest.append({
            "construct": name,
            "experiment": "pairwise_cooperativity",
            "enhancer1": enh_name,
            "enhancer2": None,
            "expected": "control",
            "length": len(result["sequence"]),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": result["features"],
        })
    
    # Pairs
    cooperativity_pairs = [
        ("HS2", "GATA1", "synergy"),
        ("HS2", "KLF1", "synergy"),
        ("HS2", "TAL1", "synergy"),
        ("GATA1", "KLF1", "synergy"),
        ("GATA1", "TAL1", "synergy"),
        ("KLF1", "TAL1", "synergy"),
        ("HS2", "HNF4A", "interference"),
        ("GATA1", "HNF4A", "interference"),
        ("KLF1", "HNF4A", "interference"),
        ("TAL1", "HNF4A", "interference"),
    ]
    
    for enh1, enh2, expected in cooperativity_pairs:
        name = f"Pair_{enh1}_{enh2}"
        print(f"Building {name} (expect {expected})...")
        result = build_pairwise_construct(
            modules[enh1],
            modules[enh2],
            promoters["HBG1"],
            filler
        )
        fasta_path = save_construct(name, result["sequence"])
        manifest.append({
            "construct": name,
            "experiment": "pairwise_cooperativity",
            "enhancer1": enh1,
            "enhancer2": enh2,
            "expected": expected,
            "length": len(result["sequence"]),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": result["features"],
        })
    
    # CTCF separator tests
    ctcf_tests = [
        ("HS2", "GATA1", "test_insulation"),
        ("HS2", "HNF4A", "control_already_independent"),
    ]
    
    for enh1, enh2, label in ctcf_tests:
        name = f"CTCF_Sep_{enh1}_{enh2}"
        print(f"Building {name}...")
        result = build_pairwise_construct(
            modules[enh1],
            modules[enh2],
            promoters["HBG1"],
            filler,
            separator=modules["CTCF"]
        )
        fasta_path = save_construct(name, result["sequence"])
        manifest.append({
            "construct": name,
            "experiment": "ctcf_separation",
            "enhancer1": enh1,
            "enhancer2": enh2,
            "separator": "CTCF",
            "expected": label,
            "length": len(result["sequence"]),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": result["features"],
        })
    
    print(f"✓ Built {len(single_enhancers) + len(cooperativity_pairs) + len(ctcf_tests)} cooperativity constructs")
    
    # ========== PART 3: SHORT-RANGE SPACING (10 constructs) ==========
    print("\n" + "=" * 80)
    print("PART 3: Short-Range Spacing Screen")
    print("=" * 80)
    
    spacings = [0, 100, 250, 500, 750, 1000, 2000, 3000, 5000, 10000]
    
    for spacing in spacings:
        name = f"Spacing_{spacing}bp_HS2_GATA1"
        print(f"Building {name}...")
        result = build_spacing_construct(
            modules["HS2"],
            modules["GATA1"],
            promoters["HBG1"],
            filler,
            spacing
        )
        fasta_path = save_construct(name, result["sequence"])
        manifest.append({
            "construct": name,
            "experiment": "short_range_spacing",
            "enhancer1": "HS2",
            "enhancer2": "GATA1",
            "spacing": spacing,
            "length": len(result["sequence"]),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": result["features"],
        })
    
    print(f"✓ Built {len(spacings)} spacing constructs")
    
    # ========== PART 4: ORIENTATION EFFECTS (16 constructs) ==========
    print("\n" + "=" * 80)
    print("PART 4: Orientation Effects")
    print("=" * 80)
    
    orientation_pairs = [
        ("HS2", "GATA1", "expect_cooperativity"),
        ("HS2", "HNF4A", "expect_independence"),
        ("GATA1", "KLF1", "expect_cooperativity"),
        ("HS2", "CTCF", "CTCF_orientation_critical"),
    ]
    
    orientations = [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]
    
    for enh1, enh2, label in orientation_pairs:
        for or1, or2 in orientations:
            name = f"Orient_{enh1}{or1}_{enh2}{or2}"
            print(f"Building {name}...")
            result = build_orientation_construct(
                modules[enh1],
                modules[enh2],
                promoters["HBG1"],
                filler,
                or1,
                or2
            )
            fasta_path = save_construct(name, result["sequence"])
            manifest.append({
                "construct": name,
                "experiment": "orientation_effects",
                "enhancer1": enh1,
                "enhancer2": enh2,
                "orientation1": or1,
                "orientation2": or2,
                "expected": label,
                "length": len(result["sequence"]),
                "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
                "features": result["features"],
            })
    
    print(f"✓ Built {len(orientation_pairs) * len(orientations)} orientation constructs")
    
    # ========== SAVE MANIFEST ==========
    print("\n" + "=" * 80)
    print(f"Saving manifest to {MANIFEST_PATH}...")
    with open(MANIFEST_PATH, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total constructs built: {len(manifest)}")
    print(f"  Cell-type specificity: {len([m for m in manifest if m['experiment'] == 'cell_type_specificity'])}")
    print(f"  Pairwise cooperativity: {len([m for m in manifest if m['experiment'] == 'pairwise_cooperativity'])}")
    print(f"  CTCF separation: {len([m for m in manifest if m['experiment'] == 'ctcf_separation'])}")
    print(f"  Short-range spacing: {len([m for m in manifest if m['experiment'] == 'short_range_spacing'])}")
    print(f"  Orientation effects: {len([m for m in manifest if m['experiment'] == 'orientation_effects'])}")
    print("\n✅ Ready for predictions!")
    print("=" * 80)


if __name__ == "__main__":
    main()
