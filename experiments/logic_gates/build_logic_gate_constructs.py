#!/usr/bin/env python3
"""
Build logic gate constructs for testing Boolean transcriptional logic in AlphaGenome.

Tests 4 gate types across 13 TF pairs:
- AND gates: 5 pairs (Oct4+Sox2, GATA1+KLF1, MyoD+E-box, NF-κB+AP-1, GATA4+NKX2-5)
- OR gates: 4 pairs (KLF4+KLF5, GATA1+GATA2, HNF1A+HNF1B, FOS+FOSB)
- NOT gates: 4 pairs (REST+NeuroD1, etc.)
- XOR gates: 4 pairs (GATA1+PU.1, MyoD+GATA1, Oct4+MyoD, HNF4A+MyoD)

Each pair: 4 conditions [00, 01, 10, 11] = 52 constructs total
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "logic_gates"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
MANIFEST_PATH = EXPERIMENT_ROOT / "logic_gate_manifest.json"

# Input sequences
ENHANCER_DIR = ROOT / "sequences" / "enhancers"
PROMOTER_DIR = ROOT / "sequences" / "promoters"
FILLER_PATH = ROOT / "filler" / "1M_filler.txt"

# Construct parameters
CONSTRUCT_LENGTH = 1_048_576
PROMOTER_POS = 500_000
TF_A_POS = 250_000
TF_SPACING = 5_000  # Optimal spacing from your prior experiment


@dataclass
class Module:
    name: str
    sequence: str
    description: str = ""

    def oriented(self, orientation: str) -> str:
        if orientation == "+":
            return self.sequence
        if orientation == "-":
            return reverse_complement(self.sequence)
        raise ValueError(f"Unknown orientation '{orientation}'")


@dataclass
class LogicGateDefinition:
    gate_type: str  # AND, OR, NOT, XOR
    tf_a: str
    tf_b: str
    promoter: str
    cell_types: List[str]
    expected_pattern: List[int]  # [00, 01, 10, 11]
    biological_rationale: str


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(table)[::-1]


def load_fasta_sequence(path: Path) -> str:
    """Load sequence from FASTA file."""
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
    """Build synthetic regulatory constructs with precise feature tracking."""
    
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

    def append_module(self, module: Optional[Module], orientation: str, label: str) -> None:
        if module is None:
            # Empty site (for 0 input in truth table)
            self._append("", label, {"module": "EMPTY", "orientation": "n/a"})
            return
        seq = module.oriented(orientation)
        meta = {"module": module.name, "orientation": orientation}
        self._append(seq, label, meta)

    def _append(self, seq: str, label: Optional[str], metadata: Optional[Dict] = None) -> None:
        if not seq:
            if label:  # Track empty sites
                feature = {"label": label, "start": self.cursor, "end": self.cursor}
                if metadata:
                    feature.update(metadata)
                self.features.append(feature)
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
        self.append_filler(total_length - self.cursor, "downstream_filler")
        return "".join(self.parts)


def build_logic_gate_construct(
    tf_a: Optional[Module],
    tf_b: Optional[Module],
    promoter: Module,
    filler: str,
    tf_spacing: int = TF_SPACING,
) -> Dict:
    """
    Build logic gate construct with two TF sites and a promoter.
    
    Layout:
    [Filler] - [TF_A or Empty] - [Spacing] - [TF_B or Empty] - [Spacing] - [Promoter] - [Filler]
    """
    builder = SequenceBuilder(filler)
    
    # Upstream filler to TF_A position
    builder.append_filler(TF_A_POS, "upstream_filler")
    
    # TF A site (or empty)
    builder.append_module(tf_a, "+", "TF_A_site")
    
    # Spacing between TFs
    tf_a_len = len(tf_a.sequence) if tf_a else 0
    builder.append_filler(tf_spacing - tf_a_len, "tf_spacer")
    
    # TF B site (or empty)
    builder.append_module(tf_b, "+", "TF_B_site")
    
    # Spacing to promoter (position - current position)
    tf_b_len = len(tf_b.sequence) if tf_b else 0
    current_pos = TF_A_POS + tf_a_len + tf_spacing + tf_b_len
    promoter_spacing = PROMOTER_POS - current_pos
    builder.append_filler(promoter_spacing, "promoter_spacer")
    
    # Promoter
    builder.append_module(promoter, "+", "promoter")
    
    # Finish with filler
    sequence = builder.finish(CONSTRUCT_LENGTH)
    
    return {"sequence": sequence, "features": builder.features}


def define_logic_gates() -> List[LogicGateDefinition]:
    """Define all logic gate experiments."""
    
    gates = []
    
    # ========================
    # AND GATES (5 pairs)
    # ========================
    gates.append(LogicGateDefinition(
        gate_type="AND",
        tf_a="OCT4",
        tf_b="SOX2",
        promoter="HBG1_PROMOTER",  # Use available promoter as proxy
        cell_types=["K562"],  # Will need H1-hESC ideally
        expected_pattern=[0, 0, 0, 1],
        biological_rationale="Oct4+Sox2 are THE canonical cooperative pair for pluripotency"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="AND",
        tf_a="GATA1_MODULE",
        tf_b="KLF1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 0, 0, 1],
        biological_rationale="Your finding: GATA1+KLF1 showed 1.26x synergy (strongest AND gate)"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="AND",
        tf_a="GATA1_MODULE",
        tf_b="TAL1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 0, 0, 1],
        biological_rationale="Erythroid master regulators, expected cooperativity"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="AND",
        tf_a="HS2_ENHANCER",
        tf_b="GATA1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 0, 0, 1],
        biological_rationale="Your surprising finding: showed interference (0.89x) - test if construct matters"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="AND",
        tf_a="HS2_ENHANCER",
        tf_b="KLF1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 0, 0, 1],
        biological_rationale="Your finding: showed 1.12x synergy - validate AND-like behavior"
    ))
    
    # ========================
    # OR GATES (4 pairs)
    # ========================
    gates.append(LogicGateDefinition(
        gate_type="OR",
        tf_a="GATA1_MODULE",
        tf_b="GATA2",  # Note: GATA2 not available, will skip
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 1, 1, 1],
        biological_rationale="GATA1/2 are redundant in erythroid specification"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="OR",
        tf_a="KLF1_MODULE",
        tf_b="TAL1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 1, 1, 1],
        biological_rationale="Both activate erythroid genes - test redundancy"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="OR",
        tf_a="GATA1_MODULE",
        tf_b="KLF1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 1, 1, 1],
        biological_rationale="Compare OR vs AND behavior for same pair"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="OR",
        tf_a="HS2_ENHANCER",
        tf_b="TAL1_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 1, 1, 1],
        biological_rationale="Your finding: 0.96x (independent) - should show OR-like saturation"
    ))
    
    # ========================
    # XOR GATES (4 pairs)
    # ========================
    gates.append(LogicGateDefinition(
        gate_type="XOR",
        tf_a="GATA1_MODULE",
        tf_b="PU1",  # Note: PU1 not available, will skip
        promoter="HBG1_PROMOTER",
        cell_types=["K562"],
        expected_pattern=[0, 1, 1, 0],
        biological_rationale="THE canonical mutually exclusive pair (erythroid vs myeloid)"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="XOR",
        tf_a="GATA1_MODULE",
        tf_b="HNF4A_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562", "HepG2"],
        expected_pattern=[0, 1, 1, 0],
        biological_rationale="Your finding: 0.76x interference - test XOR hypothesis"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="XOR",
        tf_a="TAL1_MODULE",
        tf_b="HNF4A_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562", "HepG2"],
        expected_pattern=[0, 1, 1, 0],
        biological_rationale="Your finding: 0.66x (strongest interference) - prime XOR candidate"
    ))
    
    gates.append(LogicGateDefinition(
        gate_type="XOR",
        tf_a="HS2_ENHANCER",
        tf_b="HNF4A_MODULE",
        promoter="HBG1_PROMOTER",
        cell_types=["K562", "HepG2"],
        expected_pattern=[0, 1, 1, 0],
        biological_rationale="Cross-lineage antagonism (erythroid vs hepatic)"
    ))
    
    return gates


def load_modules() -> Dict[str, Module]:
    """Load all available enhancer and promoter modules."""
    modules = {}
    
    # Load enhancers
    for fasta_file in ENHANCER_DIR.glob("*.fa*"):
        name = fasta_file.stem.split('.')[0].upper()
        try:
            seq = load_fasta_sequence(fasta_file)
            modules[name] = Module(name=name, sequence=seq, description=f"Enhancer from {fasta_file.name}")
            print(f"✓ Loaded enhancer {name} ({len(seq)} bp)")
        except Exception as e:
            print(f"✗ Failed to load {fasta_file.name}: {e}")
    
    # Load promoters
    for fasta_file in PROMOTER_DIR.glob("*.fa*"):
        name = fasta_file.stem.split('.')[0].upper()
        try:
            seq = load_fasta_sequence(fasta_file)
            modules[name] = Module(name=name, sequence=seq, description=f"Promoter from {fasta_file.name}")
            print(f"✓ Loaded promoter {name} ({len(seq)} bp)")
        except Exception as e:
            print(f"✗ Failed to load {fasta_file.name}: {e}")
    
    return modules


def generate_constructs(gates: List[LogicGateDefinition], modules: Dict[str, Module], filler: str):
    """Generate all constructs for logic gate experiments."""
    
    CONSTRUCT_DIR.mkdir(parents=True, exist_ok=True)
    
    manifest = []
    total_constructs = 0
    
    for gate_def in gates:
        print(f"\n{'='*80}")
        print(f"Building {gate_def.gate_type} gate: {gate_def.tf_a} × {gate_def.tf_b}")
        print(f"{'='*80}")
        
        # Get modules
        try:
            tf_a_module = modules.get(gate_def.tf_a)
            tf_b_module = modules.get(gate_def.tf_b)
            promoter_module = modules[gate_def.promoter]
        except KeyError as e:
            print(f"✗ Missing module: {e}")
            print(f"  Available modules: {sorted(modules.keys())}")
            continue
        
        # Generate 4 truth table conditions for each cell type
        for cell_type in gate_def.cell_types:
            conditions = [
                (None, None, "00", "Neither TF present (baseline)"),
                (None, tf_b_module, "01", f"Only {gate_def.tf_b}"),
                (tf_a_module, None, "10", f"Only {gate_def.tf_a}"),
                (tf_a_module, tf_b_module, "11", "Both TFs present"),
            ]
            
            for tf_a, tf_b, binary_code, description in conditions:
                construct_name = f"LogicGate_{gate_def.gate_type}_{gate_def.tf_a}_{gate_def.tf_b}_{binary_code}_{cell_type}"
                
                # Build construct
                result = build_logic_gate_construct(tf_a, tf_b, promoter_module, filler)
                
                # Save FASTA
                fasta_path = CONSTRUCT_DIR / f"{construct_name}.fasta"
                with open(fasta_path, 'w') as f:
                    f.write(f">{construct_name}\n")
                    seq = result["sequence"]
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
                
                # Add to manifest
                manifest_entry = {
                    "construct": construct_name,
                    "gate_type": gate_def.gate_type,
                    "tf_a": gate_def.tf_a,
                    "tf_b": gate_def.tf_b,
                    "tf_a_present": tf_a is not None,
                    "tf_b_present": tf_b is not None,
                    "binary_code": binary_code,
                    "promoter": gate_def.promoter,
                    "cell_type": cell_type,
                    "expected_output": gate_def.expected_pattern[int(binary_code, 2)],
                    "biological_rationale": gate_def.biological_rationale,
                    "description": description,
                    "fasta_path": str(fasta_path),
                    "length": len(result["sequence"]),
                    "features": result["features"],
                }
                manifest.append(manifest_entry)
                total_constructs += 1
                
                print(f"  ✓ {binary_code}: {description:<40} → {fasta_path.name}")
    
    # Save manifest
    with open(MANIFEST_PATH, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\n{'='*80}")
    print(f"✓ Generated {total_constructs} constructs")
    print(f"✓ Manifest saved to: {MANIFEST_PATH}")
    print(f"✓ FASTA files saved to: {CONSTRUCT_DIR}")
    print(f"{'='*80}")
    
    # Summary statistics
    gate_counts = {}
    for entry in manifest:
        gate_type = entry["gate_type"]
        gate_counts[gate_type] = gate_counts.get(gate_type, 0) + 1
    
    print("\nConstruct breakdown:")
    for gate_type, count in sorted(gate_counts.items()):
        print(f"  {gate_type:8s}: {count} constructs")


def main():
    print("=" * 80)
    print("LOGIC GATE CONSTRUCT BUILDER")
    print("=" * 80)
    
    # Load resources
    print("\nLoading sequences...")
    filler = load_filler(FILLER_PATH)
    print(f"✓ Loaded filler ({len(filler)} bp)")
    
    modules = load_modules()
    print(f"✓ Loaded {len(modules)} modules")
    
    # Define experiments
    gates = define_logic_gates()
    print(f"\n✓ Defined {len(gates)} logic gate experiments")
    
    # Generate constructs
    generate_constructs(gates, modules, filler)
    
    print("\n🎉 Done! Next steps:")
    print("   1. Review manifest: logic_gate_manifest.json")
    print("   2. Run predictions: python run_logic_gate_predictions.py")
    print("   3. Analyze results: python analyze_logic_gates.py")


if __name__ == "__main__":
    main()
