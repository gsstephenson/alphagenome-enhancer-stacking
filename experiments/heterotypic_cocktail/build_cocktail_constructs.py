#!/usr/bin/env python3
"""Build heterotypic enhancer cocktail constructs for AlphaGenome."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "heterotypic_cocktail"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

HS2_PATH = ROOT / "sequences" / "enhancers" / "HS2_enhancer.fa"
GATA1_MODULE_PATH = ROOT / "sequences" / "enhancers" / "GATA1_module.fa"
HNF4A_MODULE_PATH = ROOT / "sequences" / "enhancers" / "HNF4A_module.fa"
CTCF_MODULE_PATH = ROOT / "sequences" / "enhancers" / "CTCF_module.fa"
PROMOTER_PATH = ROOT / "sequences" / "promoters" / "HBG1_promoter.fa"
FILLER_PATH = ROOT / "filler" / "1M_filler.txt"

CONSTRUCT_LENGTH = 1_048_576
PROMOTER_POS = 500_000
DOMAIN_START = 250_000


@dataclass
class Module:
    name: str
    sequence: str

    def oriented(self, orientation: str) -> str:
        if orientation == "+":
            return self.sequence
        if orientation == "-":
            return reverse_complement(self.sequence)
        raise ValueError(f"Unknown orientation '{orientation}' for module {self.name}")


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(table)[::-1]


def parse_hs2(path: Path) -> str:
    text = path.read_text()
    match = re.search(r"<DNA[^>]*>(.*?)</DNA>", text, re.DOTALL)
    if not match:
        raise ValueError("HS2 sequence not found in XML FASTA")
    return "".join(match.group(1).split()).upper()


def parse_promoter(path: Path) -> str:
    text = path.read_text()
    match = re.search(r"<DNA[^>]*>(.*?)</DNA>", text, re.DOTALL)
    if not match:
        raise ValueError("Promoter DNA not found in XML FASTA")
    return "".join(match.group(1).split()).upper()


def load_filler(path: Path) -> str:
    filler = path.read_text().strip().upper()
    if not filler:
        raise ValueError("Filler sequence is empty")
    return filler


def load_plain_sequence(path: Path) -> str:
    with open(path, "r") as handle:
        lines = [line.strip() for line in handle if not line.startswith(">")]
    sequence = "".join(lines).upper()
    if not sequence:
        raise ValueError(f"Sequence empty in {path}")
    return sequence


def build_module_definitions(hs2_sequence: str) -> Dict[str, Module]:
    modules = {
        "HS2": Module("HS2", hs2_sequence),
        "GATA1": Module("GATA1", load_plain_sequence(GATA1_MODULE_PATH)),
        "HNF4A": Module("HNF4A", load_plain_sequence(HNF4A_MODULE_PATH)),
        "CTCF": Module("CTCF", load_plain_sequence(CTCF_MODULE_PATH)),
    }
    return modules


@dataclass
class CocktailConfig:
    name: str
    description: str
    module_order: List[str]
    orientation_pattern: List[str]
    module_spacing: int
    repeat_spacing: int
    repeat_count: int
    ctcf_brackets: bool = False
    repeat_separator: Optional[str] = None
    repeat_separator_orientation: str = "+"

    def __post_init__(self) -> None:
        if len(self.module_order) != len(self.orientation_pattern):
            raise ValueError(f"Orientation pattern length must match module order ({self.name})")


class SequenceBuilder:
    def __init__(self, filler: str):
        self._filler = filler
        self._filler_idx = 0
        self.parts: List[str] = []
        self.cursor = 0
        self.features: List[Dict] = []
        self.events: List[Dict] = []

    def _take_filler(self, length: int) -> str:
        if length <= 0:
            return ""
        chunks: List[str] = []
        remaining = length
        while remaining > 0:
            chunk_len = min(len(self._filler) - self._filler_idx, remaining)
            chunks.append(self._filler[self._filler_idx:self._filler_idx + chunk_len])
            self._filler_idx = (self._filler_idx + chunk_len) % len(self._filler)
            remaining -= chunk_len
        return "".join(chunks)

    def append_filler(self, length: int, label: Optional[str] = None, metadata: Optional[Dict] = None) -> None:
        if length <= 0:
            return
        seq = self._take_filler(length)
        self._append_sequence(seq, label=label, metadata=metadata)

    def append_module(self, module: Module, orientation: str, label: str, metadata: Dict) -> None:
        seq = module.oriented(orientation)
        meta = dict(metadata)
        meta.update({"module": module.name, "orientation": orientation})
        self._append_sequence(seq, label=label, metadata=meta)

    def append_sequence(self, seq: str, label: Optional[str] = None, metadata: Optional[Dict] = None) -> None:
        self._append_sequence(seq, label, metadata)

    def _append_sequence(self, seq: str, label: Optional[str], metadata: Optional[Dict]) -> None:
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

    def record_event(self, name: str, metadata: Dict) -> None:
        event = {"event": name, "position": self.cursor}
        event.update(metadata)
        self.events.append(event)

    def finish(self, total_length: int) -> str:
        if self.cursor > total_length:
            raise ValueError(f"Construct exceeds target length ({self.cursor} > {total_length})")
        self.append_filler(total_length - self.cursor)
        return "".join(self.parts)


def build_construct(config: CocktailConfig, modules: Dict[str, Module], promoter: str, filler: str) -> Dict:
    builder = SequenceBuilder(filler)
    builder.append_filler(DOMAIN_START, label="upstream_filler")

    if config.ctcf_brackets:
        builder.append_module(modules["CTCF"], "+", "ctcf_bracket", {"anchor": "left"})
        builder.record_event("ctcf_bracket_added", {"anchor": "left"})

    for repeat_idx in range(config.repeat_count):
        for order_idx, module_name in enumerate(config.module_order):
            module = modules[module_name]
            orientation = config.orientation_pattern[order_idx]
            metadata = {
                "repeat_index": repeat_idx,
                "order_index": order_idx,
            }
            builder.append_module(module, orientation, "enhancer_module", metadata)

            if order_idx < len(config.module_order) - 1:
                builder.append_filler(config.module_spacing, label="module_spacing")

        if config.repeat_separator:
            separator_module = modules[config.repeat_separator]
            builder.append_module(separator_module, config.repeat_separator_orientation, "repeat_separator", {"repeat_index": repeat_idx})

        if repeat_idx < config.repeat_count - 1:
            builder.append_filler(config.repeat_spacing, label="repeat_spacing")

    if config.ctcf_brackets:
        builder.append_module(modules["CTCF"], "-", "ctcf_bracket", {"anchor": "right"})
        builder.record_event("ctcf_bracket_added", {"anchor": "right"})

    if builder.cursor > PROMOTER_POS:
        raise ValueError(f"Construct {config.name} overflows promoter position (cursor={builder.cursor}, promoter={PROMOTER_POS})")

    builder.append_filler(PROMOTER_POS - builder.cursor, label="spacer_to_promoter")
    builder.append_sequence(promoter, label="promoter", metadata={"length": len(promoter)})

    sequence = builder.finish(CONSTRUCT_LENGTH)

    return {
        "sequence": sequence,
        "features": builder.features,
        "events": builder.events,
    }


def save_construct(name: str, sequence: str) -> Path:
    CONSTRUCT_DIR.mkdir(parents=True, exist_ok=True)
    path = CONSTRUCT_DIR / f"{name}_construct.fa"
    with open(path, "w") as handle:
        handle.write(f">{name}_construct\n")
        handle.write(sequence + "\n")
    return path


COCKTAIL_CONFIGS: List[CocktailConfig] = [
    CocktailConfig(
        name="Cocktail_1kbForward",
        description="HS2→GATA1→HNF4A repeated 12x with 1 kb intra-module and 2 kb inter-repeat spacing; CTCF brackets.",
        module_order=["HS2", "GATA1", "HNF4A"],
        orientation_pattern=["+", "+", "+"],
        module_spacing=1_000,
        repeat_spacing=2_000,
        repeat_count=12,
        ctcf_brackets=True,
    ),
    CocktailConfig(
        name="Cocktail_5kbForward",
        description="HS2→GATA1→HNF4A repeated 6x with 5 kb spacing; CTCF brackets.",
        module_order=["HS2", "GATA1", "HNF4A"],
        orientation_pattern=["+", "+", "+"],
        module_spacing=5_000,
        repeat_spacing=5_000,
        repeat_count=6,
        ctcf_brackets=True,
    ),
    CocktailConfig(
        name="Cocktail_20kbForward",
        description="HS2→GATA1→HNF4A repeated 3x with 20 kb spacing; CTCF brackets.",
        module_order=["HS2", "GATA1", "HNF4A"],
        orientation_pattern=["+", "+", "+"],
        module_spacing=20_000,
        repeat_spacing=20_000,
        repeat_count=3,
        ctcf_brackets=True,
    ),
    CocktailConfig(
        name="Cocktail_5kbReverseOrder",
        description="HNF4A→GATA1→HS2 reversed order repeated 6x with 5 kb spacing.",
        module_order=["HNF4A", "GATA1", "HS2"],
        orientation_pattern=["+", "+", "+"],
        module_spacing=5_000,
        repeat_spacing=5_000,
        repeat_count=6,
        ctcf_brackets=True,
    ),
    CocktailConfig(
        name="Cocktail_5kbAltOrientation",
        description="HS2(+), GATA1(-), HNF4A(+) repeated 6x with 5 kb spacing to test polarity sensitivity.",
        module_order=["HS2", "GATA1", "HNF4A"],
        orientation_pattern=["+", "-", "+"],
        module_spacing=5_000,
        repeat_spacing=5_000,
        repeat_count=6,
        ctcf_brackets=True,
    ),
    CocktailConfig(
        name="Cocktail_CTCFSeparated",
        description="HS2→GATA1→HNF4A repeated 5x with 2 kb spacing and CTCF separators between repeats.",
        module_order=["HS2", "GATA1", "HNF4A"],
        orientation_pattern=["+", "+", "+"],
        module_spacing=2_000,
        repeat_spacing=3_000,
        repeat_count=5,
        ctcf_brackets=True,
        repeat_separator="CTCF",
        repeat_separator_orientation="+",
    ),
]


def main() -> None:
    hs2 = parse_hs2(HS2_PATH)
    promoter = parse_promoter(PROMOTER_PATH)
    filler = load_filler(FILLER_PATH)
    modules = build_module_definitions(hs2)

    manifest: List[Dict] = []

    for config in COCKTAIL_CONFIGS:
        print(f"Building {config.name}...")
        result = build_construct(config, modules, promoter, filler)
        fasta_path = save_construct(config.name, result["sequence"])
        manifest.append(
            {
                "construct": config.name,
                "description": config.description,
                "length": len(result["sequence"]),
                "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
                "features": result["features"],
                "events": result["events"],
            }
        )
        print(f"  ✓ Saved {fasta_path.name} ({len(result['sequence']):,} bp)")

    with open(MANIFEST_PATH, "w") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"Manifest written to {MANIFEST_PATH.relative_to(EXPERIMENT_ROOT)}")


if __name__ == "__main__":
    main()
