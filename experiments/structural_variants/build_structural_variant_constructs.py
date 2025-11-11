#!/usr/bin/env python3
"""Build structural-variant enhancer constructs for AlphaGenome tests."""

import json
import re
from pathlib import Path
from typing import Callable, Dict, List

ROOT = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/AlphaGenome_EnhancerStacking")
EXPERIMENT_ROOT = ROOT / "experiments" / "structural_variants"
CONSTRUCT_DIR = EXPERIMENT_ROOT / "sequences"
MANIFEST_PATH = EXPERIMENT_ROOT / "construct_manifest.json"

ENHANCER_FILE = ROOT / "sequences" / "enhancers" / "HS2_enhancer.fa"
PROMOTER_FILE = ROOT / "sequences" / "promoters" / "HBG1_promoter.fa"
FILLER_FILE = ROOT / "filler" / "1M_filler.txt"

CONSTRUCT_LENGTH = 1_048_576
PROMOTER_POS = 500_000
ENHANCER_POS = 400_000
ANCHOR_LEFT_POS = 350_000
ANCHOR_RIGHT_POS = 450_000
RELOCATED_ENHANCER_POS = 800_000
ENHANCER_COPIES = 10

CTCF_MOTIF = "CCGCGTGGTGGCAGGAGC"  # High-affinity CTCF consensus (forward)
CTCF_LEN = len(CTCF_MOTIF)


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


CTCF_MOTIF_RC = reverse_complement(CTCF_MOTIF)


def parse_xml_fasta(path: Path) -> str:
    """Extract DNA sequence from DAS-style XML FASTA file."""
    content = path.read_text()
    match = re.search(r"<DNA[^>]*>(.*?)</DNA>", content, re.DOTALL)
    if not match:
        raise ValueError(f"No <DNA> block found in {path}")
    return "".join(match.group(1).split()).upper()


def load_filler(path: Path) -> str:
    filler = path.read_text().strip().upper()
    if not filler:
        raise ValueError("Filler sequence is empty")
    return filler


class SequenceBuilder:
    """Incrementally assemble a construct while tracking annotations."""

    def __init__(self, filler: str):
        self._parts: List[str] = []
        self._cursor = 0
        self._filler = filler
        self._filler_idx = 0
        self.features: List[Dict] = []
        self.events: List[Dict] = []

    @property
    def position(self) -> int:
        return self._cursor

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

    def append_sequence(self, sequence: str, label: str = None, metadata: Dict = None) -> None:
        if not sequence:
            return
        start = self._cursor
        self._parts.append(sequence)
        self._cursor += len(sequence)
        if label:
            feature = {"label": label, "start": start, "end": self._cursor}
            if metadata:
                feature.update(metadata)
            self.features.append(feature)

    def append_filler(self, length: int, label: str = None, metadata: Dict = None) -> None:
        if length <= 0:
            return
        self.append_sequence(self._take_filler(length), label=label, metadata=metadata)

    def append_to(self, target_position: int) -> None:
        if target_position < self._cursor:
            raise ValueError(f"Cannot move backwards (target={target_position}, cursor={self._cursor})")
        self.append_filler(target_position - self._cursor)

    def append_ctcf(self, orientation: str, anchor_label: str) -> None:
        motif = CTCF_MOTIF if orientation == "forward" else CTCF_MOTIF_RC
        self.append_sequence(motif, label="ctcf_anchor", metadata={"anchor": anchor_label, "orientation": orientation})

    def append_enhancer_block(self, enhancer: str, copies: int, block_label: str) -> None:
        block = enhancer * copies
        self.append_sequence(block, label=block_label, metadata={"copies": copies, "unit_length": len(enhancer)})

    def append_promoter(self, promoter: str) -> None:
        self.append_sequence(promoter, label="promoter", metadata={"length": len(promoter)})

    def record_event(self, name: str, metadata: Dict) -> None:
        event = {"event": name, "position": self._cursor}
        event.update(metadata)
        self.events.append(event)

    def finish(self, total_length: int) -> str:
        if self._cursor > total_length:
            raise ValueError(f"Construct exceeds target length ({self._cursor} > {total_length})")
        self.append_filler(total_length - self._cursor)
        return "".join(self._parts)


def build_loop_intact(builder: SequenceBuilder, enhancer: str, promoter: str) -> None:
    builder.append_to(ANCHOR_LEFT_POS)
    builder.append_ctcf("forward", "left")
    builder.append_to(ENHANCER_POS)
    builder.append_enhancer_block(enhancer, ENHANCER_COPIES, "hs2_block")
    builder.append_to(ANCHOR_RIGHT_POS)
    builder.append_ctcf("reverse", "right")
    spacer = PROMOTER_POS - builder.position
    if spacer > 0:
        builder.record_event("spacer_retained", {"length": spacer})
        builder.append_filler(spacer, label="loop_spacer")
    builder.append_promoter(promoter)


def build_loop_inverted(builder: SequenceBuilder, enhancer: str, promoter: str) -> None:
    builder.append_to(ANCHOR_LEFT_POS)
    builder.append_ctcf("forward", "left")
    builder.append_to(ENHANCER_POS)
    builder.append_enhancer_block(enhancer, ENHANCER_COPIES, "hs2_block")
    builder.append_to(ANCHOR_RIGHT_POS)
    builder.append_ctcf("forward", "right")
    spacer = PROMOTER_POS - builder.position
    if spacer > 0:
        builder.record_event("spacer_retained", {"length": spacer})
        builder.append_filler(spacer, label="loop_spacer")
    builder.append_promoter(promoter)
    builder.record_event("anchor_inverted", {"anchor": "right"})


def build_loop_deleted(builder: SequenceBuilder, enhancer: str, promoter: str) -> None:
    builder.append_to(ANCHOR_LEFT_POS)
    builder.append_ctcf("forward", "left")
    builder.append_to(ENHANCER_POS)
    builder.append_enhancer_block(enhancer, ENHANCER_COPIES, "hs2_block")
    builder.append_to(ANCHOR_RIGHT_POS)
    builder.append_ctcf("reverse", "right")
    deleted = max(PROMOTER_POS - builder.position, 0)
    if deleted > 0:
        builder.record_event("spacer_deleted", {"length": deleted})
    builder.append_promoter(promoter)


def build_loop_relocated(builder: SequenceBuilder, enhancer: str, promoter: str) -> None:
    builder.append_to(ANCHOR_LEFT_POS)
    builder.append_ctcf("forward", "left")
    builder.append_to(ENHANCER_POS)
    builder.record_event("enhancer_removed_from_loop", {"copies": ENHANCER_COPIES})
    builder.append_to(ANCHOR_RIGHT_POS)
    builder.append_ctcf("reverse", "right")
    builder.append_to(PROMOTER_POS)
    builder.append_promoter(promoter)
    builder.append_to(RELOCATED_ENHANCER_POS)
    builder.append_enhancer_block(enhancer, ENHANCER_COPIES, "hs2_block")
    builder.record_event("enhancer_relocated", {"to": RELOCATED_ENHANCER_POS, "copies": ENHANCER_COPIES})


VariantBuilder = Callable[[SequenceBuilder, str, str], None]

VARIANTS: List[Dict] = [
    {
        "name": "LoopIntact_10x",
        "description": "10× HS2 array with convergent CTCF anchors at 350/450 kb; promoter at 500 kb.",
        "builder": build_loop_intact,
    },
    {
        "name": "LoopInverted_10x",
        "description": "Right CTCF anchor flipped to disrupt loop polarity while keeping enhancer position constant.",
        "builder": build_loop_inverted,
    },
    {
        "name": "LoopDeleted_10x",
        "description": "Spacer between right anchor and promoter removed so promoter abuts the looped enhancer block.",
        "builder": build_loop_deleted,
    },
    {
        "name": "LoopRelocated_10x",
        "description": "Enhancer block moved to 800 kb; anchors remain around an empty loop near the promoter.",
        "builder": build_loop_relocated,
    },
]


def save_fasta(name: str, sequence: str, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    fasta_path = directory / f"{name}_construct.fa"
    with open(fasta_path, "w") as handle:
        handle.write(f">{name}_construct\n")
        handle.write(sequence + "\n")
    return fasta_path


def main() -> None:
    enhancer = parse_xml_fasta(ENHANCER_FILE)
    promoter = parse_xml_fasta(PROMOTER_FILE)
    filler = load_filler(FILLER_FILE)

    manifest: List[Dict] = []

    for variant in VARIANTS:
        print(f"Building {variant['name']}...")
        builder = SequenceBuilder(filler)
        variant_builder: VariantBuilder = variant["builder"]
        variant_builder(builder, enhancer, promoter)
        sequence = builder.finish(CONSTRUCT_LENGTH)
        fasta_path = save_fasta(variant["name"], sequence, CONSTRUCT_DIR)
        entry = {
            "construct": variant["name"],
            "description": variant["description"],
            "length": len(sequence),
            "fasta": str(fasta_path.relative_to(EXPERIMENT_ROOT)),
            "features": builder.features,
            "events": builder.events,
        }
        manifest.append(entry)
        print(f"  ✓ Saved {fasta_path.name} ({len(sequence):,} bp)")

    with open(MANIFEST_PATH, "w") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"\nManifest written to {MANIFEST_PATH.relative_to(EXPERIMENT_ROOT)}")
    print("Construct build complete.")


if __name__ == "__main__":
    main()
