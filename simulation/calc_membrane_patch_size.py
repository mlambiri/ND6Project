from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_EXCLUDE_RESNAMES = {
    # Waters
    "HOH",
    "WAT",
    "TIP",
    "TIP3",
    "TP3",
    # Common ions (3-letter PDB-ish)
    "SOD",
    "CLA",
    "POT",
    "CAL",
    "MG",
    "ZN",
    "NA",
    "CL",
    # Lipids seen in this repo
    "POP",
    "TYC",
    "CDL",
    "PEE",
    "PLX",
    "DGT",
}


@dataclass
class MinMax:
    minx: float = math.inf
    miny: float = math.inf
    minz: float = math.inf
    maxx: float = -math.inf
    maxy: float = -math.inf
    maxz: float = -math.inf

    def update(self, x: float, y: float, z: float) -> None:
        if x < self.minx:
            self.minx = x
        if y < self.miny:
            self.miny = y
        if z < self.minz:
            self.minz = z
        if x > self.maxx:
            self.maxx = x
        if y > self.maxy:
            self.maxy = y
        if z > self.maxz:
            self.maxz = z

    def valid(self) -> bool:
        return self.minx != math.inf

    def dx(self) -> float:
        return self.maxx - self.minx

    def dy(self) -> float:
        return self.maxy - self.miny

    def dz(self) -> float:
        return self.maxz - self.minz


def _split_csv(values: Iterable[str]) -> set[str]:
    out: set[str] = set()
    for item in values:
        for part in str(item).split(","):
            part = part.strip()
            if not part:
                continue
            out.add(part.upper())
    return out


def iter_pdb_atoms(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 54:
                continue
            resname = line[17:20].strip().upper()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            yield resname, x, y, z


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Compute Complex I X/Y extents and recommended membrane patch size "
            "(complex dims + 2*margin)."
        )
    )
    ap.add_argument("--pdb", required=True, help="Input PDB path (system or protein-only).")
    ap.add_argument("--margin", type=float, default=25.0, help="Margin in Å to add on each side (default: 25).")
    ap.add_argument(
        "--exclude-resnames",
        action="append",
        default=[],
        help=(
            "Comma-separated residue names to exclude from the Complex I extent calculation "
            "(can be provided multiple times)."
        ),
    )
    ap.add_argument(
        "--include-lipids",
        action="store_true",
        help="Include lipid atoms in the bounding box (by default they are excluded).",
    )
    args = ap.parse_args()

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        raise SystemExit(f"Missing PDB: {pdb_path}")

    exclude = set(DEFAULT_EXCLUDE_RESNAMES)
    exclude |= _split_csv(args.exclude_resnames)
    if args.include_lipids:
        exclude -= {"POP", "TYC", "CDL", "PEE", "PLX", "DGT"}

    mm = MinMax()
    atoms = 0
    for resname, x, y, z in iter_pdb_atoms(pdb_path):
        if resname in exclude:
            continue
        atoms += 1
        mm.update(x, y, z)

    if atoms < 1 or not mm.valid():
        raise SystemExit(
            "No atoms contributed to the bounding box. "
            "Try --include-lipids or adjust --exclude-resnames."
        )

    dx = mm.dx()
    dy = mm.dy()
    dz = mm.dz()
    patch_x = dx + 2.0 * float(args.margin)
    patch_y = dy + 2.0 * float(args.margin)

    print(f"PDB: {pdb_path}")
    print(f"Included atoms: {atoms}")
    print(f"Complex min (A): x={mm.minx:.3f} y={mm.miny:.3f} z={mm.minz:.3f}")
    print(f"Complex max (A): x={mm.maxx:.3f} y={mm.maxy:.3f} z={mm.maxz:.3f}")
    print(f"Complex extents (A): dx={dx:.3f} dy={dy:.3f} dz={dz:.3f}")
    print(f"Margin (A): {float(args.margin):.3f} per side")
    print(f"Recommended membrane patch (A): x={patch_x:.3f} y={patch_y:.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
