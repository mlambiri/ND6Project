from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


WATER_RESNAMES = {"HOH", "WAT", "TIP", "TIP3", "TP3"}
ION_RESNAMES = {"SOD", "CLA", "POT", "CAL", "MG", "ZN", "NA", "CL"}

# Lipid residue names seen in this repo (PDB/PSF may use 3- or 4-letter variants).
DEFAULT_LIPID_RESNAMES = {"POP", "POPC", "TYC", "CDL", "PEE", "PLX", "DGT"}

DEFAULT_EXCLUDE_FOR_PROTEIN_BBOX = WATER_RESNAMES | ION_RESNAMES | DEFAULT_LIPID_RESNAMES


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

    def center_xy(self) -> tuple[float, float]:
        return ((self.minx + self.maxx) * 0.5, (self.miny + self.maxy) * 0.5)


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
    """Yield (resname, atomname, x, y, z, line) for ATOM/HETATM records.

    Notes on residue names:
    - Standard PDB has 3-char residue names (cols 18-20).
    - VMD/CHARMM-style PDBs may use 4-char residue names (cols 18-21, chain blank).
    This parser reads up to 4 chars (line[17:21]) and strips whitespace.
    """
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 54:
                continue
            resname = line[17:21].strip().upper()
            atomname = line[12:16].strip().upper()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            yield resname, atomname, x, y, z, line


def compute_bbox(
    pdb_path: Path,
    *,
    exclude_resnames: set[str] | None = None,
    exclude_water: bool = False,
) -> tuple[int, MinMax]:
    mm = MinMax()
    atoms = 0
    for resname, _atomname, x, y, z, _line in iter_pdb_atoms(pdb_path):
        if exclude_resnames and resname in exclude_resnames:
            continue
        if exclude_water and resname in WATER_RESNAMES:
            continue
        atoms += 1
        mm.update(x, y, z)
    return atoms, mm


def compute_midplane_z_from_p_atoms(
    pdb_path: Path,
    *,
    lipid_resnames: set[str] | None = None,
) -> tuple[float, int]:
    mm = MinMax()
    n = 0
    for resname, atomname, x, y, z, _line in iter_pdb_atoms(pdb_path):
        if atomname != "P":
            continue
        if resname in WATER_RESNAMES or resname in ION_RESNAMES:
            continue
        if lipid_resnames is not None and resname not in lipid_resnames:
            continue
        n += 1
        mm.update(x, y, z)
    if n < 1 or not mm.valid():
        raise SystemExit(f"No lipid phosphorus atoms found in: {pdb_path}")
    return (mm.minz + mm.maxz) * 0.5, n


def translate_pdb(
    *,
    in_pdb: Path,
    out_pdb: Path,
    dx: float,
    dy: float,
    dz: float,
    strip_waters: bool,
) -> int:
    written_atoms = 0
    with in_pdb.open("r", encoding="utf-8", errors="ignore") as fin, out_pdb.open(
        "w", encoding="utf-8", newline="\n"
    ) as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if len(line) < 54:
                    continue
                resname = line[17:21].strip().upper()
                if strip_waters and resname in WATER_RESNAMES:
                    continue
                try:
                    x = float(line[30:38]) + dx
                    y = float(line[38:46]) + dy
                    z = float(line[46:54]) + dz
                except ValueError:
                    continue

                # Keep everything else unchanged (occupancy/B-factor/segname/element...).
                fout.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
                if not line.endswith("\n"):
                    fout.write("\n")
                written_atoms += 1
            else:
                fout.write(line)
                if not line.endswith("\n"):
                    fout.write("\n")

    return written_atoms


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Translate a VMD membrane patch PDB into the coordinate frame of a Complex I model.\n\n"
            "XY placement uses the protein bounding-box center.\n"
            "Z placement aligns membrane midplanes using lipid phosphorus (P) atoms.\n"
        )
    )
    ap.add_argument("--protein-pdb", required=True, help="Complex I PDB (protein-only or full system).")
    ap.add_argument("--patch-pdb", required=True, help="Membrane patch PDB produced by VMD membrane plugin.")
    ap.add_argument("--out-pdb", required=True, help="Output placed membrane PDB path.")
    ap.add_argument(
        "--reference-pdb",
        default="output/playwright/chatgpt_botprompts/models/complexI_9TI4_WT_heavy.pdb",
        help=(
            "Reference PDB used to compute the target membrane midplane Z (default: WT heavy model in this repo). "
            "Ignored if --target-midplane-z is provided."
        ),
    )
    ap.add_argument(
        "--target-midplane-z",
        type=float,
        default=None,
        help="If provided, use this Z value (A) as the desired membrane midplane.",
    )
    ap.add_argument(
        "--lipid-resnames",
        action="append",
        default=[],
        help=(
            "Comma-separated lipid residue names to use when computing midplanes from P atoms "
            "(can be provided multiple times). Default uses POP/POPC/TYC/CDL/PEE/PLX/DGT."
        ),
    )
    ap.add_argument(
        "--strip-waters",
        action="store_true",
        help="Write only non-water atoms to --out-pdb (lipids-only).",
    )
    args = ap.parse_args()

    protein_pdb = Path(args.protein_pdb)
    patch_pdb = Path(args.patch_pdb)
    out_pdb = Path(args.out_pdb)

    if not protein_pdb.exists():
        raise SystemExit(f"Missing --protein-pdb: {protein_pdb}")
    if not patch_pdb.exists():
        raise SystemExit(f"Missing --patch-pdb: {patch_pdb}")

    ref_midplane_z: float
    if args.target_midplane_z is not None:
        ref_midplane_z = float(args.target_midplane_z)
        ref_p_atoms = -1
        ref_pdb = None
    else:
        ref_pdb = Path(args.reference_pdb)
        if not ref_pdb.exists():
            raise SystemExit(
                f"Missing --reference-pdb: {ref_pdb}\n"
                "Provide --target-midplane-z or point --reference-pdb at a model that contains native lipids."
            )

        lipid_resnames = set(DEFAULT_LIPID_RESNAMES)
        lipid_resnames |= _split_csv(args.lipid_resnames)

        ref_midplane_z, ref_p_atoms = compute_midplane_z_from_p_atoms(ref_pdb, lipid_resnames=lipid_resnames)

    # Protein bbox center (XY)
    prot_atoms, prot_mm = compute_bbox(protein_pdb, exclude_resnames=DEFAULT_EXCLUDE_FOR_PROTEIN_BBOX)
    if prot_atoms < 1 or not prot_mm.valid():
        raise SystemExit(
            "Protein bounding box is empty. If your PDB is non-standard, try using a protein-only PDB."
        )
    prot_cx, prot_cy = prot_mm.center_xy()

    # Patch bbox center (XY), excluding water so waters don't influence extents.
    patch_atoms, patch_mm = compute_bbox(patch_pdb, exclude_water=True)
    if patch_atoms < 1 or not patch_mm.valid():
        raise SystemExit("Patch bounding box is empty (after excluding waters).")
    patch_cx, patch_cy = patch_mm.center_xy()

    dx = prot_cx - patch_cx
    dy = prot_cy - patch_cy

    # Patch midplane from lipid P atoms (Z)
    patch_midplane_z, patch_p_atoms = compute_midplane_z_from_p_atoms(patch_pdb, lipid_resnames=None)
    dz = ref_midplane_z - patch_midplane_z

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    written = translate_pdb(
        in_pdb=patch_pdb,
        out_pdb=out_pdb,
        dx=dx,
        dy=dy,
        dz=dz,
        strip_waters=bool(args.strip_waters),
    )

    print(f"Protein: {protein_pdb}")
    print(f"Patch:   {patch_pdb}")
    if ref_pdb is not None:
        print(f"Ref:     {ref_pdb}  (P atoms used: {ref_p_atoms})")
    else:
        print(f"Ref:     target-midplane-z={ref_midplane_z:.3f} A")
    print(f"Patch P atoms used: {patch_p_atoms}")
    print(f"Translation (A): dX={dx:.3f} dY={dy:.3f} dZ={dz:.3f}")
    print(f"Wrote: {out_pdb}  (atoms: {written})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
