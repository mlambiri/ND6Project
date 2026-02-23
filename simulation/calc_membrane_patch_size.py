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
    "POPC",
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
            chain = line[21:22].strip().upper()
            atomname = line[12:16].strip().upper()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            yield resname, chain, atomname, x, y, z


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(a, a))


def _scale(a: tuple[float, float, float], s: float) -> tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def _normalize(a: tuple[float, float, float]) -> tuple[float, float, float]:
    n = _norm(a)
    if n == 0.0:
        raise ValueError("Cannot normalize zero-length vector")
    return _scale(a, 1.0 / n)


def _jacobi_eigen_3x3(
    a: list[list[float]],
    *,
    max_iter: int = 50,
    eps: float = 1e-12,
) -> tuple[list[float], list[list[float]]]:
    """Eigen-decomposition of a real symmetric 3x3 matrix using Jacobi rotations.

    Returns (eigenvalues, eigenvectors_matrix), where eigenvectors are columns of V.
    """
    v = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]

    def max_offdiag(m: list[list[float]]):
        pairs = [(0, 1), (0, 2), (1, 2)]
        p, q = max(pairs, key=lambda ij: abs(m[ij[0]][ij[1]]))
        return p, q, abs(m[p][q])

    for _ in range(max_iter):
        p, q, off = max_offdiag(a)
        if off < eps:
            break

        app = a[p][p]
        aqq = a[q][q]
        apq = a[p][q]
        if apq == 0.0:
            continue

        phi = 0.5 * math.atan2(2.0 * apq, (aqq - app))
        c = math.cos(phi)
        s = math.sin(phi)

        for i in range(3):
            if i == p or i == q:
                continue
            aip = a[i][p]
            aiq = a[i][q]
            a[i][p] = a[p][i] = c * aip - s * aiq
            a[i][q] = a[q][i] = s * aip + c * aiq

        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        a[p][q] = a[q][p] = 0.0

        for i in range(3):
            vip = v[i][p]
            viq = v[i][q]
            v[i][p] = c * vip - s * viq
            v[i][q] = s * vip + c * viq

    eigvals = [a[0][0], a[1][1], a[2][2]]
    return eigvals, v


def _pca_basis(points: list[tuple[float, float, float]]):
    if len(points) < 3:
        raise ValueError("Need at least 3 points for PCA")

    n = float(len(points))
    mean = (sum(p[0] for p in points) / n, sum(p[1] for p in points) / n, sum(p[2] for p in points) / n)

    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for x, y, z in points:
        dx = x - mean[0]
        dy = y - mean[1]
        dz = z - mean[2]
        cov[0][0] += dx * dx
        cov[0][1] += dx * dy
        cov[0][2] += dx * dz
        cov[1][1] += dy * dy
        cov[1][2] += dy * dz
        cov[2][2] += dz * dz
    cov[0][0] /= n
    cov[0][1] /= n
    cov[0][2] /= n
    cov[1][1] /= n
    cov[1][2] /= n
    cov[2][2] /= n
    cov[1][0] = cov[0][1]
    cov[2][0] = cov[0][2]
    cov[2][1] = cov[1][2]

    eigvals, v = _jacobi_eigen_3x3(cov)
    eig = [
        (eigvals[0], (v[0][0], v[1][0], v[2][0])),
        (eigvals[1], (v[0][1], v[1][1], v[2][1])),
        (eigvals[2], (v[0][2], v[1][2], v[2][2])),
    ]
    eig.sort(key=lambda t: t[0])  # ascending

    nvec = _normalize(eig[0][1])
    uvec = _normalize(eig[2][1])
    if abs(_dot(uvec, nvec)) > 0.9:
        uvec = _normalize(eig[1][1])
    vvec = _normalize(_cross(nvec, uvec))
    uvec = _normalize(_cross(vvec, nvec))

    return mean, uvec, vvec, nvec, [e[0] for e in eig]


def _parse_chain_list(chains: str) -> set[str]:
    parts = [p.strip().upper() for p in chains.replace(",", " ").split() if p.strip()]
    return set("".join(parts))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compute recommended membrane patch size (dims + 2*margin).",
    )
    ap.add_argument("--pdb", required=True, help="Input PDB path (system or protein-only).")
    ap.add_argument("--margin", type=float, default=25.0, help="Margin in A to add on each side (default: 25).")
    ap.add_argument(
        "--mode",
        choices=["nd_plane", "xy"],
        default="nd_plane",
        help=(
            "Sizing mode: 'nd_plane' infers the membrane plane from ND subunits (recommended); "
            "'xy' uses the axis-aligned X/Y bounding box."
        ),
    )
    ap.add_argument(
        "--extent",
        choices=["nd", "protein"],
        default="nd",
        help="When --mode nd_plane, compute extents from ND subunits ('nd') or full protein projection ('protein').",
    )
    ap.add_argument(
        "--nd-chains",
        default="s,i,j,r,l,m",
        help="ND chain IDs used to infer membrane orientation (default: s,i,j,r,l,m).",
    )
    ap.add_argument(
        "--exclude-resnames",
        action="append",
        default=[],
        help=(
            "Comma-separated residue names to exclude from extent calculation (can be provided multiple times). "
            "Default excludes common waters/ions/lipids."
        ),
    )
    ap.add_argument(
        "--include-lipids",
        action="store_true",
        help="Include lipid atoms in extents (by default they are excluded).",
    )
    args = ap.parse_args()

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        raise SystemExit(f"Missing PDB: {pdb_path}")

    exclude = set(DEFAULT_EXCLUDE_RESNAMES)
    exclude |= _split_csv(args.exclude_resnames)
    if args.include_lipids:
        exclude -= {"POP", "POPC", "TYC", "CDL", "PEE", "PLX", "DGT"}

    if args.mode == "xy":
        mm = MinMax()
        atoms = 0
        for resname, _chain, _atomname, x, y, z in iter_pdb_atoms(pdb_path):
            if resname in exclude:
                continue
            atoms += 1
            mm.update(x, y, z)

        if atoms < 1 or not mm.valid():
            raise SystemExit(
                "No atoms contributed to the bounding box. Try --include-lipids or adjust --exclude-resnames."
            )

        dx = mm.dx()
        dy = mm.dy()
        dz = mm.dz()
        patch_x = dx + 2.0 * float(args.margin)
        patch_y = dy + 2.0 * float(args.margin)

        print(f"PDB: {pdb_path}")
        print("Mode: xy (axis-aligned)")
        print(f"Included atoms: {atoms}")
        print(f"Complex min (A): x={mm.minx:.3f} y={mm.miny:.3f} z={mm.minz:.3f}")
        print(f"Complex max (A): x={mm.maxx:.3f} y={mm.maxy:.3f} z={mm.maxz:.3f}")
        print(f"Complex extents (A): dx={dx:.3f} dy={dy:.3f} dz={dz:.3f}")
        print(f"Margin (A): {float(args.margin):.3f} per side")
        print(f"Recommended membrane patch (A): x={patch_x:.3f} y={patch_y:.3f}")
        return 0

    # nd_plane mode
    nd_chains = _parse_chain_list(args.nd_chains)
    nd_ca: list[tuple[float, float, float]] = []
    nd_all: list[tuple[float, float, float]] = []
    prot_all: list[tuple[float, float, float]] = []

    for resname, chain, atomname, x, y, z in iter_pdb_atoms(pdb_path):
        if resname in exclude:
            continue
        prot_all.append((x, y, z))
        if chain in nd_chains:
            nd_all.append((x, y, z))
            if atomname == "CA":
                nd_ca.append((x, y, z))

    if len(nd_ca) < 3:
        raise SystemExit(
            f"Not enough ND CA atoms to infer membrane plane (chains={sorted(c.lower() for c in nd_chains)}; "
            f"found {len(nd_ca)})."
        )

    _mean, u, v, nvec, eigs = _pca_basis(nd_ca)

    points_for_extent = prot_all if args.extent == "protein" else nd_all
    if len(points_for_extent) < 1:
        raise SystemExit("No atoms selected for extent calculation.")

    min_u = min_v = math.inf
    max_u = max_v = -math.inf
    for p in points_for_extent:
        uu = _dot(u, p)
        vv = _dot(v, p)
        min_u = min(min_u, uu)
        max_u = max(max_u, uu)
        min_v = min(min_v, vv)
        max_v = max(max_v, vv)

    du = max_u - min_u
    dv = max_v - min_v
    patch_x = du + 2.0 * float(args.margin)
    patch_y = dv + 2.0 * float(args.margin)

    print(f"PDB: {pdb_path}")
    print("Mode: nd_plane (ND PCA inferred)")
    print(f"ND chains: {','.join(sorted(c.lower() for c in nd_chains))}")
    print(f"ND CA atoms used: {len(nd_ca)}")
    print(f"ND eigvals: {', '.join(f'{x:.6f}' for x in eigs)}")
    print(f"Membrane normal (unit): nx={nvec[0]:.6f} ny={nvec[1]:.6f} nz={nvec[2]:.6f}")
    print(f"Extent selection: {args.extent}")
    print(f"Extents in membrane plane (A): du={du:.3f} dv={dv:.3f}")
    print(f"Margin (A): {float(args.margin):.3f} per side")
    print(f"Recommended membrane patch (A): x={patch_x:.3f} y={patch_y:.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
