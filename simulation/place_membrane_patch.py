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

    def center(self) -> tuple[float, float, float]:
        return (
            (self.minx + self.maxx) * 0.5,
            (self.miny + self.maxy) * 0.5,
            (self.minz + self.maxz) * 0.5,
        )

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


def _add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _normalize(a: tuple[float, float, float]) -> tuple[float, float, float]:
    n = _norm(a)
    if n == 0.0:
        raise ValueError("Cannot normalize zero-length vector")
    return _scale(a, 1.0 / n)


def _mat_vec_mul_cols(
    col0: tuple[float, float, float],
    col1: tuple[float, float, float],
    col2: tuple[float, float, float],
    v: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Multiply 3x3 matrix (given by columns) by vector v."""
    x, y, z = v
    return (
        col0[0] * x + col1[0] * y + col2[0] * z,
        col0[1] * x + col1[1] * y + col2[1] * z,
        col0[2] * x + col1[2] * y + col2[2] * z,
    )


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

        # Rotate rows/cols p and q to zero out apq.
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

        # Update eigenvectors
        for i in range(3):
            vip = v[i][p]
            viq = v[i][q]
            v[i][p] = c * vip - s * viq
            v[i][q] = s * vip + c * viq

    eigvals = [a[0][0], a[1][1], a[2][2]]
    return eigvals, v


def _pca_basis(points: list[tuple[float, float, float]]):
    if len(points) < 3:
        raise ValueError("Need at least 3 points for PCA plane fit")

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

    # columns of v are eigenvectors
    eig = [
        (eigvals[0], (v[0][0], v[1][0], v[2][0])),
        (eigvals[1], (v[0][1], v[1][1], v[2][1])),
        (eigvals[2], (v[0][2], v[1][2], v[2][2])),
    ]
    eig.sort(key=lambda t: t[0])  # ascending

    nvec = _normalize(eig[0][1])
    uvec = _normalize(eig[2][1])  # largest variance in-plane axis
    # If numerical issues make u ~ parallel to n, fall back to middle.
    if abs(_dot(uvec, nvec)) > 0.9:
        uvec = _normalize(eig[1][1])

    vvec = _normalize(_cross(nvec, uvec))
    uvec = _normalize(_cross(vvec, nvec))

    return mean, uvec, vvec, nvec, [e[0] for e in eig]


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


def compute_p_points(
    pdb_path: Path,
    *,
    lipid_resnames: set[str] | None = None,
) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for resname, atomname, x, y, z, _line in iter_pdb_atoms(pdb_path):
        if atomname != "P":
            continue
        if resname in WATER_RESNAMES or resname in ION_RESNAMES:
            continue
        if lipid_resnames is not None and resname not in lipid_resnames:
            continue
        points.append((x, y, z))
    return points


def transform_pdb(
    *,
    in_pdb: Path,
    out_pdb: Path,
    r_col0: tuple[float, float, float],
    r_col1: tuple[float, float, float],
    r_col2: tuple[float, float, float],
    t: tuple[float, float, float],
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
                    x0 = float(line[30:38])
                    y0 = float(line[38:46])
                    z0 = float(line[46:54])
                except ValueError:
                    continue

                x1, y1, z1 = _mat_vec_mul_cols(r_col0, r_col1, r_col2, (x0, y0, z0))
                x = x1 + t[0]
                y = y1 + t[1]
                z = z1 + t[2]

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


def _parse_chain_list(chains: str) -> set[str]:
    parts: list[str] = []
    for part in chains.replace(",", " ").split():
        part = part.strip()
        if not part:
            continue
        parts.append(part.lower())
    return set("".join(parts))


def compute_chain_ca_points(pdb_path: Path, *, chains: set[str]) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue
            chain = line[21:22].strip().lower()
            if chain not in chains:
                continue
            atomname = line[12:16].strip().upper()
            if atomname != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            points.append((x, y, z))
    return points


def compute_chain_atom_points(pdb_path: Path, *, chains: set[str]) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue
            chain = line[21:22].strip().lower()
            if chain not in chains:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            points.append((x, y, z))
    return points


def _median(values: list[float]) -> float:
    if not values:
        raise ValueError("median of empty list")
    vs = sorted(values)
    mid = len(vs) // 2
    if len(vs) % 2 == 1:
        return float(vs[mid])
    return 0.5 * (float(vs[mid - 1]) + float(vs[mid]))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Rotate + translate a VMD membrane patch PDB into the coordinate frame of a Complex I model.\n\n"
            "Membrane orientation is inferred from PCA of CA atoms in the ND subunits (default chains: s,i,j,r,l,m).\n"
            "The patch is then centered on the ND CA centroid.\n"
        )
    )
    ap.add_argument("--protein-pdb", required=True, help="Complex I PDB (protein-only or full system).")
    ap.add_argument("--patch-pdb", required=True, help="Membrane patch PDB produced by VMD membrane plugin.")
    ap.add_argument("--out-pdb", required=True, help="Output placed membrane PDB path.")
    ap.add_argument(
        "--strip-waters",
        action="store_true",
        help="Write only non-water atoms to --out-pdb (lipids-only).",
    )
    ap.add_argument(
        "--nd-chains",
        default="s,i,j,r,l,m",
        help="Comma-separated chain IDs for ND1..ND6 used to infer membrane orientation (default: s,i,j,r,l,m).",
    )
    ap.add_argument(
        "--margin",
        type=float,
        default=25.0,
        help="Margin in A per side (used only for reporting recommended patch size).",
    )
    ap.add_argument(
        "--no-rotate",
        action="store_true",
        help="Only translate the patch (do not rotate it to match the reference membrane plane).",
    )
    args = ap.parse_args()

    protein_pdb = Path(args.protein_pdb)
    patch_pdb = Path(args.patch_pdb)
    out_pdb = Path(args.out_pdb)

    if not protein_pdb.exists():
        raise SystemExit(f"Missing --protein-pdb: {protein_pdb}")
    if not patch_pdb.exists():
        raise SystemExit(f"Missing --patch-pdb: {patch_pdb}")

    nd_chains = _parse_chain_list(args.nd_chains)
    nd_ca_points = compute_chain_ca_points(protein_pdb, chains=nd_chains)
    if len(nd_ca_points) < 3:
        raise SystemExit(
            f"Not enough ND CA atoms to infer membrane orientation in {protein_pdb} "
            f"(chains={sorted(nd_chains)}; points={len(nd_ca_points)})."
        )

    nd_center, u_ref, v_ref, n_ref, ref_eigs = _pca_basis(nd_ca_points)

    # Patch: use its lipid phosphorus atoms as the patch "midplane center".
    patch_p_points = compute_p_points(patch_pdb, lipid_resnames=None)
    if len(patch_p_points) < 3:
        raise SystemExit(f"Not enough P atoms found in patch: {patch_pdb} (found {len(patch_p_points)})")
    patch_plane_point, u_patch, v_patch, n_patch, patch_eigs = _pca_basis(patch_p_points)

    if args.no_rotate:
        r0, r1, r2 = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
        rotated_patch_center = patch_plane_point
    else:
        # Rotate patch so its normal matches the ND-derived membrane normal.
        # We build a rotation that maps patch XYZ axes to (u_ref,v_ref,n_ref) in the current coordinate frame.
        r0, r1, r2 = u_ref, v_ref, n_ref
        rotated_patch_center = _mat_vec_mul_cols(r0, r1, r2, patch_plane_point)

    # Center target:
    # - in-plane center uses the ND bounding box center in the inferred membrane plane (u/v),
    #   so the requested margin applies symmetrically.
    # - midplane (along n) uses the median of ND CA projections for robustness.
    nd_all_points = compute_chain_atom_points(protein_pdb, chains=nd_chains)
    mins_u = math.inf
    maxs_u = -math.inf
    mins_v = math.inf
    maxs_v = -math.inf
    for p in nd_all_points:
        uu = _dot(u_ref, p)
        vv = _dot(v_ref, p)
        mins_u = min(mins_u, uu)
        maxs_u = max(maxs_u, uu)
        mins_v = min(mins_v, vv)
        maxs_v = max(maxs_v, vv)

    center_u = 0.5 * (mins_u + maxs_u)
    center_v = 0.5 * (mins_v + maxs_v)
    center_n = _median([_dot(n_ref, p) for p in nd_ca_points])
    target_center = _add(_add(_scale(u_ref, center_u), _scale(v_ref, center_v)), _scale(n_ref, center_n))

    t = _sub(target_center, rotated_patch_center)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    written = transform_pdb(
        in_pdb=patch_pdb,
        out_pdb=out_pdb,
        r_col0=r0,
        r_col1=r1,
        r_col2=r2,
        t=t,
        strip_waters=bool(args.strip_waters),
    )

    # Report recommended patch size in the inferred membrane plane.
    nd_u = maxs_u - mins_u
    nd_v = maxs_v - mins_v
    rec_x = nd_u + 2.0 * float(args.margin)
    rec_y = nd_v + 2.0 * float(args.margin)

    print(f"Protein: {protein_pdb}")
    print(f"Patch:   {patch_pdb}")
    print(f"ND chains:     {','.join(sorted(nd_chains))}  (CA atoms: {len(nd_ca_points)})")
    print(f"ND eigvals:    {', '.join(f'{x:.6f}' for x in ref_eigs)}")
    print(f"Patch P atoms: {len(patch_p_points)}")
    print(f"Patch eigvals: {', '.join(f'{x:.6f}' for x in patch_eigs)}")
    print(f"ND normal:     nx={n_ref[0]:.6f} ny={n_ref[1]:.6f} nz={n_ref[2]:.6f}")
    print(f"Rotation:      {'disabled' if args.no_rotate else 'enabled'}")
    print(f"Translate (A): dX={t[0]:.3f} dY={t[1]:.3f} dZ={t[2]:.3f}")
    print(f"Recommended patch (A): x={rec_x:.3f} y={rec_y:.3f}  (margin={float(args.margin):.1f}/side; ND-only)")
    print(f"Wrote: {out_pdb}  (atoms: {written})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
