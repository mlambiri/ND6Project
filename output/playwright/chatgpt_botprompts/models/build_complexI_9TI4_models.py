#!/usr/bin/env python3
"""
Build "clean" starting coordinate files for human mitochondrial Complex I from RCSB PDB 9TI4.

Outputs (by default):
  - complexI_9TI4_WT_heavy.pdb
  - complexI_9TI4_ND1_A52T_heavy.pdb
  - complexI_9TI4_ND6_M64V_heavy.pdb
  - nd1_chain_s_WT_heavy.pdb
  - nd1_chain_s_A52T_heavy.pdb
  - nd6_chain_m_WT_heavy.pdb
  - nd6_chain_m_M64V_heavy.pdb

Notes
-----
* This script intentionally exports *heavy atoms only* (no hydrogens) and drops HOH,
  to avoid protonation/hydrogen-placement issues before CHARMM/psfgen builds.
* The LHON variant m.14484T>C in MT-ND6 maps to p.Met64Val (M64V). In 9TI4 this is:
    auth_asym_id (chain) = "m", auth_seq_id (residue) = 64, MET -> VAL
* Side-chain coordinates are approximate for the new VAL CG2; downstream minimization
  (or rebuilding the residue with psfgen + guesscoord) is recommended.
* The LHON variant m.3460G>A in MT-ND1 maps to p.Ala52Thr (A52T). In 9TI4 this is:
    auth_asym_id (chain) = "s", auth_seq_id (residue) = 52, ALA -> THR
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


ALLOWED_CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


@dataclass(frozen=True)
class AtomRecord:
    group: str  # ATOM / HETATM
    atom_id: int  # _atom_site.id
    element: str  # _atom_site.type_symbol
    resname: str  # _atom_site.auth_comp_id
    chain_auth: str  # _atom_site.auth_asym_id
    chain_pdb: str  # mapped 1-char chain id for PDB
    segid: str  # 1-4 chars (CHARMM/VMD extension)
    resseq: int  # _atom_site.auth_seq_id
    icode: str  # _atom_site.pdbx_PDB_ins_code ('' if missing)
    atomname: str  # _atom_site.auth_atom_id
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float


def _tokenize_cif_row(line: str) -> List[str]:
    # CIF quoting is only relevant if the token *starts* with a quote.
    tokens: List[str] = []
    i = 0
    n = len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            tokens.append(line[start:i])
            i += 1  # skip closing quote
            continue
        start = i
        while i < n and not line[i].isspace():
            i += 1
        tokens.append(line[start:i])
    return tokens


def _iter_atom_site_loop(cif_path: Path) -> Tuple[List[str], Iterator[List[str]]]:
    """
    Returns (headers, rows_iterator) for the first _atom_site loop.
    """
    f = cif_path.open("r", encoding="utf-8", errors="replace")

    def _rows() -> Iterator[List[str]]:
        in_loop = False
        headers: List[str] = []
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line == "loop_":
                in_loop = True
                headers = []
                continue
            if not in_loop:
                continue
            if line.startswith("_atom_site."):
                headers.append(line)
                continue
            if headers and headers[0].startswith("_atom_site."):
                # Data section for the atom_site loop.
                if line.startswith("#"):
                    break
                yield _tokenize_cif_row(line)
        f.close()

    # We need to read headers first; easiest is to scan once and then re-open for rows.
    headers: List[str] = []
    with cif_path.open("r", encoding="utf-8", errors="replace") as scan:
        in_loop = False
        candidate: List[str] = []
        for raw in scan:
            line = raw.strip()
            if not line:
                continue
            if line == "loop_":
                in_loop = True
                candidate = []
                continue
            if not in_loop:
                continue
            if line.startswith("_atom_site."):
                candidate.append(line)
                continue
            if candidate:
                headers = candidate
                break
    if not headers:
        raise RuntimeError(f"Could not find _atom_site loop in {cif_path}")
    return headers, _rows()


def _altloc_rank(altloc: str) -> int:
    # Prefer 'A' > '' > others
    if altloc == "A":
        return 3
    if altloc == "":
        return 2
    if altloc == ".":
        return 2
    if altloc == "?":
        return 1
    return 0


def _safe_float(value: str, default: float) -> float:
    if value in (".", "?"):
        return default
    return float(value)


def _safe_int(value: str) -> int:
    if value in (".", "?"):
        raise ValueError("Missing integer value")
    return int(value)


def _build_chain_map(chain_ids: List[str]) -> Dict[str, str]:
    used = set()
    mapping: Dict[str, str] = {}
    for auth in chain_ids:
        if auth in mapping:
            continue
        if len(auth) == 1 and auth in ALLOWED_CHAIN_IDS and auth not in used:
            mapping[auth] = auth
            used.add(auth)
            continue
        for candidate in ALLOWED_CHAIN_IDS:
            if candidate not in used:
                mapping[auth] = candidate
                used.add(candidate)
                break
        else:
            raise RuntimeError(
                "Ran out of PDB chain IDs; consider exporting mmCIF instead of PDB."
            )
    return mapping


def parse_cif_atoms(
    cif_path: Path,
    *,
    heavy_only: bool = True,
    drop_hoh: bool = True,
    model_num: int = 1,
) -> Tuple[List[AtomRecord], Dict[str, str]]:
    headers, rows = _iter_atom_site_loop(cif_path)
    col = {h: i for i, h in enumerate(headers)}

    required = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_alt_id",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_PDB_model_num",
    ]
    missing = [h for h in required if h not in col]
    if missing:
        raise RuntimeError(f"Missing expected _atom_site columns: {missing}")

    # First pass: choose one altloc per atom key (highest occupancy).
    chosen: Dict[Tuple[str, str, int, str, str, str], AtomRecord] = {}

    # Keep chain appearance order (for mapping).
    chain_order: List[str] = []
    seen_chains = set()

    for fields in rows:
        # Skip malformed lines.
        if len(fields) < len(headers):
            continue

        pdb_model = _safe_int(fields[col["_atom_site.pdbx_PDB_model_num"]])
        if pdb_model != model_num:
            continue

        element = fields[col["_atom_site.type_symbol"]].upper()
        if heavy_only and element == "H":
            continue

        resname = fields[col["_atom_site.auth_comp_id"]].upper()
        if drop_hoh and resname == "HOH":
            continue

        group = fields[col["_atom_site.group_PDB"]].upper()
        atom_id = _safe_int(fields[col["_atom_site.id"]])
        altloc_raw = fields[col["_atom_site.label_alt_id"]]
        altloc = "" if altloc_raw in (".", "?") else altloc_raw
        resseq = _safe_int(fields[col["_atom_site.auth_seq_id"]])
        chain_auth = fields[col["_atom_site.auth_asym_id"]]
        atomname = fields[col["_atom_site.auth_atom_id"]]
        icode_raw = fields[col["_atom_site.pdbx_PDB_ins_code"]]
        icode = "" if icode_raw in (".", "?") else icode_raw
        x = _safe_float(fields[col["_atom_site.Cartn_x"]], default=0.0)
        y = _safe_float(fields[col["_atom_site.Cartn_y"]], default=0.0)
        z = _safe_float(fields[col["_atom_site.Cartn_z"]], default=0.0)
        occupancy = _safe_float(fields[col["_atom_site.occupancy"]], default=1.0)
        bfactor = _safe_float(fields[col["_atom_site.B_iso_or_equiv"]], default=0.0)

        if chain_auth not in seen_chains:
            chain_order.append(chain_auth)
            seen_chains.add(chain_auth)

        # Key that identifies the same atom across altlocs.
        key = (group, chain_auth, resseq, icode, resname, atomname)

        # chain_pdb assigned later via mapping; segid derived from auth id.
        record = AtomRecord(
            group=group,
            atom_id=atom_id,
            element=element,
            resname=resname,
            chain_auth=chain_auth,
            chain_pdb="?",
            segid=chain_auth[:4],
            resseq=resseq,
            icode=icode,
            atomname=atomname,
            x=x,
            y=y,
            z=z,
            occupancy=occupancy,
            bfactor=bfactor,
        )

        best = chosen.get(key)
        if best is None:
            chosen[key] = record
            continue
        if record.occupancy > best.occupancy + 1e-6:
            chosen[key] = record
            continue
        if abs(record.occupancy - best.occupancy) <= 1e-6:
            if _altloc_rank(altloc) > _altloc_rank(""):
                chosen[key] = record

    # Build chain id mapping for PDB output.
    chain_map = _build_chain_map(chain_order)

    records = [
        AtomRecord(
            **{
                **r.__dict__,
                "chain_pdb": chain_map[r.chain_auth],
            }
        )
        for r in chosen.values()
    ]
    records.sort(key=lambda r: r.atom_id)
    return records, chain_map


def _format_atom_name(atomname: str, element: str) -> str:
    name = atomname.strip()
    if len(name) >= 4:
        return name[:4]
    if not name:
        return "    "
    el = element.strip().upper()
    if name[0].isdigit():
        return name.rjust(4)
    if len(el) == 1 and name[0].upper() == el:
        return name.rjust(4)
    return name.ljust(4)


def format_pdb_line(serial: int, r: AtomRecord) -> str:
    record = "HETATM" if r.group == "HETATM" else "ATOM"
    atomname = _format_atom_name(r.atomname, r.element)
    altloc = " "
    resname = (r.resname or "UNK")[:3].rjust(3)
    chain = (r.chain_pdb or " ")[:1]
    resseq = r.resseq
    icode = (r.icode or " ")[:1]
    segid = (r.segid or "")[:4].ljust(4)
    element = (r.element or "").strip().upper()[:2].rjust(2)

    # PDB columns (CHARMM/VMD-friendly; places segid in cols 73-76 and element in cols 77-78)
    return (
        f"{record:<6}{serial:5d} "
        f"{atomname}{altloc}{resname} {chain}{resseq:4d}{icode}   "
        f"{r.x:8.3f}{r.y:8.3f}{r.z:8.3f}"
        f"{r.occupancy:6.2f}{r.bfactor:6.2f}"
        f"{'':6}{segid}{element}{'':2}"
    )


def write_pdb(path: Path, records: Iterable[AtomRecord]) -> None:
    write_pdb_with_remarks(path, records, remarks=None)


def _format_remark_lines(text: str, remark_num: int = 1) -> List[str]:
    prefix = f"REMARK {remark_num:3d} "
    width = 80 - len(prefix)
    if width < 10:
        raise RuntimeError("Unexpected PDB remark prefix width")

    # Keep URLs intact by not splitting tokens unless they exceed the width.
    words = text.split()
    if not words:
        return [prefix.rstrip().ljust(80)]

    lines: List[str] = []
    current = ""
    for w in words:
        if not current:
            if len(w) <= width:
                current = w
            else:
                # Hard-split long tokens (rare; e.g., very long URLs).
                for i in range(0, len(w), width):
                    chunk = w[i : i + width]
                    lines.append((prefix + chunk)[:80].ljust(80))
                current = ""
            continue

        if len(current) + 1 + len(w) <= width:
            current = f"{current} {w}"
        else:
            lines.append((prefix + current)[:80].ljust(80))
            if len(w) <= width:
                current = w
            else:
                for i in range(0, len(w), width):
                    chunk = w[i : i + width]
                    lines.append((prefix + chunk)[:80].ljust(80))
                current = ""

    if current:
        lines.append((prefix + current)[:80].ljust(80))
    return lines


def write_pdb_with_remarks(
    path: Path, records: Iterable[AtomRecord], *, remarks: Optional[List[str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        if remarks:
            for r in remarks:
                for line in _format_remark_lines(r):
                    out.write(line)
                    out.write("\n")
        serial = 1
        for r in records:
            out.write(format_pdb_line(serial, r))
            out.write("\n")
            serial += 1
        out.write("END\n")


def _vec_sub(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_add(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_scale(a: Tuple[float, float, float], s: float) -> Tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def _vec_dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vec_cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vec_norm(a: Tuple[float, float, float]) -> float:
    return math.sqrt(_vec_dot(a, a))


def _vec_unit(a: Tuple[float, float, float]) -> Tuple[float, float, float]:
    n = _vec_norm(a)
    if n < 1e-12:
        raise ValueError("Zero-length vector")
    return _vec_scale(a, 1.0 / n)


def _rodrigues_rotate(
    v: Tuple[float, float, float],
    axis_unit: Tuple[float, float, float],
    angle_rad: float,
) -> Tuple[float, float, float]:
    # v_rot = v*cos + (k x v)*sin + k*(k·v)*(1-cos)
    cos_t = math.cos(angle_rad)
    sin_t = math.sin(angle_rad)
    k = axis_unit
    kxv = _vec_cross(k, v)
    kdotv = _vec_dot(k, v)
    term1 = _vec_scale(v, cos_t)
    term2 = _vec_scale(kxv, sin_t)
    term3 = _vec_scale(k, kdotv * (1.0 - cos_t))
    return _vec_add(_vec_add(term1, term2), term3)


def mutate_nd6_m64v(
    records: List[AtomRecord],
    *,
    chain_auth: str = "m",
    resseq: int = 64,
    icode: str = "",
) -> List[AtomRecord]:
    """
    Apply MET->VAL at (chain_auth, resseq, icode) using heavy atoms only.
    """
    target = [r for r in records if r.chain_auth == chain_auth and r.resseq == resseq and r.icode == icode]
    if not target:
        raise RuntimeError(f"Could not find target residue {chain_auth}:{resseq}{icode} in records")

    # Collect coordinates for CA/CB/CG (MET).
    atom_by_name = {r.atomname.strip().upper(): r for r in target}
    if "CA" not in atom_by_name or "CB" not in atom_by_name or "CG" not in atom_by_name:
        # Still do a minimal rename/remove; psfgen can rebuild missing atoms later.
        mutated: List[AtomRecord] = []
        for r in records:
            if r.chain_auth != chain_auth or r.resseq != resseq or r.icode != icode:
                mutated.append(r)
                continue
            if r.atomname.strip().upper() in {"SD", "CE"}:
                continue
            atomname = "CG1" if r.atomname.strip().upper() == "CG" else r.atomname
            mutated.append(
                AtomRecord(
                    **{**r.__dict__, "resname": "VAL", "atomname": atomname}
                )
            )
        return mutated

    ca = atom_by_name["CA"]
    cb = atom_by_name["CB"]
    cg = atom_by_name["CG"]

    cb_xyz = (cb.x, cb.y, cb.z)
    ca_xyz = (ca.x, ca.y, ca.z)
    cg_xyz = (cg.x, cg.y, cg.z)

    axis = _vec_sub(ca_xyz, cb_xyz)  # CB -> CA
    v = _vec_sub(cg_xyz, cb_xyz)  # CB -> CG
    axis_unit = _vec_unit(axis)

    # Approximate VAL geometry: rotate the MET CB->CG vector around CB->CA.
    v_rot = _rodrigues_rotate(v, axis_unit, angle_rad=2.0 * math.pi / 3.0)
    cg2_xyz = _vec_add(cb_xyz, v_rot)

    cg2 = AtomRecord(
        group="ATOM",
        atom_id=max(r.atom_id for r in records) + 1,
        element="C",
        resname="VAL",
        chain_auth=cg.chain_auth,
        chain_pdb=cg.chain_pdb,
        segid=cg.segid,
        resseq=resseq,
        icode=icode,
        atomname="CG2",
        x=cg2_xyz[0],
        y=cg2_xyz[1],
        z=cg2_xyz[2],
        occupancy=cg.occupancy,
        bfactor=cg.bfactor,
    )

    mutated_records: List[AtomRecord] = []
    inserted_cg2 = False
    for r in records:
        if r.chain_auth != chain_auth or r.resseq != resseq or r.icode != icode:
            mutated_records.append(r)
            continue

        atom = r.atomname.strip().upper()
        if atom in {"SD", "CE"}:
            continue
        if atom == "CG":
            mutated_records.append(
                AtomRecord(**{**r.__dict__, "resname": "VAL", "atomname": "CG1"})
            )
            continue

        mutated_records.append(AtomRecord(**{**r.__dict__, "resname": "VAL"}))

        # Insert CG2 right after CB if possible; if CB missing, it lands after CA/C/O etc.
        if not inserted_cg2 and atom == "CB":
            mutated_records.append(cg2)
            inserted_cg2 = True

    if not inserted_cg2:
        mutated_records.append(cg2)

    return mutated_records


def mutate_nd1_a52t(
    records: List[AtomRecord],
    *,
    chain_auth: str = "s",
    resseq: int = 52,
    icode: str = "",
) -> List[AtomRecord]:
    """
    Apply ALA->THR at (chain_auth, resseq, icode) using heavy atoms only.

    Adds OG1 and CG2 atoms with approximate tetrahedral geometry around CB.
    Downstream minimization (or rebuilding the residue with psfgen + guesscoord) is recommended.
    """
    target = [r for r in records if r.chain_auth == chain_auth and r.resseq == resseq and r.icode == icode]
    if not target:
        raise RuntimeError(f"Could not find target residue {chain_auth}:{resseq}{icode} in records")

    atom_by_name = {r.atomname.strip().upper(): r for r in target}
    if "CA" not in atom_by_name or "CB" not in atom_by_name:
        mutated: List[AtomRecord] = []
        for r in records:
            if r.chain_auth != chain_auth or r.resseq != resseq or r.icode != icode:
                mutated.append(r)
                continue
            mutated.append(AtomRecord(**{**r.__dict__, "resname": "THR"}))
        return mutated

    ca = atom_by_name["CA"]
    cb = atom_by_name["CB"]
    ca_xyz = (ca.x, ca.y, ca.z)
    cb_xyz = (cb.x, cb.y, cb.z)

    axis = _vec_sub(ca_xyz, cb_xyz)  # CB -> CA
    axis_unit = _vec_unit(axis)

    # Choose a reference direction not parallel to the CB->CA axis.
    ref: Optional[Tuple[float, float, float]] = None
    for cand in ("N", "C", "O"):
        if cand in atom_by_name:
            a = atom_by_name[cand]
            ref = _vec_sub((a.x, a.y, a.z), ca_xyz)  # CA -> cand
            break
    if ref is None:
        ref = (1.0, 0.0, 0.0) if abs(axis_unit[0]) < 0.9 else (0.0, 1.0, 0.0)

    perp = _vec_cross(axis_unit, ref)
    if _vec_norm(perp) < 1e-6:
        ref = (0.0, 1.0, 0.0) if abs(axis_unit[1]) < 0.9 else (0.0, 0.0, 1.0)
        perp = _vec_cross(axis_unit, ref)
    perp1 = _vec_unit(perp)
    perp2 = _rodrigues_rotate(perp1, axis_unit, angle_rad=2.0 * math.pi / 3.0)

    # Tetrahedral directions: cos(109.47) = -1/3.
    along = -1.0 / 3.0
    perp_scale = 2.0 * math.sqrt(2.0) / 3.0

    dir_og1 = _vec_add(_vec_scale(axis_unit, along), _vec_scale(perp1, perp_scale))
    dir_cg2 = _vec_add(_vec_scale(axis_unit, along), _vec_scale(perp2, perp_scale))

    # Approx bond lengths (heavy atoms).
    cb_og = 1.43  # C-O
    cb_cg = 1.53  # C-C

    og1_xyz = _vec_add(cb_xyz, _vec_scale(dir_og1, cb_og))
    cg2_xyz = _vec_add(cb_xyz, _vec_scale(dir_cg2, cb_cg))

    next_id = max(r.atom_id for r in records) + 1
    og1 = AtomRecord(
        group="ATOM",
        atom_id=next_id,
        element="O",
        resname="THR",
        chain_auth=cb.chain_auth,
        chain_pdb=cb.chain_pdb,
        segid=cb.segid,
        resseq=resseq,
        icode=icode,
        atomname="OG1",
        x=og1_xyz[0],
        y=og1_xyz[1],
        z=og1_xyz[2],
        occupancy=cb.occupancy,
        bfactor=cb.bfactor,
    )

    cg2 = AtomRecord(
        group="ATOM",
        atom_id=next_id + 1,
        element="C",
        resname="THR",
        chain_auth=cb.chain_auth,
        chain_pdb=cb.chain_pdb,
        segid=cb.segid,
        resseq=resseq,
        icode=icode,
        atomname="CG2",
        x=cg2_xyz[0],
        y=cg2_xyz[1],
        z=cg2_xyz[2],
        occupancy=cb.occupancy,
        bfactor=cb.bfactor,
    )

    mutated_records: List[AtomRecord] = []
    inserted = False
    for r in records:
        if r.chain_auth != chain_auth or r.resseq != resseq or r.icode != icode:
            mutated_records.append(r)
            continue

        atom = r.atomname.strip().upper()
        mutated_records.append(AtomRecord(**{**r.__dict__, "resname": "THR"}))

        if not inserted and atom == "CB":
            mutated_records.append(og1)
            mutated_records.append(cg2)
            inserted = True

    if not inserted:
        mutated_records.append(og1)
        mutated_records.append(cg2)

    return mutated_records


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cif",
        default=str(Path("output/playwright/chatgpt_botprompts/pdb/9TI4.cif")),
        help="Input mmCIF (default: output/playwright/chatgpt_botprompts/pdb/9TI4.cif)",
    )
    ap.add_argument(
        "--outdir",
        default=str(Path("output/playwright/chatgpt_botprompts/models")),
        help="Output directory (default: output/playwright/chatgpt_botprompts/models)",
    )
    args = ap.parse_args(argv)

    cif_path = Path(args.cif)
    outdir = Path(args.outdir)

    pdb_id = cif_path.stem.upper()
    source_url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    records, chain_map = parse_cif_atoms(cif_path, heavy_only=True, drop_hoh=True, model_num=1)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "complexI_9TI4_chain_map.json").write_text(
        json.dumps(chain_map, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    common_remarks = [
        "Template mmCIF downloaded from RCSB PDB:",
        source_url,
        f"Generated from {pdb_id} (heavy atoms only; HOH removed).",
    ]

    wt_full = outdir / "complexI_9TI4_WT_heavy.pdb"
    write_pdb_with_remarks(wt_full, records, remarks=common_remarks)

    wt_protein_only = outdir / "complexI_9TI4_WT_heavy_proteinOnly.pdb"
    write_pdb_with_remarks(
        wt_protein_only,
        [r for r in records if r.group == "ATOM"],
        remarks=common_remarks,
    )

    wt_nd6 = outdir / "nd6_chain_m_WT_heavy.pdb"
    write_pdb_with_remarks(
        wt_nd6,
        [r for r in records if r.chain_auth == "m"],
        remarks=common_remarks,
    )

    wt_nd1 = outdir / "nd1_chain_s_WT_heavy.pdb"
    write_pdb_with_remarks(
        wt_nd1,
        [r for r in records if r.chain_auth == "s"],
        remarks=common_remarks,
    )

    nd1_records = mutate_nd1_a52t(records, chain_auth="s", resseq=52, icode="")
    nd1_full = outdir / "complexI_9TI4_ND1_A52T_heavy.pdb"
    write_pdb_with_remarks(
        nd1_full,
        nd1_records,
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND1 chain s resid 52 ALA->THR.",
            "Variant: m.3460G>A (p.Ala52Thr).",
        ],
    )

    nd1_protein_only = outdir / "complexI_9TI4_ND1_A52T_heavy_proteinOnly.pdb"
    write_pdb_with_remarks(
        nd1_protein_only,
        [r for r in nd1_records if r.group == "ATOM"],
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND1 chain s resid 52 ALA->THR.",
            "Variant: m.3460G>A (p.Ala52Thr).",
        ],
    )

    nd1_chain = outdir / "nd1_chain_s_A52T_heavy.pdb"
    write_pdb_with_remarks(
        nd1_chain,
        [r for r in nd1_records if r.chain_auth == "s"],
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND1 chain s resid 52 ALA->THR.",
            "Variant: m.3460G>A (p.Ala52Thr).",
        ],
    )

    mut_records = mutate_nd6_m64v(records, chain_auth="m", resseq=64, icode="")
    mut_full = outdir / "complexI_9TI4_ND6_M64V_heavy.pdb"
    write_pdb_with_remarks(
        mut_full,
        mut_records,
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND6 chain m resid 64 MET->VAL.",
            "Variant: m.14484T>C (p.Met64Val).",
        ],
    )

    mut_protein_only = outdir / "complexI_9TI4_ND6_M64V_heavy_proteinOnly.pdb"
    write_pdb_with_remarks(
        mut_protein_only,
        [r for r in mut_records if r.group == "ATOM"],
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND6 chain m resid 64 MET->VAL.",
            "Variant: m.14484T>C (p.Met64Val).",
        ],
    )

    mut_nd6 = outdir / "nd6_chain_m_M64V_heavy.pdb"
    write_pdb_with_remarks(
        mut_nd6,
        [r for r in mut_records if r.chain_auth == "m"],
        remarks=common_remarks
        + [
            "Mutation applied: MT-ND6 chain m resid 64 MET->VAL.",
            "Variant: m.14484T>C (p.Met64Val).",
        ],
    )

    print(f"Wrote: {wt_full}")
    print(f"Wrote: {wt_protein_only}")
    print(f"Wrote: {wt_nd1}")
    print(f"Wrote: {nd1_full}")
    print(f"Wrote: {nd1_protein_only}")
    print(f"Wrote: {nd1_chain}")
    print(f"Wrote: {wt_nd6}")
    print(f"Wrote: {mut_full}")
    print(f"Wrote: {mut_protein_only}")
    print(f"Wrote: {mut_nd6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
