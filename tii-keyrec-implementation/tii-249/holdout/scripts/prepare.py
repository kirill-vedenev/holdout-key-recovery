#!/usr/bin/env python3
"""Prepare a single-position TII-249 holdout instance from public data only."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def rref(rows, n):
    a = list(rows)
    pivots = []
    for c in range(n):
        found = next((i for i in range(len(pivots), len(a)) if (a[i] >> c) & 1), None)
        if found is None:
            continue
        p = len(pivots)
        a[p], a[found] = a[found], a[p]
        for i in range(len(a)):
            if i != p and (a[i] >> c) & 1:
                a[i] ^= a[p]
        pivots.append(c)
    return a[: len(pivots)], pivots


def nullspace(rows, n):
    a, pivots = rref(rows, n)
    pivot_set = set(pivots)
    free = [j for j in range(n) if j not in pivot_set]
    basis = []
    for j in free:
        v = 1 << j
        for row, pivot in zip(a, pivots):
            v |= ((row >> j) & 1) << pivot
        basis.append(v)
    return basis, free


def restrict(v, keep):
    return sum(((v >> j) & 1) << i for i, j in enumerate(keep))


def read_public(path):
    lines = [s.strip() for s in Path(path).read_text().splitlines() if s.strip()]
    if len(lines) != 129:
        raise ValueError("TII-249 must contain 128 matrix rows and one field-modulus line")
    rows = []
    for s in lines[:-1]:
        if not s.startswith("[") or not s.endswith("]"):
            raise ValueError("invalid public matrix row")
        bits = s[1:-1].split()
        if len(bits) != 235 or any(b not in ("0", "1") for b in bits):
            raise ValueError("expected a binary 128-by-235 matrix")
        rows.append(sum(int(b) << j for j, b in enumerate(bits)))
    modulus = json.loads(lines[-1])
    if len(modulus) != 9 or modulus[-1] != 1 or any(type(b) is not int or b not in (0, 1) for b in modulus):
        raise ValueError("invalid GF(256) modulus")
    if len(rref(rows, 235)[1]) != 128:
        raise ValueError("public parity check does not have the published rank 128")
    return rows, modulus


def retained_orders(d, multiplicity):
    keep = []
    for u in range(multiplicity - 1, -1, -1):
        if not any(math.comb(d - u, v - u) & 1 for v in keep):
            keep.append(u)
    return sorted(keep)


def kappa(k, d, p):
    coeffs = [
        sum((-1) ** j * math.comb(p, j) * math.comb(k, h - 2 * j) for j in range(min(p, h // 2) + 1))
        for h in range(d + 1)
    ]
    return (coeffs[-1] if min(coeffs) > 0 else 0), coeffs


def make_saarinen_request(public, spec, output, block_bits=512):
    levels = spec["retained_derivative_orders"]
    if block_bits % 64:
        raise ValueError("Saarinen worker block width must be a multiple of 64")
    block_words = block_bits // 64
    lines = [
        "mceliecex-holdout-cuda-request-v1",
        f"public_sha256 {public['public_sha256']}",
        f"k {public['parameters']['k']}",
        f"degree {public['parameters']['degree']}",
        f"words {block_words}",
        "repetitions 1",
        "seed 24920260930",
        "emit_output 0",
        f"expected_columns {spec['monomials']}",
        f"expected_rows {spec['reduced_constraint_rows']}",
        f"expected_nonzeros {spec['exact_nonzero_entries']}",
        f"point_count {spec['active_constraint_points']}",
    ]
    for point in public["points"]:
        if point["mask"].bit_count() < public["parameters"]["degree"] - point["multiplicity"] + 1:
            continue
        lines.append("point {} {} {}".format(point["mask"], len(levels), " ".join(map(str, levels))))
    output.write_text("\n".join(lines) + "\n")


def prepare(output, heldout=0, k=50, default_elements=512):
    src = ROOT.parents[2] / "tii-results" / "tii-249" / "pk_McEliece_249.txt"
    provenance = json.loads((ROOT.parents[2] / "tii-results" / "tii-249" / "source.json").read_text())
    if sha256(src) != provenance["sha256"]:
        raise ValueError("official public input SHA-256 mismatch")
    H, modulus = read_public(src)
    if not 0 <= heldout < 235 or not 8 <= k <= 107:
        raise ValueError("require one original position in [0,234] and 8<=k<=107")

    _, full_free = nullspace(H, 235)
    shortened = [j for j in full_free if j != heldout][: 107 - k]
    if len(shortened) != 107 - k:
        raise ValueError("could not shorten to the requested dimension")
    keep = [j for j in range(235) if j not in shortened]
    hs = [restrict(row, keep) for row in H]
    G, free = nullspace(hs, len(keep))
    if len(G) != k:
        raise ValueError("shortening did not produce the requested dimension")

    held_local = keep.index(heldout)
    order = [j for j in free if j != held_local]
    order += [j for j in range(len(keep)) if j not in order and j != held_local]
    _, pivots_in_order = rref([restrict(row, order) for row in G], len(order))
    info = [order[j] for j in pivots_in_order]
    if len(info) != k:
        raise ValueError("no information set avoids the held-out position")

    augmented = [restrict(row, info) | (row << k) for row in G]
    aug, ap = rref(augmented, k + len(keep))
    if ap != list(range(k)):
        raise ValueError("information-set normalization failed")
    Y = [row >> k for row in aug]
    for i, row in enumerate(Y):
        if restrict(row, info) != 1 << i:
            raise ValueError("systematic generator invariant failed")
        if any((row & h).bit_count() & 1 for h in hs):
            raise ValueError("shortened generator is outside the public code")

    columns = [sum(((row >> j) & 1) << i for i, row in enumerate(Y)) for j in range(len(keep))]
    infinity = 0
    for p in columns:
        infinity ^= p
    if not columns[held_local] or columns[held_local].bit_count() == 1:
        raise ValueError("held-out column is degenerate in this chart; choose another public position explicitly")

    m, t, d, s = 8, 16, 5, 4
    minimum = max(10 * math.comb(m, 2), 5 * k)
    planned = max(default_elements, ((minimum + 63) // 64) * 64)
    points = [
        {"original": keep[j], "mask": columns[j], "multiplicity": s}
        for j in range(len(keep))
        if j != held_local
    ]
    points.append({"original": -1, "mask": infinity, "multiplicity": s})
    active = [p for p in points if p["mask"].bit_count() >= d - s + 1]
    levels = retained_orders(d, s)
    nr = len(active) * sum(math.comb(k, u) for u in levels)
    nmon = math.comb(k, d)
    nnz = 0
    for p in active:
        w = p["mask"].bit_count()
        for u in levels:
            for a in range(u + 1):
                if a <= w and u - a <= k - w and d - u <= w - a:
                    nnz += math.comb(w, a) * math.comb(k - w, u - a) * math.comb(w - a, d - u)

    n = len(keep)
    forced = s * n
    prediction, coeffs = kappa(k, d, n - k)
    public = dict(
        parameters=dict(m=m, n=n, t=t, k=k, degree=d),
        challenge="TII-249",
        original_parameters=dict(m=8, n=235, t=16, k=107),
        original_indices=keep,
        shortened_original_indices=shortened,
        information_set=info,
        information_set_original=[keep[j] for j in info],
        heldout=[held_local],
        heldout_original=heldout,
        heldout_mask=columns[held_local],
        multiplicities=[0 if j == held_local else s for j in range(n)],
        infinity_multiplicity=s,
        infinity_mask=infinity,
        field_modulus=modulus,
        generator=[[(row >> j) & 1 for j in range(n)] for row in Y],
        random_control=False,
        public_sha256=provenance["sha256"],
        minimum_elements=minimum,
        default_elements=planned,
        points=points,
    )
    identity = hashlib.sha256(canonical(public)).hexdigest()
    public["instance_id"] = identity
    spec = dict(
        instance_id=identity,
        single_position=True,
        heldout_original=heldout,
        minimum_elements=minimum,
        default_elements=planned,
        target_formula="max(10*binom(m,2),5*k_short)",
        k_short=k,
        original_k=107,
        monomial_order="colex",
        monomials=nmon,
        reduced_constraint_rows=nr,
        packed_matrix_bytes=nr * ((nmon + 63) // 64) * 8,
        exact_nonzero_entries=nnz,
        retained_derivative_orders=levels,
        all_derivative_orders=list(range(s)),
        active_constraint_points=len(active),
        implicit_zero_constraint_points=len(points) - len(active),
        full_constraint_rows=len(points) * sum(math.comb(k, u) for u in range(s)),
        kappa_prediction=prediction,
        kappa_coefficients=coeffs,
        forced_zeros=forced,
        composition_degree_bound=d * (n - 2 * t - 1),
        shortfall=d * (n - 2 * t - 1) + 1 - forced,
        ordinary_zero_count_suffices=forced > d * (n - 2 * t - 1),
        weis_zero_count_assumption_suffices=d * (n - 2 * t - 1) + 1 - forced <= t,
        kernel_membership_will_be_exact=True,
        curve_validity_is_not_certified_by_kernel_membership=True,
        sparse_linear_algebra="Saarinen holdout CUDA worker / CADO lingen path",
        dense_elimination_feasible=False,
    )
    if not spec["weis_zero_count_assumption_suffices"] or prediction < minimum:
        raise ValueError("chosen chart parameters fail the configured parameter criteria")

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    expected = {
        "public.json",
        "matrix-spec.json",
        "operator.txt",
        "monomials-colex.u64",
        "checksums.json",
        "saarinen-holdout-request.txt",
    }
    if any((output / name).exists() for name in expected):
        raise FileExistsError("prepared files already exist; choose a new output directory")
    (output / "public.json").write_text(json.dumps(public, indent=2) + "\n")
    (output / "matrix-spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    text = (
        f"TII_HOLDOUT_V1\ninstance_id {identity}\nk {k}\ndegree {d}\n"
        f"minimum {minimum}\nheldout_original {heldout}\nheldout_mask {columns[held_local]}\npoints {len(points)}\n"
    )
    text += "".join(f"{p['original']} {p['mask']} {p['multiplicity']}\n" for p in points)
    (output / "operator.txt").write_text(text)
    mons = sorted(sum(1 << j for j in a) for a in itertools.combinations(range(k), d))
    with (output / "monomials-colex.u64").open("wb") as f:
        for mon in mons:
            f.write(struct.pack("<Q", mon))
    make_saarinen_request(public, spec, output / "saarinen-holdout-request.txt")
    checksums = {name: sha256(output / name) for name in sorted(expected - {"checksums.json"})}
    (output / "checksums.json").write_text(json.dumps(checksums, indent=2) + "\n")
    return spec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=ROOT / "prepared")
    ap.add_argument("--heldout-original", type=int, default=0)
    ap.add_argument("--shortened-dimension", type=int, default=50)
    ap.add_argument("--default-elements", type=int, default=512)
    args = ap.parse_args()
    try:
        print(json.dumps(prepare(args.output, args.heldout_original, args.shortened_dimension, args.default_elements), indent=2))
    except (ValueError, OSError) as exc:
        ap.exit(1, f"Preparation failed: {exc}\n")


if __name__ == "__main__":
    main()
