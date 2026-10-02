"""Summarizes the JSON-lines results of product_rank into markdown tables.

Usage: python3 summarize.py RESULTS_DIR

Reads mceliece.jsonl, edge.jsonl, toy.jsonl, moderate.jsonl and full.jsonl
(or their .gz forms) and crosscheck.json when present, and writes summary.md
and summary.json to the same directory. Uses only the Python standard library.
"""

import collections
import gzip
import json
import math
import os
import statistics
import sys


def load(path):
    """Reads PATH or its gzip-compressed form PATH.gz."""
    if os.path.exists(path):
        f = open(path)
    elif os.path.exists(path + ".gz"):
        f = gzip.open(path + ".gz", "rt")
    else:
        return []
    with f:
        return [json.loads(line) for line in f if line.strip()]


def fmt_time(xs):
    return "%.2f" % statistics.median(xs) if xs else "-"


def classify(r):
    """Compares the deficit of a key above the counting threshold with the
    number predicted by special points of the curve, including infinity:
    one per extra copy of a point (equal columns of the parity-extended
    generator matrix), two per base point (zero column), and one per
    degenerate tangent."""
    if r.get("predicted_deficit", -1) < 0:
        return "unknown"
    return "explained" if r["predicted_deficit"] == r["deficit"] else "unexplained"


def global_checks(rows):
    bad = collections.Counter()
    for r in rows:
        c = r["checks"]
        bad["inconsistent"] += not c["consistent"]
        bad["goppa_square_failed"] += not c["goppa_square"]
        bad["f_degree_failed"] += not c["f_degree"]
        bad["product_not_in_S"] += not c["product_in_S"]
        bad["rank_above_bound"] += r["rank"] > r["bound"]
        bad["control_in_W"] += r["control_in_W"] > 0
        if r["reached_bound"]:
            bad["missing_R_or_S_at_bound"] += r["R_and_S"] != r["positions_tested"]
    return {k: bad[k] for k in sorted(bad)}


def mceliece_table(rows):
    groups = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["scheme"], r["full_n"], r["label"])):
        groups.setdefault(r["label"], []).append(r)
    lines = [
        "| Parameter set | Shortening | (m, t) | n | k | D | 2D-t+1 | C(k+1,2) | Keys | dim W = 2D-t+1 | Positions with R and S | Median s |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    out = []
    for label, rs in groups.items():
        r0 = rs[0]
        ks = sorted({r["k"] for r in rs})
        kstr = str(ks[0]) if len(ks) == 1 else "%d-%d" % (ks[0], ks[-1])
        reached = sum(r["reached_bound"] for r in rs)
        pos = sum(r["positions_tested"] for r in rs)
        both = sum(r["R_and_S"] for r in rs)
        name = label.rsplit("-", 1)[0]
        scheme = "Weis" if "weis" in label else "GIAJS"
        lines.append(
            "| `%s` | %s | (%d, %d) | %d | %s | %d | %d | %d | %d | %d | %d / %d | %s |"
            % (name, scheme, r0["m"], r0["t"], r0["n"], kstr, r0["D"], r0["bound"],
               r0["products"], len(rs), reached, both, pos,
               fmt_time([r["seconds"]["total"] for r in rs])))
        out.append({"label": label, "keys": len(rs), "reached": reached,
                    "positions": pos, "R_and_S": both, "k_values": ks})
    return lines, out


def edge_table(rows):
    groups = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["m"], r["t"], r["k"])):
        groups.setdefault((r["m"], r["t"], r["k"]), []).append(r)
    lines = [
        "| (m, t) | k | n | 2D-t+1 | C(k+1,2) | Keys | dim W = min(C(k+1,2), 2D-t+1) | Observed dim W |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    out = []
    for (m, t, k), rs in groups.items():
        r0 = rs[0]
        expected = min(r0["products"], r0["bound"])
        good = sum(r["rank"] == expected for r in rs)
        dims = sorted({r["rank"] for r in rs})
        lines.append("| (%d, %d) | %d | %d | %d | %d | %d | %d | %s |"
                     % (m, t, k, r0["n"], r0["bound"], r0["products"], len(rs), good,
                        ", ".join(map(str, dims))))
        out.append({"m": m, "t": t, "k": k, "keys": len(rs), "generic": good})
    return lines, out


def toy_table(rows):
    groups = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["m"], r["t"], r["k"])):
        groups.setdefault((r["m"], r["t"], r["n"] - r["m"] * r["t"]), []).append(r)
    lines = [
        "| (m, t) | k | n | 2D-t+1 | C(k+1,2) | Slack | Keys | Deficient | Explained | Unexplained | Estimate | Deficient keys with R and S somewhere |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    out = []
    for (m, t, kl), rs in groups.items():
        r0 = rs[0]
        above = r0["products"] >= r0["bound"]
        n = r0["n"]
        # Probability that the n + 1 columns of the parity-extended generator
        # matrix (uniform in F_2^k) contain a repeated or a zero column.
        birthday = -math.expm1(-((n + 1) * n / 2 + (n + 1)) / 2.0 ** kl)
        if not above:
            indep = sum(r["rank"] == r["products"] for r in rs)
            lines.append("| (%d, %d) | %d | %d | %d | %d | %d | %d | below threshold: %d keys with independent products | | | | |"
                         % (m, t, kl, n, r0["bound"], r0["products"], r0["products"] - r0["bound"], len(rs), indep))
            out.append({"m": m, "t": t, "k": kl, "keys": len(rs), "below_threshold": True,
                        "independent_products": indep})
            continue
        defic = [r for r in rs if r["deficit"] > 0]
        cls = collections.Counter(classify(r) for r in defic)
        some = sum(r["R_and_S"] > 0 for r in defic)
        lines.append("| (%d, %d) | %d | %d | %d | %d | %d | %d | %d (%.1f%%) | %d | %d%s | %.1f%% | %d |"
                     % (m, t, kl, n, r0["bound"], r0["products"], r0["products"] - r0["bound"], len(rs), len(defic),
                        100.0 * len(defic) / len(rs), cls["explained"], cls["unexplained"],
                        " (+%d unknown)" % cls["unknown"] if cls["unknown"] else "",
                        100 * birthday, some))
        out.append({"m": m, "t": t, "k": kl, "n": n, "keys": len(rs), "deficient": len(defic),
                    "explained": cls["explained"], "unexplained": cls["unexplained"],
                    "unknown": cls["unknown"], "birthday_estimate": birthday,
                    "deficient_with_R_and_S": some})
    return lines, out


def moderate_table(rows):
    groups = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["m"], r["t"], r["k"])):
        groups.setdefault((r["m"], r["t"], r["k"]), []).append(r)
    lines = [
        "| (m, t) | k | n | 2D-t+1 | C(k+1,2) | Slack | Keys | dim W = 2D-t+1 | Deficient | Explained | Positions with R and S |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    out = []
    for (m, t, k), rs in groups.items():
        r0 = rs[0]
        defic = [r for r in rs if r["deficit"] > 0]
        cls = collections.Counter(classify(r) for r in defic)
        reached = sum(r["reached_bound"] for r in rs)
        pos = sum(r["positions_tested"] for r in rs if r["reached_bound"])
        both = sum(r["R_and_S"] for r in rs if r["reached_bound"])
        lines.append("| (%d, %d) | %d | %d | %d | %d | %d | %d | %d | %d | %d | %d / %d |"
                     % (m, t, k, r0["n"], r0["bound"], r0["products"],
                        r0["products"] - r0["bound"], len(rs), reached, len(defic),
                        cls["explained"], both, pos))
        out.append({"m": m, "t": t, "k": k, "keys": len(rs), "reached": reached,
                    "deficient": len(defic), "explained": cls["explained"],
                    "unexplained": cls["unexplained"], "unknown": cls["unknown"]})
    return lines, out


def full_table(rows):
    groups = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["m"], r["n"])):
        groups.setdefault(r["label"], []).append(r)
    lines = [
        "| Parameters | (m, t) | n | k | D | 2D-t+1 | C(k+1,2) | Keys | dim W = 2D-t+1 | Deficient | Explained | Positions with R and S |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    out = []
    for label, rs in groups.items():
        r0 = rs[0]
        ks = sorted({r["k"] for r in rs})
        kstr = str(ks[0]) if len(ks) == 1 else "%d-%d" % (ks[0], ks[-1])
        defic = [r for r in rs if r["deficit"] > 0]
        cls = collections.Counter(classify(r) for r in defic)
        reached = sum(r["reached_bound"] for r in rs)
        pos = sum(r["positions_tested"] for r in rs if r["reached_bound"])
        both = sum(r["R_and_S"] for r in rs if r["reached_bound"])
        lines.append("| %s | (%d, %d) | %d | %s | %d | %d | %d | %d | %d | %d | %d | %d / %d |"
                     % (label.replace("full-", "").upper(), r0["m"], r0["t"], r0["n"], kstr, r0["D"],
                        r0["bound"], r0["products"], len(rs), reached, len(defic),
                        cls["explained"], both, pos))
        out.append({"label": label, "keys": len(rs), "reached": reached, "deficient": len(defic),
                    "explained": cls["explained"], "unexplained": cls["unexplained"]})
    return lines, out


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "results"
    mce = load(os.path.join(d, "mceliece.jsonl"))
    edge = load(os.path.join(d, "edge.jsonl"))
    toy = load(os.path.join(d, "toy.jsonl"))
    moderate = load(os.path.join(d, "moderate.jsonl"))
    full = load(os.path.join(d, "full.jsonl"))
    everything = mce + edge + toy + moderate + full
    summary = {"checks": global_checks(everything),
               "instances": {"mceliece": len(mce), "edge": len(edge), "toy": len(toy),
                             "moderate": len(moderate), "full": len(full)}}
    md = ["# Summary of results", ""]
    md += ["Generated by `summarize.py` from the JSON-lines files in this directory.", ""]
    md += ["## Global checks", "",
           "Counts of keys failing a check (all should be zero) over %d keys." % len(everything), ""]
    md += ["| Check | Failures |", "|---|---:|"]
    md += ["| %s | %d |" % (k, v) for k, v in summary["checks"].items()]
    if mce:
        lines, summary["mceliece"] = mceliece_table(mce)
        md += ["", "## Classic McEliece keys", ""] + lines
    if edge:
        lines, summary["edge"] = edge_table(edge)
        md += ["", "## Counting threshold at Classic McEliece (m, t)", ""] + lines
    if moderate:
        lines, summary["moderate"] = moderate_table(moderate)
        md += ["", "## Moderate shortened codes", "",
               "Positions are counted on keys that reach the bound."] + [""] + lines
    if full:
        lines, summary["full"] = full_table(full)
        md += ["", "## Unshortened codes with the TII challenge parameters", "",
               "The products are computed in coefficient form, since 2D+1 exceeds the field size."] + [""] + lines
    if toy:
        lines, summary["toy"] = toy_table(toy)
        md += ["", "## Small codes, all products used", ""] + lines
    cc_path = os.path.join(d, "crosscheck.json")
    if os.path.exists(cc_path):
        with open(cc_path) as f:
            cc = json.load(f)
        summary["crosscheck"] = {k: cc[k] for k in ("instances", "deficient", "all_agree")}
        md += ["", "## SageMath cross-check", "",
               "%d keys (%d deficient); agreement on rank, memberships and the cross-product rank: %s."
               % (cc["instances"], cc["deficient"], "yes" if cc["all_agree"] else "NO")]
    with open(os.path.join(d, "summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(d, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
