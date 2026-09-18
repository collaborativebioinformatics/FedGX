#!/usr/bin/env python3
"""
build_fedx_dashboard.py - turn real FedGX pipeline output into fedx_data.js
for the FedGX Explorer dashboard (fedx_dashboard.html).

Inputs (from your existing pipeline):
  - One `<COHORT>.GWAMA.txt` or `<COHORT>.cleaned.txt` per site, produced by
    gwas_cohort_qc_with_gwama.R. Binary OR/confidence-interval input and
    quantitative BETA/SE input are both normalised to the log-effect scale.
  - One combined meta-analysis output file from GWAMA (whitespace-delimited
    by default; column names default to GWAMA's usual header - rs_number,
    beta, se, n_samples - overridable via --meta-col-* since GWAMA's exact
    header can vary by version/invocation).

Sites are assigned generic labels in the order given on the command line:
the first --site is "site1", the second "site2", and so on - no cohort
names are baked into the dashboard.

Because gwas_cohort_qc_with_gwama.R builds its VARIANT_ID (and therefore
GWAMA's MARKERNAME) as "CHR:POS:allele:allele", this script recovers
chromosome and position directly from the marker string - no need to
cross-reference back to the site files for genomic coordinates.

Example:
  python3 build_fedx_dashboard.py \\
    --site COHORT1.cleaned.txt \\
    --site COHORT2.cleaned.txt \\
    --site COHORT3.cleaned.txt \\
    --meta META.out.txt \\
    --host "NVIDIA Brev" \\
    --out-dir .

Writes: <out-dir>/fedx_data.js
Then just open fedx_dashboard.html (from the same build as this script)
with fedx_data.js sitting next to it in the same folder - no server needed.
"""

import argparse
import json
import math
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

DEFAULT_SITE_COLORS = ["#0064A6", "#D54E12", "#86A10B", "#8AC2E6", "#F0B600", "#66757E"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Build fedx_data.js for FedGX Explorer from real site + GWAMA output files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--site", dest="sites", action="append", required=True, metavar="PATH",
                    help="Path to one site's .GWAMA.txt or .cleaned.txt file. Repeat in site order.")
    p.add_argument("--meta", required=True, metavar="PATH",
                    help="Path to the combined GWAMA meta-analysis output file.")
    p.add_argument("--model", choices=["fixed", "random"], default="fixed",
                    help="Model represented by --meta [default: %(default)s]")
    p.add_argument("--trait-type", choices=["auto", "binary", "quantitative"],
                    default="auto", help="Input trait type [default: %(default)s]")
    p.add_argument("--method", default="regenie",
                    help="Site GWAS method recorded in dashboard metadata")

    p.add_argument("--host", default="NVIDIA Brev",
                    help="Host label applied to every site unless overridden with --site-host [default: %(default)s]")
    p.add_argument("--site-host", dest="site_hosts", action="append", default=[], metavar="TEXT",
                    help="Per-site host override, same order as --site. May be passed fewer times than --site.")

    # site (.cleaned.txt) column mapping - defaults match gwas_cohort_qc_with_gwama.R's own output
    p.add_argument("--site-delim", default="\t", help="Site file delimiter [default: tab]")
    p.add_argument("--site-col-marker", default="VARIANT_ID")
    p.add_argument("--site-col-chr", default="CHR")
    p.add_argument("--site-col-pos", default="POS")
    p.add_argument("--site-col-rsid", default="SNP")
    p.add_argument("--site-col-eaf", default="EAF")
    p.add_argument("--site-col-beta", default="BETA")
    p.add_argument("--site-col-se", default="SE")
    p.add_argument("--site-col-n", default="N")
    p.add_argument("--site-col-z", default="Z", help="Set to NULL to derive Z = BETA/SE instead of reading a column")

    # meta (GWAMA output) column mapping - defaults match GWAMA's typical header;
    # override these if your GWAMA invocation produced different column names
    p.add_argument("--meta-delim", default=None,
                    help="Meta file delimiter; default auto-detects whitespace-or-tab")
    p.add_argument("--meta-col-marker", default="rs_number")
    p.add_argument("--meta-col-beta", default="beta")
    p.add_argument("--meta-col-se", default="se")
    p.add_argument("--meta-col-n", default="n_samples")
    p.add_argument("--meta-col-z", default="NULL", help="Set to a column name to read Z directly instead of deriving BETA/SE")

    p.add_argument("--out-dir", default=".", help="Output directory [default: current directory]")
    p.add_argument("--out-file", default="fedx_data.js", help="Output filename [default: %(default)s]")
    p.add_argument("--max-variants", type=int, default=50000,
                    help="Maximum variants embedded in the browser payload; 0 keeps all [default: %(default)s]")
    p.add_argument("--html-template", default=None, metavar="PATH",
                    help="Copy this dashboard HTML template to <out-dir>/index.html")
    return p.parse_args()


def na_to_null(x):
    if x is None or str(x).strip().upper() == "NULL" or str(x).strip() == "":
        return None
    return x


def require_columns(df, cols, label):
    missing = [c for c in cols if c and c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {label} is missing required column(s) {missing}. "
                  f"Available columns: {list(df.columns)}")


def first_present(columns, candidates):
    names = set(columns)
    return next((candidate for candidate in candidates if candidate in names), None)


def numeric(df, column):
    return pd.to_numeric(df[column], errors="coerce")


def log_positive(values):
    return values.map(lambda value: math.log(value) if pd.notna(value) and value > 0 else math.nan)


def effect_and_se(df, label, beta_hint=None, se_hint=None):
    """Return log-scale effect and SE from quantitative or binary columns."""
    beta_col = first_present(
        df.columns,
        [beta_hint, "beta", "BETA"] if beta_hint else ["beta", "BETA"],
    )
    se_col = first_present(
        df.columns,
        [se_hint, "beta_se", "BETA_se", "SE", "se"] if se_hint
        else ["beta_se", "BETA_se", "SE", "se"],
    )
    if beta_col and se_col:
        return numeric(df, beta_col), numeric(df, se_col), "quantitative"

    or_col = first_present(df.columns, ["OR", "or"])
    low_col = first_present(df.columns, ["OR_95L", "or_95l"])
    high_col = first_present(df.columns, ["OR_95U", "or_95u"])
    if or_col and low_col and high_col:
        odds = numeric(df, or_col)
        lower = log_positive(numeric(df, low_col))
        upper = log_positive(numeric(df, high_col))
        return log_positive(odds), (upper - lower) / (2 * 1.96), "binary"

    raise ValueError(
        f"{label} has neither BETA/SE nor OR/OR_95L/OR_95U columns; "
        f"available columns: {list(df.columns)}"
    )


def parse_marker_chr_pos(marker: str):
    # VARIANT_ID / MARKERNAME convention from gwas_cohort_qc_with_gwama.R:
    # "CHR:POS:allele:allele" (alleles alphabetically sorted, not ref/alt).
    parts = str(marker).split(":")
    if len(parts) < 2:
        return None, None
    chr_raw, pos_raw = parts[0], parts[1]
    try:
        pos = int(float(pos_raw))
    except ValueError:
        return None, None
    return chr_raw, pos


def load_site(path, args, site_index):
    delim = args.site_delim
    df = pd.read_csv(path, sep=delim, dtype=str, low_memory=False)
    label = f"site file {path}"
    marker_col = first_present(
        df.columns,
        [args.site_col_marker, "VARIANT_ID", "MARKERNAME", "rs_number"],
    )
    eaf_col = first_present(df.columns, [args.site_col_eaf, "EAF", "eaf"])
    n_col = first_present(df.columns, [args.site_col_n, "N", "n_samples"])
    if not marker_col or not eaf_col or not n_col:
        raise ValueError(
            f"{label} needs marker, EAF, and N columns; available: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["marker"] = df[marker_col].astype(str)
    chr_col = first_present(df.columns, [args.site_col_chr, "CHR", "chr"])
    pos_col = first_present(df.columns, [args.site_col_pos, "POS", "pos"])
    if chr_col and pos_col:
        out["chr"] = df[chr_col].astype(str)
        out["pos"] = numeric(df, pos_col)
    else:
        parsed = out["marker"].map(parse_marker_chr_pos)
        out["chr"] = parsed.map(lambda item: item[0])
        out["pos"] = parsed.map(lambda item: item[1])

    rsid_col = first_present(df.columns, [args.site_col_rsid, "SNP", "rsid"])
    out["rsid"] = df[rsid_col].astype(str) if rsid_col else out["marker"]
    out["eaf"] = numeric(df, eaf_col)
    out["n"] = numeric(df, n_col)
    out["beta"], out["se"], detected_trait = effect_and_se(
        df, label, args.site_col_beta, args.site_col_se
    )

    if args.trait_type != "auto" and detected_trait != args.trait_type:
        raise ValueError(
            f"{label} looks {detected_trait}, but --trait-type is {args.trait_type}"
        )

    col_z = first_present(df.columns, [na_to_null(args.site_col_z), "Z", "z"])
    if col_z:
        out["z"] = numeric(df, col_z)
    else:
        out["z"] = out["beta"] / out["se"]

    before = len(out)
    out = out.dropna(subset=["marker", "eaf", "beta", "se", "n", "z"])
    out = out[out["se"] > 0]
    dropped = before - len(out)
    if dropped:
        print(f"  site{site_index}: dropped {dropped} row(s) with missing/invalid values")

    out = out.drop_duplicates(subset="marker", keep="first")
    return out.set_index("marker")


def load_meta(path, args):
    delim = args.meta_delim
    if delim is None:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")
    else:
        df = pd.read_csv(path, sep=delim, dtype=str)

    label = f"meta file {path}"
    marker_col = first_present(
        df.columns,
        [args.meta_col_marker, "rs_number", "MARKERNAME", "VARIANT_ID"],
    )
    if not marker_col:
        raise ValueError(f"{label} has no recognizable marker column")

    out = pd.DataFrame()
    out["marker"] = df[marker_col].astype(str)
    out["beta"], out["se"], detected_trait = effect_and_se(
        df, label, args.meta_col_beta, args.meta_col_se
    )
    if args.trait_type != "auto" and detected_trait != args.trait_type:
        raise ValueError(
            f"{label} looks {detected_trait}, but --trait-type is {args.trait_type}"
        )

    n_col = first_present(df.columns, [args.meta_col_n, "n_samples", "N"])
    if n_col:
        out["n"] = numeric(df, n_col)
        out.loc[out["n"] <= 0, "n"] = math.nan
    else:
        out["n"] = None

    col_z = first_present(df.columns, [na_to_null(args.meta_col_z), "z", "Z"])
    if col_z:
        out["z"] = numeric(df, col_z)
    else:
        out["z"] = out["beta"] / out["se"]

    p_col = first_present(df.columns, ["p-value", "P", "p", "p_value"])
    out["p"] = numeric(df, p_col) if p_col else math.nan
    q_col = first_present(df.columns, ["q_statistic", "Q"])
    qp_col = first_present(df.columns, ["q_p-value", "Q_P"])
    i2_col = first_present(df.columns, ["i2", "I2"])
    out["q"] = numeric(df, q_col) if q_col else math.nan
    out["q_p"] = numeric(df, qp_col) if qp_col else math.nan
    out["i2"] = numeric(df, i2_col) if i2_col else math.nan

    before = len(out)
    out = out.dropna(subset=["marker", "beta", "se", "z"])
    out = out[out["se"] > 0]
    dropped = before - len(out)
    if dropped:
        print(f"  meta: dropped {dropped} row(s) with missing/invalid values")

    out = out.drop_duplicates(subset="marker", keep="first")
    return out.set_index("marker")


def build_chr_info(all_chr_pos):
    """Chromosome layout from the data actually present (not an assumed
    full-genome assembly) - each chromosome's plotted span runs from its
    min to max observed position, in the order chromosomes first appear
    sorted numerically (falling back to string order for non-numeric
    contigs, e.g. X/Y)."""
    spans = {}
    for chr_raw, pos in all_chr_pos:
        lo, hi = spans.get(chr_raw, (pos, pos))
        spans[chr_raw] = (min(lo, pos), max(hi, pos))

    def sort_key(c):
        try:
            return (0, int(c))
        except ValueError:
            return (1, str(c))

    ordered = sorted(spans.keys(), key=sort_key)
    lengths = {c: max(1, spans[c][1] - spans[c][0]) for c in ordered}
    total = sum(lengths.values())

    chr_info = []
    chr_offset = {}
    chr_min = {}
    cum = 0
    for c in ordered:
        start_frac = cum / total
        cum += lengths[c]
        end_frac = cum / total
        chr_info.append({"chr": c, "startFrac": start_frac, "endFrac": end_frac})
        chr_offset[c] = (start_frac, end_frac)
        chr_min[c] = spans[c][0]
    return chr_info, chr_offset, chr_min, lengths


def main():
    args = parse_args()
    if args.max_variants < 0:
        sys.exit("ERROR: --max-variants must be zero or positive")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_sites = len(args.sites)
    site_keys = [f"site{i+1}" for i in range(n_sites)]
    site_labels = [f"Site {i+1}" for i in range(n_sites)]
    site_hosts = list(args.site_hosts) + [args.host] * (n_sites - len(args.site_hosts))
    site_colors = [DEFAULT_SITE_COLORS[i % len(DEFAULT_SITE_COLORS)] for i in range(n_sites)]

    print("========================================")
    print("FedGX dashboard data build")
    print("========================================")
    print(f"Sites  : {n_sites}")
    for i, path in enumerate(args.sites):
        print(f"  {site_keys[i]:8s} <- {path}  (host: {site_hosts[i]})")
    print(f"Meta   : {args.meta}")
    print("========================================")

    print("\nReading site files...")
    site_frames = {}
    for i, path in enumerate(args.sites):
        key = site_keys[i]
        df = load_site(path, args, i + 1)
        site_frames[key] = df
        print(f"  {key}: {len(df):,} variants")

    print("\nReading meta file...")
    meta_df = load_meta(args.meta, args)
    print(f"  meta: {len(meta_df):,} variants")

    all_markers = set(meta_df.index)
    for df in site_frames.values():
        all_markers |= set(df.index)
    all_markers = sorted(all_markers)
    print(f"\nTotal distinct markers across meta + all sites: {len(all_markers):,}")

    # chromosome/position layout, parsed from the marker string itself
    chr_pos_by_marker = {}
    unparsed = 0
    for m in all_markers:
        c, pos = parse_marker_chr_pos(m)
        if c is None:
            unparsed += 1
            continue
        chr_pos_by_marker[m] = (c, pos)
    if unparsed:
        print(f"WARNING: {unparsed} marker(s) could not be parsed as CHR:POS:... and will be dropped "
              f"(expected GWAMA MARKERNAME format from gwas_cohort_qc_with_gwama.R, e.g. '7:130520000:C:T')")

    markers = [m for m in all_markers if m in chr_pos_by_marker and m in meta_df.index]
    skipped_no_meta = len(all_markers) - unparsed - len(markers)
    if skipped_no_meta > 0:
        print(f"NOTE: {skipped_no_meta} marker(s) present in a site but absent from the meta file were skipped "
              f"(the dashboard shows the meta-analysis view, so a variant needs a combined result to appear)")

    total_meta_markers = len(markers)
    if args.max_variants and len(markers) > args.max_variants:
        # Preserve the strongest 20% of displayed associations, then sample
        # the remaining capacity evenly across genomic order. This keeps lead
        # signals while bounding JavaScript/SVG size for full GWAS results.
        top_count = max(1, args.max_variants // 5)
        strongest = sorted(
            markers,
            key=lambda marker: abs(float(meta_df.loc[marker, "z"])),
            reverse=True,
        )[:top_count]
        strongest_set = set(strongest)
        genome_order = sorted(
            (marker for marker in markers if marker not in strongest_set),
            key=lambda marker: (
                int(chr_pos_by_marker[marker][0])
                if str(chr_pos_by_marker[marker][0]).isdigit() else 10_000,
                str(chr_pos_by_marker[marker][0]),
                chr_pos_by_marker[marker][1],
            ),
        )
        needed = args.max_variants - len(strongest)
        evenly_spaced = [
            genome_order[min((index * len(genome_order)) // needed,
                             len(genome_order) - 1)]
            for index in range(needed)
        ] if needed and genome_order else []
        selected = strongest_set | set(evenly_spaced)
        markers = [marker for marker in markers if marker in selected]
        print(
            f"Dashboard payload limited to {len(markers):,} of "
            f"{total_meta_markers:,} meta-analysis variants"
        )

    chr_info, chr_offset, chr_min, chr_len = build_chr_info(chr_pos_by_marker[m] for m in markers)
    total_len = sum(chr_len.values())

    print("\nAssembling per-SNP records...")
    coverage_counts = {k: 0 for k in site_keys}
    snps = []
    for m in markers:
        chr_raw, pos = chr_pos_by_marker[m]
        start_frac, _ = chr_offset[chr_raw]
        cum_frac = start_frac + (pos - chr_min[chr_raw]) / total_len

        site_obj = {}
        rsid = m
        for key in site_keys:
            df = site_frames[key]
            if m in df.index:
                row = df.loc[m]
                site_obj[key] = {
                    "eaf": round(float(row["eaf"]), 6),
                    "beta": round(float(row["beta"]), 6),
                    "se": round(float(row["se"]), 6),
                    "n": int(row["n"]),
                    "z": round(float(row["z"]), 6),
                }
                coverage_counts[key] += 1
                if row["rsid"] and row["rsid"] != "nan":
                    rsid = row["rsid"]
            else:
                site_obj[key] = None

        mrow = meta_df.loc[m]
        meta_obj = {
            "beta": round(float(mrow["beta"]), 6),
            "se": round(float(mrow["se"]), 6),
            "z": round(float(mrow["z"]), 6),
        }
        for source, target in (("p", "p"), ("q", "q"),
                               ("q_p", "qP"), ("i2", "i2")):
            if pd.notna(mrow[source]):
                meta_obj[target] = float(mrow[source])
        if pd.notna(mrow["n"]):
            meta_obj["n"] = int(mrow["n"])

        snps.append({
            "marker": m,
            "rsid": rsid,
            "chr": chr_raw,
            "posMb": pos / 1e6,
            "cumFrac": cum_frac,
            "site": site_obj,
            "meta": meta_obj,
        })

    snps.sort(key=lambda r: r["cumFrac"])

    sites_meta = []
    for i, key in enumerate(site_keys):
        n_for_site = site_frames[key]["n"].median() if len(site_frames[key]) else 0
        sites_meta.append({
            "key": key,
            "label": site_labels[i],
            "host": site_hosts[i],
            "n": int(n_for_site) if not math.isnan(n_for_site) else 0,
            "color": site_colors[i],
        })

    total_n = sum(s["n"] for s in sites_meta)
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    source_info = (
        f"Loaded {generated_at} from {n_sites} site summary-statistics file(s) "
        f"({', '.join(Path(p).name for p in args.sites)}) and GWAMA meta-analysis output "
        f"({Path(args.meta).name})."
    )

    model_label = (
        "fixed-effect inverse-variance"
        if args.model == "fixed"
        else "random-effects sensitivity analysis"
    )
    data = {
        "generatedAt": generated_at,
        "sourceInfo": source_info,
        "analysis": {
            "method": args.method.upper(),
            "traitType": args.trait_type,
            "model": args.model,
        },
        "sampling": {
            "totalVariants": total_meta_markers,
            "displayedVariants": len(markers),
            "sampled": len(markers) < total_meta_markers,
        },
        "meta": {
            "label": f"Combined {args.model}-effects meta-analysis",
            "sub": f"GWAMA, {model_label}",
            "n": total_n,
        },
        "sites": sites_meta,
        "chrInfo": chr_info,
        "snps": snps,
    }

    out_path = out_dir / args.out_file
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("// Generated by build_fedx_dashboard.py - do not edit by hand.\n")
        f.write("window.FEDX_DATA = ")
        json.dump(data, f, allow_nan=False)
        f.write(";\n")

    html_output = None
    if args.html_template:
        template = Path(args.html_template)
        if not template.is_file():
            sys.exit(f"ERROR: dashboard HTML template not found: {template}")
        html_output = out_dir / "index.html"
        shutil.copyfile(template, html_output)

    print("\n========================================")
    print("Done")
    print("========================================")
    print(f"Variants written : {len(snps):,}")
    for key in site_keys:
        pct = 100 * coverage_counts[key] / len(snps) if snps else 0
        print(f"  {key} coverage   : {coverage_counts[key]:,} / {len(snps):,} ({pct:.1f}%)")
    print(f"Output           : {out_path}")
    if html_output:
        print(f"Dashboard        : {html_output}")
    print("\nOpen fedx_dashboard.html from the same folder as this script (with "
          f"{args.out_file} sitting next to it) to view it - no server needed.")


if __name__ == "__main__":
    main()
