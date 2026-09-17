#!/usr/bin/env python3
"""Multi-ancestry genotypes for FedGen sites, to pair with LDAK --make-phenos.

  reference  shared SNPs, per-ancestry allele frequencies (Balding-Nichols tree), causal SNPs
             shared by all ancestries (effects correlated --rg) or unique to one ancestry
             -> causals_<ANC>.txt, effects_<ANC>.txt, causal_types.tsv
  site       one site's .bed/.bim/.fam/.covar with its ancestry's frequencies and LD
  gwas       covariate-adjusted linear GWAS of an LDAK .pheno (optional)
  meta       fixed-effect meta-analysis of site GWAS files with I^2
"""
import argparse
import math
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import NormalDist

import numpy as np

# name: (parent, Fst from parent, LD autocorrelation, mean LD-block length in SNPs); parents first
TREE = {"AFR": ("ANC", 0.07, 0.80, 25), "OOA": ("ANC", 0.14, 0, 1),
        "EUR": ("OOA", 0.08, 0.92, 60), "EAS": ("OOA", 0.12, 0.93, 70)}
ANCESTRIES = ["AFR", "EUR", "EAS"]
MAGIC = bytes([0x6C, 0x1B, 0x01])
ENCODE = np.array([3, 2, 0], np.uint8)       # A1 count 0/1/2 -> PLINK 2-bit code
DECODE = np.array([2, 0, 1, 0], np.float64)  # PLINK code -> A1 count (missing not expected)
SHIFTS = np.array([0, 2, 4, 6], np.uint8)
norm_ppf = np.vectorize(NormalDist().inv_cdf, otypes=[float])
_erfc = np.frompyfunc(math.erfc, 1, 1)


def pvalue(z):
    return _erfc(np.abs(z) / math.sqrt(2)).astype(float)


def reference(a):
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    m = a.num_snps
    chrom = np.sort(rng.integers(1, 23, m))
    pos = np.concatenate([np.cumsum(rng.integers(1000, 20000, (chrom == c).sum())) for c in range(1, 23)])
    alleles = np.array(list("ACGT"))[np.argsort(rng.random((m, 4)), 1)[:, :2]]
    maf = rng.uniform(0.01, 0.5, m)
    freq = {"ANC": np.where(rng.random(m) < 0.5, maf, 1 - maf)}
    for name, (parent, fst, _, _) in TREE.items():
        p, k = np.clip(freq[parent], 1e-9, 1 - 1e-9), (1 - fst) / fst
        freq[name] = rng.beta(p * k, (1 - p) * k)

    # Shared causals: MAF >= 1% in every ancestry, effects correlated rg across ancestries.
    # Unique causals: MAF >= 1% in their own ancestry, causal (non-zero effect) only there.
    maf = {x: np.minimum(freq[x], 1 - freq[x]) for x in ANCESTRIES}
    pool = np.flatnonzero(np.min(list(maf.values()), 0) >= 0.01)
    causal = {"shared": rng.choice(pool, a.num_shared_causals, replace=False)}
    for x in ANCESTRIES:
        pool = np.setdiff1d(np.flatnonzero(maf[x] >= 0.01), np.concatenate(list(causal.values())))
        causal[x] = rng.choice(pool, a.num_unique_causals[x], replace=False)
    cov = np.full((3, 3), a.rg) + (1 - a.rg) * np.eye(3)
    z_shared = rng.multivariate_normal(np.zeros(3), cov, (a.num_phenos, a.num_shared_causals))

    names = np.array([f"rs{c}_{b}" for c, b in zip(chrom, pos)])
    for k, x in enumerate(ANCESTRIES):
        idx = np.concatenate([causal["shared"], causal[x]])
        z = np.concatenate([z_shared[:, :, k], rng.standard_normal((a.num_phenos, len(causal[x])))], 1)
        p = freq["ANC"][idx]
        np.savetxt(out / f"effects_{x}.txt", z * (2 * p * (1 - p)) ** (a.power / 2), fmt="%.6g")
        (out / f"causals_{x}.txt").write_text((" ".join(names[idx]) + "\n") * a.num_phenos)
    with open(out / "causal_types.tsv", "w") as fh:
        fh.write("SNP\tTYPE\n" + "".join(f"{names[i]}\t{t}\n" for t, idx in causal.items() for i in idx))
    np.savez(out / "reference.npz", chrom=chrom, pos=pos, names=names, alleles=alleles,
             causal=np.sort(np.concatenate(list(causal.values()))), **{x: freq[x] for x in ANCESTRIES})


def per_ancestry(text):
    """'5' -> 5 for every ancestry; 'EUR=5,AFR=10' -> those counts, 0 for the rest."""
    if "=" not in text:
        return dict.fromkeys(ANCESTRIES, int(text))
    counts = dict.fromkeys(ANCESTRIES, 0)
    for item in text.split(","):
        name, value = item.split("=")
        if name.strip() not in counts:
            raise argparse.ArgumentTypeError(f"unknown ancestry {name!r}; use {ANCESTRIES}")
        counts[name.strip()] = int(value)
    return counts


def _chromosome(job):
    path, p, rho, block, n, seed = job
    rng = np.random.default_rng(seed)
    thr = norm_ppf(np.clip(p, 1e-12, 1 - 1e-12)).astype(np.float32)  # P(z < thr) = p
    z, noise = np.empty(2 * n, np.float32), np.empty(2 * n, np.float32)
    with open(path, "wb") as fh:
        for s in range(0, len(p), 256):
            g = np.empty((min(256, len(p) - s), n), np.int8)
            for i in range(len(g)):
                if s + i == 0 or rng.random() < 1 / block:  # new LD block
                    rng.standard_normal(dtype=np.float32, out=z)
                else:                                     # AR(1) along the chromosome
                    rng.standard_normal(dtype=np.float32, out=noise)
                    noise *= np.float32(math.sqrt(1 - rho ** 2))
                    z *= np.float32(rho)
                    z += noise
                hap = z < thr[s + i]
                g[i] = hap[:n]
                g[i] += hap[n:]  # int8 add (bool + bool would be OR)
            codes = np.pad(ENCODE[g], ((0, 0), (0, -n % 4))).reshape(len(g), -1, 4)
            fh.write((codes << SHIFTS).sum(2, dtype=np.uint8).tobytes())


def site(a):
    ref = np.load(Path(a.reference) / "reference.npz")
    rng = np.random.default_rng(a.seed)
    causal = ref["causal"]
    others = np.setdiff1d(np.arange(len(ref["names"])), causal)
    keep = np.sort(np.concatenate([causal, rng.choice(others, a.num_snps - len(causal), replace=False)]))
    chrom, p = ref["chrom"][keep], ref[a.ancestry][keep]
    _, _, rho, block = TREE[a.ancestry]

    stem, n = Path(a.out), a.num_samples
    stem.parent.mkdir(parents=True, exist_ok=True)
    seeds = np.random.SeedSequence([a.seed, ANCESTRIES.index(a.ancestry)]).spawn(22)
    jobs = [(f"{stem}.part{c}", p[chrom == c], rho, block, n, seeds[c - 1]) for c in range(1, 23)]
    with ProcessPoolExecutor(a.threads) as ex:
        list(ex.map(_chromosome, jobs))
    with open(f"{stem}.bed", "wb") as bed:
        bed.write(MAGIC)
        for job in jobs:
            with open(job[0], "rb") as part:
                shutil.copyfileobj(part, bed, 16 << 20)
            os.remove(job[0])

    pos = ref["pos"][keep]
    np.savetxt(f"{stem}.bim", np.column_stack([chrom, ref["names"][keep], pos / 1e6, pos, ref["alleles"][keep]]),
               fmt="%s", delimiter="\t")
    ids = np.char.add(f"{stem.name}_{a.ancestry}_", np.arange(1, n + 1).astype(str))
    sex, age = rng.integers(1, 3, n), rng.normal(55, 8, n).clip(40, 70).round(1)
    np.savetxt(f"{stem}.fam", np.column_stack([ids, ids, [0] * n, [0] * n, sex, [-9] * n]), fmt="%s")
    np.savetxt(f"{stem}.covar", np.column_stack([ids, ids, sex, age]), fmt="%s", header="FID IID sex age", comments="")


def read_table(path):
    """LDAK/PLINK phenotype or covariate file -> {IID: values} (NA -> nan)."""
    rows = [line.split() for line in open(path) if line.strip()]
    rows = rows[1:] if rows[0][0] in ("FID", "ID1") else rows
    return {r[1]: [np.nan if v == "NA" else float(v) for v in r[2:]] for r in rows}


def gwas(a):
    ids = np.array([line.split()[1] for line in open(f"{a.bfile}.fam")])
    bim = np.loadtxt(f"{a.bfile}.bim", dtype=str, ndmin=2)
    n, m = len(ids), len(bim)
    pheno, covar = read_table(a.pheno), read_table(a.covar)
    y = np.array([pheno.get(i, [np.nan])[0] for i in ids])
    keep = np.isfinite(y)
    C = np.column_stack([np.ones(keep.sum()), [covar[i] for i in ids[keep]]])
    Q = np.linalg.qr(C)[0]
    y = y[keep] - Q @ (Q.T @ y[keep])
    dof = len(y) - C.shape[1] - 1

    bed = np.memmap(f"{a.bfile}.bed", np.uint8, "r", 3, (m, (n + 3) // 4))
    beta, se, freq = np.empty(m), np.empty(m), np.empty(m)
    for s in range(0, m, 256):
        codes = (np.asarray(bed[s:s + 256])[..., None] >> SHIFTS) & 3
        x = DECODE[codes.reshape(len(codes), -1)[:, :n][:, keep]].T  # (samples, snps)
        freq[s:s + 256] = x.mean(0) / 2
        x -= x.mean(0)          # invariant SNPs become exactly 0 -> NaN results
        x -= Q @ (Q.T @ x)
        xx = (x * x).sum(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            b = (x.T @ y) / xx
            beta[s:s + 256], se[s:s + 256] = b, np.sqrt((y @ y - b * b * xx) / dof / xx)

    causal = set(Path(a.causals).read_text().split()) if a.causals else set()
    is_causal = np.isin(bim[:, 1], list(causal)).astype(int)
    cols = [bim[:, 1], bim[:, 0], bim[:, 3], bim[:, 4], bim[:, 5], freq.round(5), [keep.sum()] * m,
            beta, se, pvalue(beta / se), is_causal]
    np.savetxt(a.out, np.column_stack(cols), fmt="%s", delimiter="\t", comments="",
               header="SNP\tCHR\tBP\tA1\tA2\tA1_FREQ\tN\tBETA\tSE\tP\tCAUSAL")


def summarise(name, z, causal, types):
    ok = np.isfinite(z)
    hit = ok & (pvalue(np.nan_to_num(z)) < 5e-8)
    lam = np.median(z[ok & ~causal] ** 2) / 0.4549
    found = "".join(f"{f'{(hit & m).sum()}/{m.sum()}':>9}" for m in types.values())
    print(f"{name:<12}{hit.sum():>9}{(hit & causal).sum():>8}{lam:>10.3f}{found}")


def meta(a):
    tabs = [np.genfromtxt(f, names=True, dtype=None, encoding=None) for f in a.gwas]
    snps = np.unique(np.concatenate([t["SNP"] for t in tabs]))
    B, S = np.full((2, len(tabs), len(snps)), np.nan)
    causal = np.zeros(len(snps), bool)
    for k, t in enumerate(tabs):
        i = np.searchsorted(snps, t["SNP"])
        ok = np.minimum(t["A1_FREQ"], 1 - t["A1_FREQ"]) >= 0.01
        B[k, i[ok]], S[k, i[ok]] = t["BETA"][ok], t["SE"][ok]
        causal[i] |= t["CAUSAL"] == 1

    types = {}  # causal SNPs found at p < 5e-8, per causal type (shared / unique to an ancestry)
    if a.truth:
        snp, kind = np.array([line.split() for line in open(a.truth)][1:]).T
        causal |= np.isin(snps, snp)
        types = {t: np.isin(snps, snp[kind == t]) for t in dict.fromkeys(kind)}
    print(f"{'study':<12}{'GWS hits':>9}{'causal':>8}{'lambdaGC':>10}" + "".join(f"{t:>9}" for t in types))
    for k, f in enumerate(a.gwas):
        summarise(Path(f).name.split(".")[0], B[k] / S[k], causal, types)

    with np.errstate(divide="ignore", invalid="ignore"):
        W, b = np.nan_to_num(1 / S ** 2), np.nan_to_num(B)
        beta, se = (W * b).sum(0) / W.sum(0), 1 / np.sqrt(W.sum(0))
        q, df = (W * (b - beta) ** 2).sum(0), (W > 0).sum(0) - 1
        i2 = np.nan_to_num(np.clip((q - df) / q, 0, 1))
    p = pvalue(beta / se)
    summarise("META", beta / se, causal, types)
    print(f"mean I^2: causal SNPs {i2[causal & (df > 0)].mean():.2f}, other SNPs {i2[~causal & (df > 0)].mean():.2f}")
    np.savetxt(a.out, np.column_stack([snps, df + 1, beta, se, p, q, i2, causal.astype(int)]), fmt="%s",
               delimiter="\t", comments="", header="SNP\tK\tBETA\tSE\tP\tQ\tI2\tCAUSAL")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reference")
    r.add_argument("--out", required=True)
    r.add_argument("--num-snps", type=int, default=520000)
    r.add_argument("--num-shared-causals", type=int, default=20, help="causal in every ancestry")
    r.add_argument("--num-unique-causals", type=per_ancestry, default="0",
                   help="causal in one ancestry only: one number for each ancestry, or e.g. EUR=5,EAS=5,AFR=10")
    r.add_argument("--num-phenos", type=int, default=1)
    r.add_argument("--power", type=float, default=-0.25)
    r.add_argument("--rg", type=float, default=0.8, help="cross-ancestry effect correlation")
    r.add_argument("--seed", type=int, default=2024)
    s = sub.add_parser("site")
    s.add_argument("--reference", required=True)
    s.add_argument("--ancestry", required=True, choices=ANCESTRIES)
    s.add_argument("--num-samples", type=int, required=True)
    s.add_argument("--num-snps", type=int, required=True)
    s.add_argument("--seed", type=int, required=True)
    s.add_argument("--threads", type=int, default=os.cpu_count())
    s.add_argument("--out", required=True)
    g = sub.add_parser("gwas")
    for opt in ("--bfile", "--pheno", "--covar", "--out"):
        g.add_argument(opt, required=True)
    g.add_argument("--causals", help="reference causals_<ANC>.txt (adds the CAUSAL column)")
    t = sub.add_parser("meta")
    t.add_argument("--gwas", nargs="+", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--truth", help="reference causal_types.tsv: count found causal SNPs by type")
    a = ap.parse_args()
    {"reference": reference, "site": site, "gwas": gwas, "meta": meta}[a.cmd](a)


if __name__ == "__main__":
    main()