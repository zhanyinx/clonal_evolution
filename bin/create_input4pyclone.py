#!/usr/bin/env python
"""Build the PyClone-VI input table for one patient (one or more tumour samples).

Rules (one row per SNV per tumour sample):
- Variants: SNVs only (MAF Variant_Type == SNP), taken from every sample's MAF;
  MAF columns are selected by header name, never by position.
- Allele counts: the MAF t_ref_count / t_alt_count of the sample where the
  variant was called. For a sample whose MAF lacks the variant, the counts are
  re-counted from its CRAM/BAM with `samtools mpileup -B -A -Q20 -q20 -d0`
  (BAQ off, no depth cap); read-start (^x) and read-end ($) marks and indel
  sequences are stripped before counting.
- Purity: tumour_content = the sample's ASCAT purity file, which the ASCAT step
  forces to the samplesheet cellularity when one is given. A purity that is
  missing, 0 or outside (0, 1] stops the run.
- Copy number: the ASCAT segment that contains the position (segment
  coordinates are 1-based inclusive, as are MAF positions). A position that no
  segment contains takes, in order: the flanking segments' state when both
  flanks agree; else the nearest segment on the same chromosome within
  --max_gap bp. Otherwise the variant is dropped for the whole patient (never
  a silent 1/1). Positions with major_cn == 0 are dropped too.
A variant dropped in any sample is dropped from every sample of the patient.
Every decision is written to --report.
"""

import argparse
import math
import os
import re
import subprocess
import sys

import pandas as pd

from utils import read_maf

MAF_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "t_ref_count",
    "t_alt_count",
]
PYCLONE_COLUMNS = [
    "mutation_id",
    "sample_id",
    "ref_counts",
    "alt_counts",
    "normal_cn",
    "major_cn",
    "minor_cn",
    "tumour_content",
]
MPILEUP_ARGS = ["-B", "-A", "-Q", "20", "-q", "20", "-d", "0"]
_INDEL = re.compile(r"[+-](\d+)")


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samples", required=True,
                   help="TSV with header: sample_id, maf, alignment, cnvs, purityploidy, sex")
    p.add_argument("--fasta", default=None, help="Reference fasta (needed to re-count from CRAM)")
    p.add_argument("--output", required=True, help="PyClone-VI input TSV")
    p.add_argument("--report", required=True, help="Per variant x sample decision report TSV")
    p.add_argument("--max_gap", type=int, default=1_000_000,
                   help="Max distance (bp) to the nearest ASCAT segment for an uncovered position")
    p.add_argument("--samtools", default="samtools")
    return p.parse_args()


def norm_chrom(c) -> str:
    c = str(c)
    return c if c.startswith("chr") else "chr" + c


def read_snvs(maf_path: str) -> pd.DataFrame:
    """SNVs of one MAF, columns selected by header name."""
    maf = read_maf(maf_path)
    missing = [c for c in MAF_COLUMNS if c not in maf.columns]
    if missing:
        raise ValueError(f"MAF {maf_path} lacks columns {missing}")
    maf = maf.loc[maf["Variant_Type"] == "SNP", MAF_COLUMNS].copy()
    maf["Chromosome"] = maf["Chromosome"].map(norm_chrom)
    maf["Start_Position"] = maf["Start_Position"].astype(int)
    maf["key"] = list(zip(maf["Chromosome"], maf["Start_Position"],
                          maf["Reference_Allele"], maf["Tumor_Seq_Allele2"]))
    if maf["key"].duplicated().any():
        raise ValueError(f"MAF {maf_path} has duplicated SNVs")
    return maf.set_index("key")


def count_bases(bases: str, ref: str, alt: str):
    """(ref, alt) counts from an mpileup base string; ^x, $ and indel sequences stripped."""
    ref, alt = ref.upper(), alt.upper()
    n_ref = n_alt = 0
    i = 0
    while i < len(bases):
        b = bases[i]
        if b == "^":  # read start: skip the mapping-quality char
            i += 2
            continue
        if b == "$":
            i += 1
            continue
        if b in "+-":
            m = _INDEL.match(bases, i)
            i = m.end() + int(m.group(1))
            continue
        if b in ".,":
            n_ref += 1
        elif b.upper() == alt:
            n_alt += 1
        elif b.upper() == ref:  # mpileup without a reference prints the base itself
            n_ref += 1
        i += 1
    return n_ref, n_alt


def recount(alignment: str, fasta, keys, samtools: str, workdir: str, tag: str):
    """Re-count (ref, alt) at SNV keys from an alignment, BAQ off, no depth cap."""
    out = {k: (0, 0) for k in keys}
    if not keys:
        return out
    pos_file = os.path.join(workdir, f"{tag}.recount.pos")
    with open(pos_file, "w") as fh:
        for c, p in sorted({(k[0], k[1]) for k in keys}):
            fh.write(f"{c}\t{p}\n")
    cmd = [samtools, "mpileup", *MPILEUP_ARGS, "-l", pos_file]
    if fasta:
        cmd += ["-f", fasta]
    cmd.append(alignment)
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    piles = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        piles[(norm_chrom(f[0]), int(f[1]))] = f[4]
    for k in keys:
        bases = piles.get((k[0], k[1]))
        if bases is not None:
            out[k] = count_bases(bases, k[2], k[3])
    return out


def read_segments(path: str) -> pd.DataFrame:
    seg = pd.read_csv(path, sep="\t")
    seg.columns = ["chr", "startpos", "endpos", "major_cn", "minor_cn"]
    seg["chr"] = seg["chr"].map(norm_chrom)
    return seg


def assign_cn(seg: pd.DataFrame, chrom: str, pos: int, max_gap: int):
    """(major, minor, source) for a 1-based position, or (None, None, reason)."""
    s = seg[seg["chr"] == chrom]
    hit = s[(s["startpos"] <= pos) & (s["endpos"] >= pos)]
    if len(hit) > 1:
        return None, None, "overlapping_segments"
    if len(hit) == 1:
        r = hit.iloc[0]
        return int(r.major_cn), int(r.minor_cn), "segment"
    if s.empty:
        return None, None, "no_segment_on_chromosome"
    left = s[s["endpos"] < pos].sort_values("endpos").tail(1)
    right = s[s["startpos"] > pos].sort_values("startpos").head(1)
    if len(left) and len(right):
        l, r = left.iloc[0], right.iloc[0]
        if (l.major_cn, l.minor_cn) == (r.major_cn, r.minor_cn):
            return int(l.major_cn), int(l.minor_cn), "flanks_agree"
    cand = []
    if len(left):
        cand.append((pos - left.iloc[0].endpos, left.iloc[0]))
    if len(right):
        cand.append((right.iloc[0].startpos - pos, right.iloc[0]))
    dist, best = min(cand, key=lambda x: x[0])
    if dist <= max_gap:
        return int(best.major_cn), int(best.minor_cn), "nearest_segment"
    return None, None, "no_segment_within_max_gap"


def read_purity(path: str) -> float:
    pp = pd.read_csv(path, sep="\t")
    val = float(pp["AberrantCellFraction"].values[0])
    if math.isnan(val) or not (0 < val <= 1):
        raise ValueError(f"Purity {val} in {path} is not in (0, 1]; stopping.")
    return val


def main():
    args = _parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    need = ["sample_id", "maf", "alignment", "cnvs", "purityploidy", "sex"]
    if [c for c in need if c not in samples.columns]:
        raise ValueError(f"--samples needs columns {need}")
    if samples["sample_id"].duplicated().any():
        raise ValueError("duplicated sample_id in --samples")

    mafs = {r.sample_id: read_snvs(r.maf) for r in samples.itertuples()}
    keys = sorted({k for m in mafs.values() for k in m.index})
    genes = {}
    for m in mafs.values():
        for k, g in m["Hugo_Symbol"].items():
            genes.setdefault(k, g)

    workdir = os.path.dirname(os.path.abspath(args.output))
    rows, report = [], []
    for r in samples.itertuples():
        if r.sex not in ("XX", "XY"):
            raise ValueError(f"sex must be XX or XY, sample {r.sample_id} has {r.sex}")
        purity = read_purity(r.purityploidy)
        seg = read_segments(r.cnvs)
        maf = mafs[r.sample_id]
        absent = [k for k in keys if k not in maf.index]
        rec = recount(r.alignment, args.fasta, absent, args.samtools, workdir, r.sample_id)
        for k in keys:
            if k in maf.index:
                ref_c, alt_c = int(maf.at[k, "t_ref_count"]), int(maf.at[k, "t_alt_count"])
                csrc = "maf"
            else:
                (ref_c, alt_c), csrc = rec[k], "recount"
            major, minor, cn_src = assign_cn(seg, k[0], k[1], args.max_gap)
            if major == 0:
                major, minor, cn_src = None, None, "major_cn_0"
            normal_cn = 1 if (r.sex == "XY" and k[0] in ("chrX", "chrY")) else 2
            mid = f"{k[0]}:{k[1]}:{k[2]}:{k[3]}"
            rows.append(dict(mutation_id=mid, sample_id=r.sample_id, ref_counts=ref_c, alt_counts=alt_c,
                             normal_cn=normal_cn, major_cn=major, minor_cn=minor, tumour_content=purity))
            report.append(dict(mutation_id=mid, gene=genes[k], chrom=k[0], pos=k[1], ref=k[2], alt=k[3],
                               sample_id=r.sample_id, count_source=csrc, ref_counts=ref_c, alt_counts=alt_c,
                               major_cn=major, minor_cn=minor, cn_source=cn_src, tumour_content=purity))

    out = pd.DataFrame(rows, columns=PYCLONE_COLUMNS)
    rep = pd.DataFrame(report)
    dropped = set(out.loc[out["major_cn"].isna(), "mutation_id"])
    rep["in_pyclone"] = ~rep["mutation_id"].isin(dropped)
    out = out[~out["mutation_id"].isin(dropped)].copy()
    out[["major_cn", "minor_cn"]] = out[["major_cn", "minor_cn"]].astype(int)
    if out.empty:
        raise ValueError("No SNV left for PyClone after the copy-number rule.")
    out.to_csv(args.output, sep="\t", index=False)
    rep.to_csv(args.report, sep="\t", index=False)
    print(f"{len(keys)} SNVs, {len(dropped)} dropped by the copy-number rule, "
          f"{out['mutation_id'].nunique()} to PyClone across {len(samples)} sample(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
