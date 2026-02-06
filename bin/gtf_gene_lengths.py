#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict
import sys

def parse_args():
    ap = argparse.ArgumentParser(
        description="Extract gene info (gene_id, gene_name, gene_type, chr, start, end, strand) and merged exon length from a GTF."
    )
    ap.add_argument("--gtf", required=True, help="GTF annotation file (e.g., GENCODE).")
    ap.add_argument("--out", required=True, help="Output TSV file path.")
    ap.add_argument("--fallback_gene_span", action="store_true",
                    help="If a gene has no exon records, use gene feature span (end-start+1) as merged length.")
    return ap.parse_args()

def parse_gtf_attributes(attr_str):
    d = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        k, v = item.split(" ", 1)
        d[k] = v.strip().strip('"')
    return d

def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals, key=lambda x: x[0])
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            total += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s + 1)
    return total

def main():
    args = parse_args()

    gene_meta = {}   # gid -> (gene_name, gene_type, chr, strand)
    gene_span = {}   # gid -> (min_start, max_end)
    exons = defaultdict(list)

    with open(args.gtf, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            start = int(start); end = int(end)
            a = parse_gtf_attributes(attrs)

            gid = a.get("gene_id") or a.get("geneId") or a.get("gene")
            if not gid:
                continue

            if gid not in gene_meta:
                gname = a.get("gene_name") or a.get("Name") or ""
                gtype = a.get("gene_type") or a.get("gene_biotype") or ""
                gene_meta[gid] = (gname, gtype, chrom, strand)

            if feature == "gene":
                if gid not in gene_span:
                    gene_span[gid] = (start, end)
                else:
                    s0, e0 = gene_span[gid]
                    gene_span[gid] = (min(s0, start), max(e0, end))

            if feature == "exon":
                exons[gid].append((start, end))

    n_out = 0
    with open(args.out, "w") as fo:
        fo.write("GID\tSymbol\tGene_Type\tChr\tStart\tEnd\tStrand\tGene_Length\n")
        for gid, (gname, gtype, chrom, strand) in gene_meta.items():
            if exons.get(gid):
                merged_len = merge_intervals(exons[gid])
                s = min([x[0] for x in exons[gid]])
                e = max([x[1] for x in exons[gid]])
            else:
                if args.fallback_gene_span and gid in gene_span:
                    s, e = gene_span[gid]
                    merged_len = max(0, e - s + 1)
                else:
                    merged_len = 0
                    s, e = gene_span.get(gid, (0, 0))

            fo.write(f"{gid}\t{gname}\t{gtype}\t{chrom}\t{s}\t{e}\t{strand}\t{merged_len}\n")
            n_out += 1

    print(f"[done] Wrote {n_out} genes to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
