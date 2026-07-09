#!/usr/bin/env python
# coding: utf-8

"""
Python-based tool for that scans an NT/AA alignment and searchs for the insertions.

Parametrs:
    - alignment_file: input file (FASTA)
    - output_fasta: output file with the insertions (FASTA)
    - output_summary: output file, table with stats (TSV)
    - min_len: minimal insertion length (int)
    - gap_thr: gap rate threshold in a column (float)

Usage:
python insertions_search.py \
        -i FILENAME.fasta \
        -o_f FILENAME_insertions.fasta \
        -o_s FILENAME_summary.txt \
        --min_len 1 \
        --gap_thr 0.7
"""

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import argparse

def find_insertions(alignment_file, output_fasta, output_summary, min_len=1, gap_thr=0.7, domains_file=None, hmm_file=None):

#if domains_file exists (non-mandatory section)
    domains = {}
    if domains_file is not None:
        with open(domains_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.strip().split()
                if len(p) < 23:
                    continue
                try:
                    sid = p[3]
                    start = int(p[19])
                    end = int(p[20])
                except:
                    continue

                if sid not in domains:
                    domains[sid] = []
                domains[sid].append((start, end))

        for sid in domains:
            domains[sid].sort()

#if hmm_file exists (non-mandatory section)
    domain_names = []
    if hmm_file is not None:
        with open(hmm_file, "r") as f:
            for line in f:
                if line.startswith("NAME"):
                    name = line.strip().split()[1]
                    name = name.split("_i=")[0]
                    domain_names.append(name)

#Read the alignment
    aln = AlignIO.read(alignment_file, "fasta")
    arr = np.array([list(r.seq) for r in aln])

    found = []

#Scan for the insertions
    for i, rec in enumerate(aln):
        sid = rec.id
        L = arr.shape[1]
        pos = 0

        while pos < L:
            if arr[i, pos] != "-":

                other = np.concatenate([arr[:i, pos], arr[i+1:, pos]])
                if np.sum(other == "-") / len(other) >= gap_thr:
                    s = pos
                    internal_gap_count = 0
                    pos += 1

                    while pos < L:
                        other2 = np.concatenate([arr[:i, pos], arr[i+1:, pos]])
                        if np.sum(other2 == "-") / len(other2) >= gap_thr:
                            if arr[i, pos] == "-":
                                internal_gap_count += 1
                                if internal_gap_count > 500: #4     # можно вообще убрать. если внутри вставки есть гэп (например, в дупликации NS3 с убиквитином - считаем за одну вставку) 
                                    break
                            else:
                                internal_gap_count = 0

                            pos += 1
                        else:
                            break

                    e = pos
                    seq = "".join(arr[i, s:e]).replace("-", "")
                    real_len = len(seq)

                    if real_len >= min_len:
                        found.append((sid, seq, s + 1, e, real_len))

                else:
                    pos += 1
            else:
                pos += 1

#Group the insertions
    records = []
    for sid, seq, aln_start, aln_end, length in found:
        description = f"len={length}|pos={aln_start}-{aln_end}"
        records.append(SeqRecord(seq=Seq(seq), id=sid, description=description))
    SeqIO.write(records, output_fasta, "fasta")

#Cluster the insertions
    intervals = sorted([(s, e, l) for sid, seq, s, e, l in found], key=lambda x: x[0])
    clusters = []
    if intervals:
        cur_start, cur_end, l = intervals[0]
        count = 1

        for s, e, l in intervals[1:]:
            if s <= cur_end:  #перекрытие
                cur_end = max(cur_end, e)
                count += 1
            else:
                clusters.append((cur_start, cur_end, count))
                cur_start, cur_end = s, e
                count = 1

        clusters.append((cur_start, cur_end, count))

    merged_clusters = [c for c in clusters if c[2] > 1]

#The most variable cluster search
    most_variable = None
    if merged_clusters:
        max_density = 0
        for s, e, count in merged_clusters:
            region_len = e - s + 1
            density = count / region_len
            if density > max_density:
                max_density = density
                most_variable = (s, e, region_len, count, density)

    else:
        if found:
            longest = max(found, key=lambda x: x[4])
            s = longest[2]
            e = longest[3]
            l = longest[4]
            most_variable = (s, e, l, 1, 1.0)

#Save to output
    with open(output_summary, "w") as out:
        out.write("aln_start\taln_end\tlength\tsid\n")
        for sid, seq, s, e, l in found:
            out.write(f"{s}\t{e}\t{l}\t{sid}\n")

        if merged_clusters:
            out.write("\n# Insertions clusters\n")
            out.write("cluster_start\tcluster_end\tregion_length\tinsertions_count\n")
            for s, e, count in merged_clusters:
                region_len = e - s + 1
                out.write(f"{s}\t{e}\t{region_len}\t{count}\n")

        out.write("\n# Overall record\n")
        out.write(f"Insertions count \t{len(found)}\n")

        if merged_clusters:
            out.write(f"Clusters count \t{len(merged_clusters)}\n")

        if most_variable:
            s, e, region_len, count, density = most_variable
            out.write("\n# most variable region\n")
            out.write(f"start\t{s}\n")
            out.write(f"end\t{e}\n")
            out.write(f"length\t{region_len}\n")
            out.write(f"insertions\t{count}\n")
            out.write(f"density\t{density:.4f}\n")

    return found

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insertion search in an NT/AA alignment")
    parser.add_argument("-i", "--input", required=True, help="Input file (FASTA)")
    parser.add_argument("-o_f", "--output_fasta", required=True, help="Output file with the insertions (FASTA)")
    parser.add_argument("-o_s", "--output_summary", required=True, help="Output file, table with stats (TSV)")
    parser.add_argument("--min_len", type=int, default=1, help="Minimal insertion length")
    parser.add_argument("--gap_thr", type=float, default=0.7, help="Gap rate threshold in a column")
    parser.add_argument("--domains", help="output includes information from hmmscan table with domains, (non-mandatory)")
    parser.add_argument("--hmm", help="output includes information from HMM profiles (non-mandatory)")

    args = parser.parse_args()

    find_insertions(
        alignment_file=args.input,
        output_fasta=args.output_fasta,
        output_summary=args.output_summary,
        min_len=args.min_len,
        gap_thr=args.gap_thr,
        domains_file=args.domains,
        hmm_file=args.hmm)

