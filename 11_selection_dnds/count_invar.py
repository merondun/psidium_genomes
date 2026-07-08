#!/usr/bin/env python3
import sys
from itertools import zip_longest

def read_fasta(path):
    seqs = []
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    seqs.append("".join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            seqs.append("".join(seq))
    return seqs

def count_sites(seqs):
    invariant = 0
    variant = 0
    aln_len = len(seqs[0])
    for col in zip(*seqs):
        bases = set(col) - set("-")   # ignore gaps
        if len(bases) == 1:
            invariant += 1
        else:
            variant += 1
    return invariant, variant

if __name__ == "__main__":
    for fa in sys.argv[1:]:
        seqs = read_fasta(fa)
        inv, var = count_sites(seqs)
        print(fa)
        print("Invariant sites:", inv)
        print("Variant sites:", var)
        print()
