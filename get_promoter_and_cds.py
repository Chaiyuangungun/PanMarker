#!/usr/bin/env python
import sys


def load_fa(fa_db, in_fa):
    with open(in_fa, 'r') as fin:
        id = ""
        seq = ""
        for line in fin:
            if line[0] == '>':
                if seq != "":
                    fa_db[id] = seq
                id = line.strip().split()[0][1:]
                seq = ""
            else:
                seq += line.strip().upper()
    fa_db[id] = seq


def reverse_seq(seq):
    rev_seq = ""
    base_db = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for base in seq[::-1]:
        if base in base_db:
            rev_seq += base_db[base]
        else:
            rev_seq += base
    return rev_seq


def extract_promoter_and_cds(in_gff3, in_fa, prm_len, out_pre):
    print("Loading genome")
    fa_db = {}
    load_fa(fa_db, in_fa)

    print("Loading and extracting promoter and cds")
    out_cds = out_pre + '.cds'
    out_promoter = out_pre + '.prm'
    fcds = open(out_cds, 'w')
    fprm = open(out_promoter, 'w')

    with open(in_gff3, 'r') as fin:
        region_db = {}
        for line in fin:
            if line.strip() == '' or line[0] ==  '#':
                continue
            data = line.strip().split()
            if data[2] == 'gene':
                for info in data[8].split(';'):
                    if info.startswith("Name"):
                        id = info.split('=')[1]
                chrn = data[0]
                if id not in region_db:
                    region_db[id] = {'chrn': chrn, 'rna': [], 'best_score': 0, 'direct': "", 'cds': []}

            elif data[2] == 'mRNA':
                score = 1
                for info in data[8].split(';'):
                    if info.startswith('cov') or info.startswith('iden'):
                        score *= float(info.split('=')[1])
                if score > region_db[id]['best_score']:
                    region_db[id]['chrn'] = data[0]
                    region_db[id]['best_score'] = score
                    region_db[id]['rna'] = [int(data[3]), int(data[4])]
                    region_db[id]['direct'] = data[6]
                    region_db[id]['cds'] = []
            elif data[2] == 'CDS':
                if score == region_db[id]['best_score']:
                    region_db[id]['cds'].append([int(data[3]), int(data[4])])

        for id in sorted(region_db):
            chrn = region_db[id]['chrn']
            sp, ep = region_db[id]['rna']
            if region_db[id]['direct'] == '+':
                psp = sp-prm_len-1
                if psp >= 0:
                    fprm.write(">%s::%d::%d\n%s\n"%(id, psp+1, sp-1, fa_db[chrn][psp: sp]))
                    cds = ""
                    for sp, ep in region_db[id]['cds']:
                        cds += fa_db[chrn][sp-1: ep]
                    fcds.write('>%s\n%s\n'%(id, cds))
            else:
                pep = ep+prm_len
                if pep <= len(fa_db[chrn]):
                    fprm.write(">%s::%d::%d\n%s\n"%(id, ep, pep-1, fa_db[chrn][ep: pep]))
                    cds = ""
                    for sp, ep in region_db[id]['cds']:
                        cds += reverse_seq(fa_db[chrn][sp-1: ep])
                    fcds.write('>%s\n%s\n'%(id, cds))
    
    print("Finished")


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python %s <in_gff3> <in_fa> <promotor_len> <out_prefix>"%sys.argv[0])
    else:
        in_gff3, in_fa, prm_len, out_pre = sys.argv[1:]
        prm_len = int(float(prm_len))
        extract_promoter_and_cds(in_gff3, in_fa, prm_len, out_pre)
