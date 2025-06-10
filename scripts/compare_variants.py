import argparse
import os
import pandas as pd
import vcfpy
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--cnvnator', nargs='+')
parser.add_argument('--gatk')
parser.add_argument('--delly', nargs='+')
parser.add_argument('--lumpy', nargs='+')
parser.add_argument('--manta', nargs='+')
parser.add_argument('--canvas', nargs='+')
parser.add_argument('--meerkat', nargs='+')
parser.add_argument('--output')
args = parser.parse_args()

WINDOW = 20

def parse_cnvnator(file):
    variants = []
    with open(file) as f:
        for line in f:
            if line.startswith('deletion') or line.startswith('duplication'):
                parts = line.strip().split()
                typ = 'DEL' if 'deletion' in line else 'DUP'
                chrom, rest = parts[1].split(':')
                start, end = map(int, rest.split('-'))
                variants.append((chrom, start, end, typ))
    return variants

def parse_meerkat(file):
    variants = []
    with open(file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2]) if parts[2].isdigit() else start
            svtype = parts[3] if len(parts) > 3 else 'INV'
            variants.append((chrom, start, end, svtype))
    return variants

def parse_vcf_files(files, source):
    records = []
    for f in files:
        reader = vcfpy.Reader.from_path(f)
        for rec in reader:
            if rec.INFO.get('SVTYPE') in ['DEL', 'DUP', 'INV']:
                records.append((rec.CHROM, rec.POS, rec.INFO.get('END', rec.POS), rec.INFO['SVTYPE'], source))
    return records

def build_variant_table(all_variants):
    bins = defaultdict(dict)
    for tool, entries in all_variants.items():
        for (chrom, start, end, svtype) in entries:
            bins[(chrom, start, svtype)][tool] = 1

    output = []
    for (chrom, pos, svtype), tools_called in bins.items():
        row = {
            'chrom': chrom,
            'pos': pos,
            'svtype': svtype,
            **{tool: int(tool in tools_called) for tool in all_variants}
        }
        output.append(row)

    return pd.DataFrame(output)

variant_calls = {
    'cnvnator': sum([parse_cnvnator(f) for f in args.cnvnator], []),
    'gatk': parse_vcf_files([args.gatk], 'gatk'),
    'delly': parse_vcf_files(args.delly, 'delly'),
    'lumpy': parse_vcf_files(args.lumpy, 'lumpy'),
    'manta': parse_vcf_files(args.manta, 'manta'),
    'canvas': parse_vcf_files(args.canvas, 'canvas'),
    'meerkat': sum([parse_meerkat(f) for f in args.meerkat], [])
}

df = build_variant_table(variant_calls)
df.sort_values(['chrom', 'pos'], inplace=True)
df.to_csv(args.output, sep='\t', index=False)
