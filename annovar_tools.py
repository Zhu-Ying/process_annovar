#! /usr/bin/env python3
import argparse

from process_annovar import avinput_to_vcf, split_annovar_by_gene


def vcf_parser(subparsers):
    sub_parser = subparsers.add_parser('vcf', help='convert avinput to vcf')
    sub_parser.add_argument('--avinput', '-i', help='avinput infile')
    sub_parser.add_argument('--reference', '-r', help='reference fasta')
    sub_parser.add_argument('--vcf', '-o', help='vcf outfile')
    sub_parser.set_defaults(func=lambda args: avinput_to_vcf(avinput=args.avinput, reference=args.reference, vcf=args.vcf))


def split_parser(subparsers):
    sub_parser = subparsers.add_parser('split', help='split annovar result by gene')
    sub_parser.add_argument('--avoutput', '-i', help='avoutput infile')
    sub_parser.add_argument('--refgenes', '-r', action="append", help='refgene files')
    sub_parser.add_argument('--gene_db', '-g', help='the gene database, this execute will split annovar result by this gene_db')
    sub_parser.add_argument('--outfile', '-o', help='the split outfile')
    sub_parser.set_defaults(func=lambda args: split_annovar_by_gene(avoutput=args.avoutput, refgenes=args.refgenes, gene_db=args.gene_db, outfile=args.outfile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Process ANNOVAR result')
    subparsers = parser.add_subparsers(help='Process ANNOVAR result')
    vcf_parser(subparsers)
    split_parser(subparsers)
    args = parser.parse_args()
    print(args.refgenes)
    args.func(args)
