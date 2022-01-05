#! /usr/bin/env python3
import argparse

from process_annovar import avinput_to_vcf, split_annovar_by_gene, check


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
    sub_parser.add_argument('--gene_based', '-g', help='the gene_based name of annovar, this execute will split annovar result by this gene_based')
    sub_parser.add_argument('--outfile', '-o', help='the split outfile')
    sub_parser.set_defaults(func=lambda args: split_annovar_by_gene(avoutput=args.avoutput, refgenes=args.refgenes, gene_db=args.gene_db, outfile=args.outfile))


def check_parser(subparsers):
    sub_parser = subparsers.add_parser('check', help='check avoutput format')
    sub_parser.add_argument('--input', '-i', help='avoutput infile')
    sub_parser.add_argument('--output', '-o', help='avoutput outfile')
    sub_parser.set_defaults(func=lambda args: check(input=args.input, output=args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('ANNOVAR tools')
    subparsers = parser.add_subparsers(help='ANNOVAR tools')
    vcf_parser(subparsers)
    split_parser(subparsers)
    check_parser(subparsers)
    args = parser.parse_args()
    args.func(args)
