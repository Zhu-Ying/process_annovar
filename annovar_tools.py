#! /usr/bin/env python3
import argparse

from process_annovar import avinput_to_vcf, split_annovar_by_gene, add_cnv_entrez_id, check


def vcf_parser(subparsers):
    sub_parser = subparsers.add_parser('vcf', help='convert avinput to vcf')
    sub_parser.add_argument('--avinput', '-i', help='avinput infile')
    sub_parser.add_argument('--reference', '-r', help='reference fasta')
    sub_parser.add_argument('--vcf', '-o', help='vcf outfile')
    sub_parser.set_defaults(func=lambda args: avinput_to_vcf(avinput=args.avinput, reference=args.reference, vcf=args.vcf))


def split_parser(subparsers):
    sub_parser = subparsers.add_parser('split', help='split annovar result by gene')
    sub_parser.add_argument('--avoutput', '-i', required=True, help='avoutput infile')
    sub_parser.add_argument('--gene_based', '-g', required=True, help='the gene_based name of annovar, this execute will split annovar result by this gene_based')
    sub_parser.add_argument('--outfile', '-o', required=True, help='the split outfile')
    sub_parser.add_argument('--refgenes', '-r', required=True, action="append", help='refgene files')
    sub_parser.add_argument('--mane_select', '-m', help='mane_select file')
    sub_parser.add_argument('--ncbi_gene_info', '-n', help='ncbi gene info file')
    sub_parser.set_defaults(func=lambda args: split_annovar_by_gene(
        avoutput=args.avoutput,
        gene_based=args.gene_based,
        outfile=args.outfile,
        refgenes=args.refgenes,
        ncbi_gene_info=args.ncbi_gene_info,
        mane_select=args.mane_select
    ))


def cnv_parser(subparsers):
    sub_parser = subparsers.add_parser('cnv', help='split annovar result by gene')
    sub_parser.add_argument('--avoutput', '-i', required=True, help='avoutput infile')
    sub_parser.add_argument('--colnames', '-c', required=True, action="append", help='the colnames  in annovar output need add entrez_id')
    sub_parser.add_argument('--outfile', '-o', required=True, help='the split outfile')
    sub_parser.add_argument('--refgenes', '-r', required=True, action="append", help='refgene files')
    sub_parser.add_argument('--mane_select', '-m', help='mane_select file')
    sub_parser.add_argument('--ncbi_gene_info', '-n', help='ncbi gene info file')
    sub_parser.set_defaults(func=lambda args: add_cnv_entrez_id(
        avoutput=args.avoutput,
        colnames=args.colnames,
        outfile=args.outfile,
        refgenes=args.refgenes,
        ncbi_gene_info=args.ncbi_gene_info,
        mane_select=args.mane_select
    ))


def check_parser(subparsers):
    sub_parser = subparsers.add_parser('check', help='check avoutput format')
    sub_parser.add_argument('--input', '-i', required=True, help='avoutput infile')
    sub_parser.add_argument('--output', '-o', required=True, help='avoutput outfile')
    sub_parser.set_defaults(func=lambda args: check(input=args.input, output=args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('ANNOVAR tools')
    subparsers = parser.add_subparsers(help='ANNOVAR tools')
    vcf_parser(subparsers)
    split_parser(subparsers)
    cnv_parser(subparsers)
    check_parser(subparsers)
    args = parser.parse_args()
    args.func(args)
