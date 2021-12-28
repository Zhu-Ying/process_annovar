import csv
import re

from collections import namedtuple
from typing import Union

GeneAnno = namedtuple('GeneAnno', ['gene', 'region', 'detail', 'event'])
Snv = namedtuple('Snv', ['chrom', 'start', 'end', 'ref', 'alt'])
TRANS_TO_GENES: dict[str:set] = dict()


def split_gene_anno(func: str, gene: str, exonic_func, gene_detail, aa_change) -> list:
    funcs = func.split(';') if func else []
    genes = gene.split(';') if gene else []
    gene_details = re.split(r',|;', gene_detail) if gene_detail else []
    aa_changes = re.split(r',|;', aa_change) if aa_change else []
    gene_anno_dict = dict()
    in_gene = bool(set(funcs) - {'intergenic', 'upstream', 'downstream'})
    if in_gene:
        if len(genes) >= len(funcs):
            if exonic_func != '.':
                funcs += ['exonic'] * (len(genes) - len(funcs))
            else:
                funcs += [funcs[0]] * (len(genes) - len(funcs))
        else:
            raise Exception('the gene number lt the func number')
        for i in range(len(genes)):
            region, gene = funcs[i], genes[i]
            if region.find('exonic') != -1:
                event = exonic_func
                detail = ','.join([ac for ac in aa_changes if ac.split(':')[0] == gene]) or '.'
            else:
                event = 'splicing' if region.find('splic') != -1 else '.'
                tmp_details = list()
                for tmp_detail in gene_details:
                    if tmp_detail != '.':
                        trans_id = tmp_detail.split(':')[0]
                        genes_set = TRANS_TO_GENES.get(trans_id)
                        tmp_gene = list((genes_set & set(genes) or genes_set))[0]
                        if tmp_gene == gene:
                            tmp_details.append(f'{gene}:{tmp_detail}')
                detail = ','.join(tmp_details)
            if gene_anno_dict.get(gene):
                old: GeneAnno = gene_anno_dict.get(gene)
                if old.region.find(region) == -1:
                    region = f'{old.region},{region}'
                if old.event.find(event) == -1:
                    event = f'{old.event},{event}'
                if old.detail.find(detail) == -1:
                    detail = f'{old.detail},{detail}'
            gene_anno_dict[gene] = GeneAnno(gene=gene, region=region, detail=detail, event=event)
    else:
        gene_anno_dict.setdefault('.', GeneAnno(gene='.', region='Non-gene', detail='.', event='.'))
    return list(gene_anno_dict.values())


def read_refgene(refgene: str):
    fi = open(refgene)
    for line in fi:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        TRANS_TO_GENES.setdefault(fields[1], set()).add(fields[12])
    fi.close()


def parse_row(row: dict, gene_db: str) -> Union[Snv, dict, dict]:
    snv = Snv(chrom=row.get('Chr'), start=row.get('Start'), end=row.get('End'), ref=row.get('Ref'), alt=row.get('Alt'))
    info = dict()
    for key, val in row.items():
        if key in ['Chr', 'Start', 'End', 'Ref', 'Alt']:
            continue
        keys = key.split('.')
        if len(keys) > 1 and keys[0] in ['Func', 'Gene', 'ExonicFunc', 'GeneDetail', 'AAChange']:
            continue
        info[key] = val
    func, gene, exonic_func = row.get(f'Func.{gene_db}', ''), row.get(f'Gene.{gene_db}', ''), row.get(f'ExonicFunc.{gene_db}', '')
    gene_detail, aa_change = row.get(f'GeneDetail.{gene_db}', ''), row.get(f'AAChange.{gene_db}', '')
    gene_annos = split_gene_anno(func=func, gene=gene, exonic_func=exonic_func, gene_detail=gene_detail, aa_change=aa_change)
    return snv, gene_annos, info


def split_annovar_by_gene(avoutput: str, refgenes: list[str], gene_db: str, outfile: str):
    for refgene in refgenes:
        read_refgene(refgene)
    fi = open(avoutput)
    fo = open(outfile, 'w')
    reader = csv.DictReader(fi, delimiter='\t')
    head = 'Chr\tStart\tEnd\tRef\tAlt\tGene\tEvent\tRegion\tDetail\t'
    info_keys = list()
    for row in reader:
        snv, gene_annos, info = parse_row(row, gene_db)
        if not info_keys:
            info_keys = list(info.keys())
            head += '\t'.join(info_keys)
            fo.write(f'{head}\n')
        info_text = '\t'.join([info.get(key, '.') for key in info_keys])
        for gene_anno in gene_annos:
            fo.write(f'{snv.chrom}\t{snv.start}\t{snv.end}\t{snv.ref}\t{snv.alt}\t'
                     f'{gene_anno.gene}\t{gene_anno.event}\t{gene_anno.region}\t{gene_anno.detail}\t{info_text}\n')
    fi.close()
