import csv
import re
from collections import namedtuple
from .data import TRANS_TO_GENE, GENE_SYMBOL_TO_ID, set_data
GeneAnno = namedtuple('GeneAnno', ['gene', 'entrez_id', 'region', 'detail', 'event'])
Snv = namedtuple('Snv', ['chrom', 'start', 'end', 'ref', 'alt'])


def get_gene_details(gene_detail, aa_change) -> dict:
    gene_details = re.split(r',|;', gene_detail) if gene_detail else []
    aa_changes = re.split(r',|;', aa_change) if aa_change else []
    trans_ids = set()
    detail_dict = dict()
    for detail in gene_details:
        if detail != '.' and ':' in detail:
            trans_id = detail.split(':')[0]
            gene = TRANS_TO_GENE.get(trans_id)
            if trans_id not in trans_ids:
                detail_dict.setdefault(gene, list()).append(detail)
            trans_ids.add(trans_id)
    for detail in aa_changes:
        if detail != '.' and ':' in detail:
            gene, trans_id = detail.split(':')[0:2]
            if trans_id not in trans_ids:
                detail_dict.setdefault(gene, list()).append(detail)
            trans_ids.add(trans_id)
    return detail_dict


def split_gene_anno(func: str, gene: str, exonic_func, gene_detail, aa_change) -> list:
    funcs = func.split(';') if func else []
    genes = gene.split(';') if gene else []
    detail_dict = get_gene_details(gene_detail, aa_change)
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
            entrez_id = GENE_SYMBOL_TO_ID.get(gene, '.')
            detail = ','.join(detail_dict.get(gene, []))
            is_exonic_splicing = re.findall(r'c.\d+[+-]\d+_\d+del|c.\d+_\d+[+-]\d+del', detail)
            if is_exonic_splicing:
                event = 'splicing'
                region = 'exonic_splicing'
            elif region.find('exonic') != -1:
                event = exonic_func
            else:
                event = 'splicing' if region.find('splic') != -1 else '.'
            if region.startswith('exon') and detail == '.':
                region, event, detail = '.', '.', '.'
            if gene_anno_dict.get(gene):
                old: GeneAnno = gene_anno_dict.get(gene)
                if old.region.find(region) == -1 and old.region != '.':
                    region = f'{old.region},{region}'
                if old.event.find(event) == -1 and old.event != '.':
                    event = f'{old.event},{event}'
                if old.detail.find(detail) == -1 and old.detail != '.':
                    detail = f'{old.detail},{detail}'
                else:
                    detail = old.detail
            gene_anno_dict[gene] = GeneAnno(gene=gene, entrez_id=entrez_id, region=region, detail=detail, event=event)
    else:
        gene_anno_dict.setdefault('.', GeneAnno(gene='.', entrez_id='.', region='Non-gene', detail='.', event='.'))
    return list(gene_anno_dict.values())


def parse_row(row: dict, gene_based: str) -> [Snv, dict, dict]:
    snv = Snv(chrom=row.get('Chr'), start=row.get('Start'), end=row.get('End'), ref=row.get('Ref'), alt=row.get('Alt'))
    info = dict()
    for key, val in row.items():
        if key in ['Chr', 'Start', 'End', 'Ref', 'Alt']:
            continue
        keys = key.split('.')
        if len(keys) > 1 and keys[0] in ['Func', 'Gene', 'ExonicFunc', 'GeneDetail', 'AAChange']:
            continue
        info[key] = val
    func, gene, exonic_func = row.get(f'Func.{gene_based}', ''), row.get(f'Gene.{gene_based}', ''), row.get(f'ExonicFunc.{gene_based}', '')
    gene_detail, aa_change = row.get(f'GeneDetail.{gene_based}', ''), row.get(f'AAChange.{gene_based}', '')
    gene_annos = split_gene_anno(func=func, gene=gene, exonic_func=exonic_func, gene_detail=gene_detail, aa_change=aa_change)
    return snv, gene_annos, info


def split_annovar_by_gene(avoutput: str, gene_based: str, outfile: str, refgenes: list[str], ncbi_gene_info: str = None, gene2refseq: str = None):
    set_data(refgenes=refgenes, ncbi_gene_info=ncbi_gene_info, gene2refseq=gene2refseq)
    fi = open(avoutput)
    fo = open(outfile, 'w')
    reader = csv.DictReader(fi, delimiter='\t')
    head = 'Chr\tStart\tEnd\tRef\tAlt\tGene\tEntrezID\tEvent\tRegion\tDetail\t'
    info_keys = list()
    for row in reader:
        snv, gene_annos, info = parse_row(row, gene_based)
        if not info_keys:
            info_keys = list(info.keys())
            head += '\t'.join(info_keys)
            fo.write(f'{head}\n')
        info_text = '\t'.join([info.get(key, '.') for key in info_keys])
        for gene_anno in gene_annos:
            fo.write(f'{snv.chrom}\t{snv.start}\t{snv.end}\t{snv.ref}\t{snv.alt}\t'
                     f'{gene_anno.gene}\t{gene_anno.entrez_id}\t'
                     f'{gene_anno.event}\t{gene_anno.region}\t{gene_anno.detail}\t{info_text}\n')
    fi.close()
    fo.close()
