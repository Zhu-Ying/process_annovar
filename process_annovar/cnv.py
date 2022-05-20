import csv
from .data import GENE_SYMBOL_TO_ID, set_data

csv.field_size_limit(500 * 1024 * 1024)


def add_cnv_entrez_id(avoutput: str, colnames: list[str], outfile: str, refgenes: list[str], ncbi_gene_info: str = None, gene2refseq: str = None):
    set_data(refgenes=refgenes, ncbi_gene_info=ncbi_gene_info, gene2refseq=gene2refseq)
    fi = open(avoutput)
    fo = open(outfile, 'w')
    reader = csv.DictReader(fi, delimiter='\t')
    writer = csv.DictWriter(fo, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    for row in reader:
        for colname in colnames:
            if not row.get(colname) or row.get(colname) == '.':
                continue
            exon_annos = row.get(colname).replace('Name=', '').split(',')
            for i in range(len(exon_annos)):
                symbol = exon_annos[i].split(':')[0]
                entrez_id = GENE_SYMBOL_TO_ID.get(symbol, '.')
                exon_annos[i] = f'{entrez_id}:{exon_annos[i]}'
            row[colname] = 'Name=%s' % ','.join(exon_annos)
        writer.writerow(row)
    fi.close()
    fo.close()
