import csv
import gzip

TRANS_TO_GENES: dict[str:set] = dict()
GENE_SYMBOL_TO_ID: dict[str:str] = dict()


def read_mane_select(mane_select: str) -> dict:
    trans_to_id = dict()
    with gzip.open(mane_select, 'rt') as fi:
        reader = csv.DictReader(fi, delimiter='\t')
        for row in reader:
            trans = row.get('RefSeq_nuc').split('.')[0]
            entrez_id = row.get('#NCBI_GeneID').split(':')[1]
            trans_to_id.setdefault(trans, entrez_id)
    return trans_to_id


def read_ncbi_gene_info(ncbi_gene_info: str) -> [dict, dict]:
    symbol_to_id = dict()
    synonyms_to_id = dict()
    with gzip.open(ncbi_gene_info, 'rt') as fi:
        reader = csv.DictReader(fi, delimiter='\t')
        for row in reader:
            symbol = row.get('Symbol')
            entrez_id = row.get('GeneID')
            synonyms = row.get('Synonyms').split('|')
            symbol_to_id.setdefault(symbol, entrez_id)
            for name in synonyms:
                synonyms_to_id.setdefault(name, entrez_id)
    return symbol_to_id, synonyms_to_id


def read_refgene(refgene: str):
    fi = open(refgene)
    for line in fi:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        trans, symbol = fields[1], fields[12]
        TRANS_TO_GENES.setdefault(trans, set()).add(symbol)
    fi.close()


def set_data(refgenes: list[str], ncbi_gene_info: str, mane_select: str):
    for refgene in refgenes:
        read_refgene(refgene)
    if ncbi_gene_info:
        symbol_to_id, synonyms_to_id = read_ncbi_gene_info(ncbi_gene_info)
    if mane_select:
        trans_to_id = read_mane_select(mane_select)
    if symbol_to_id or synonyms_to_id or trans_to_id:
        for trans, symbols in TRANS_TO_GENES.items():
            trans_short = trans.split('.')[0]
            for symbol in symbols:
                if GENE_SYMBOL_TO_ID.get(symbol):
                    continue
                entrez_id = trans_to_id.get(trans_short) or symbol_to_id.get(symbol) or synonyms_to_id.get(symbol) or '.'
                GENE_SYMBOL_TO_ID[symbol] = entrez_id
