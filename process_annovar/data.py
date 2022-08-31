import pandas as pd

TRANS_TO_GENE: dict[str:set] = dict()
GENE_SYMBOL_TO_ID: dict[str:str] = dict()
SYMBOL_TO_HGNC_NAME: dict[str:str] = dict()


def read_gene2refseq(gene2refseq: str) -> dict:
    usecols = ['#tax_id', 'RNA_nucleotide_accession.version', 'GeneID']
    df = pd.read_table(gene2refseq, compression='gzip', usecols=usecols, dtype={"#tax_id": str, "GeneID": str})
    df.rename(columns={'#tax_id': 'TaxID', 'RNA_nucleotide_accession.version': 'Transcript'}, inplace=True)
    df = df.query('TaxID == "9606"')
    df.Transcript = df.Transcript.apply(lambda x: x.split('.')[0])
    return df.set_index(['Transcript'])['GeneID'].to_dict()


def read_ncbi_gene_info(ncbi_gene_info: str) -> [dict, dict]:
    symbol_to_id = dict()
    synonyms_to_id = dict()
    df = pd.read_table(ncbi_gene_info, compression='gzip', usecols=['GeneID', 'Symbol', 'Synonyms'], dtype={"GeneID": str})
    for row in df.iloc:
        symbol_to_id.setdefault(row.Symbol, row.GeneID)
        for name in row.get('Synonyms').split('|'):
            synonyms_to_id.setdefault(name, row.GeneID)
    return symbol_to_id, synonyms_to_id


def read_refgene(refgene: str) -> dict:
    df = pd.read_table(refgene, header=None)
    df.columns = [
        'Bin', 'Name', 'Chrom', 'Strand', 'TxStart', 'TxEnd', 'CdsStart', 'CdsEnd',
        'ExonCount', 'ExonStarts', 'ExonEnds', 'Score', 'Gene', 'CdsStartStat', 'CdsEndStat', 'ExonFrames'
    ]
    TRANS_TO_GENE.update(df.set_index(['Name'])['Gene'].to_dict())


def read_gene_hgnc_name(gene_hgnc_name: str) -> dict:
    df = pd.read_table(gene_hgnc_name)
    for row in df.iloc:
        symbol, name = row.get('symbol'), row.get('name')
        SYMBOL_TO_HGNC_NAME[symbol] = name


def set_data(refgenes: list[str], ncbi_gene_info: str, gene2refseq: str, gene_hgnc_name: str):
    [read_refgene(refgene) for refgene in refgenes]
    if ncbi_gene_info:
        symbol_to_id, synonyms_to_id = read_ncbi_gene_info(ncbi_gene_info)
    if gene2refseq:
        trans_to_id = read_gene2refseq(gene2refseq)
    for trans, symbol in TRANS_TO_GENE.items():
        trans_short = trans.split('.')[0]
        entrez_id = trans_to_id.get(trans_short) or symbol_to_id.get(symbol) or synonyms_to_id.get(symbol)
        if entrez_id:
            GENE_SYMBOL_TO_ID.setdefault(symbol, entrez_id)
    if gene_hgnc_name:
        read_gene_hgnc_name(gene_hgnc_name)
