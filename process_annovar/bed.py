from collections import namedtuple

BEDCnv = namedtuple('VCFSnv', ['chrom', 'start', 'end', 'type'])


def parse_info(info: str) -> dict:
    data = dict()
    for item in info.split(';'):
        if item.find('=') != -1:
            key, val = item.split('=')[0:2]
            data.setdefault(key.upper(), val)
    return data


def read_avinput(infile: str) -> list:
    cnvs = list()
    fi = open(infile)
    for line in fi:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        info = fields[5] if len(fields) > 5 else '.'
        alt = parse_info(info).get('ALT')
        if alt:
            for typo in alt.split('/'):
                cnvs.append(BEDCnv(chrom=fields[0], start=int(fields[1]), end=int(fields[2]), type=typo))
    return cnvs


def avinput_to_bed(avinput: str, bed: str):
    bed_cnvs = read_avinput(avinput)
    fo = open(bed, 'w')
    for bed_cnv in bed_cnvs:
        fo.write(f'{bed_cnv.chrom}\t{bed_cnv.start}\t{bed_cnv.end}\t{bed_cnv.type}\n')
    fo.close()
