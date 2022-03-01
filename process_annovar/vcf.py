from pysam import FastaFile
from collections import namedtuple

AVSnv = namedtuple('AVSnv', ['id', 'chrom', 'start', 'end', 'ref', 'alt', 'info'])
VCFSnv = namedtuple('VCFSnv', ['chrom', 'pos', 'ref', 'alt'])


def fetch_seq(fasta: FastaFile, chrom: str, start: int, end: int):
    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'
    if chrom == 'chrMT':
        chrom = 'chrM'
    return fasta.fetch(reference=chrom, start=start, end=end)


def recovery_ins(av_snv: AVSnv, fasta: FastaFile) -> VCFSnv:
    pos = int(av_snv.start)
    ref = fetch_seq(fasta=fasta, chrom=av_snv.chrom, start=pos - 1, end=pos)
    alt = f'{ref}{av_snv.alt}'
    return VCFSnv(chrom=av_snv.chrom, pos=pos, ref=ref, alt=alt)


def recovery_del(av_snv: AVSnv, fasta: FastaFile) -> VCFSnv:
    pos = int(av_snv.start) - 1
    alt = fetch_seq(fasta=fasta, chrom=av_snv.chrom, start=pos - 1, end=pos)
    ref = f'{alt}{av_snv.ref}'
    return VCFSnv(chrom=av_snv.chrom, pos=pos, ref=ref, alt=alt)


def recovery_snp(av_snv: AVSnv) -> VCFSnv:
    return VCFSnv(chrom=av_snv.chrom,
                  pos=av_snv.start,
                  ref=av_snv.ref,
                  alt=av_snv.alt)


def read_avinput(infile: str) -> list:
    snvs = list()
    fi = open(infile)
    for line in fi:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        info = fields[5] if len(fields) > 5 else '.'
        snvs.append(
            AVSnv(chrom=fields[0],
                  start=int(fields[1]),
                  end=int(fields[2]),
                  ref=fields[3],
                  alt=fields[4],
                  id=':'.join(fields[0:5]),
                  info=info))
    return snvs


def avinput_to_vcf(avinput: str, reference: str, vcf: str):
    fasta = FastaFile(reference)
    av_snvs = read_avinput(avinput)
    fo = open(vcf, 'w')
    fo.write('##fileformat=VCFv4.0\n')
    fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for av_snv in av_snvs:
        if av_snv.ref == '-':
            vcf_snv = recovery_ins(av_snv, fasta)
        elif av_snv.alt == '-':
            vcf_snv = recovery_del(av_snv, fasta)
        else:
            vcf_snv = recovery_snp(av_snv)
        fo.write(
            f'{vcf_snv.chrom}\t{vcf_snv.pos}\t{av_snv.id}\t{vcf_snv.ref}\t{vcf_snv.alt}\t.\t.\t{av_snv.info}\n'
        )
    fo.close()
