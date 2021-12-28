import re
import csv


def load_info(info: str) -> dict:
    data = dict()
    for field in info.split(';'):
        match = re.match(re.compile(r'^(\S+)=(\S+)$'))
        if match:
            data.setdefault(match.group(1), match.group(2))
    return data


def dump_info(data: dict) -> str:
    return ';'.join([f'{k}={v}' for k, v in data.items()])


def get_header(avoutput: str):
    with open(avoutput) as fi:
        return fi.readline().strip().split('\t')


def format_cnv_avouput(input: str, output: str):
    fieldnames = get_header(input)
    if 'Otherinfo1' not in fieldnames:
        raise Exception('ERROR: the Otherinfo1 not found')
    fi = open(input)
    fo = open(output, 'w')
    reader = csv.DictReader(fi, delimiter='\t')
    writer = csv.DictWriter(fo, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    for row in reader:
        info = load_info(row['Otherinfo1'])
        ref, alt = info.pop('ref', None), info.pop('alt', None)
        if not ref or not alt:
            raise Exception(
                f'ERROR: ref or alt not found in otherinfo: {row["Chr"]}:{row["Start"]}-{row["End"]} {row["Otherinfo1"]}'
            )
        row.update({'Ref': ref, 'Alt': alt})
        writer.writerow(row)
    fi.close()
    fo.close()
