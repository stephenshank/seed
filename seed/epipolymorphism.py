import sys
import json
import argparse
from bisect import bisect_left, bisect_right

import pysam


def load_cpgs(input_cpgs, output_cpgs=None):
    cpg_table = open(input_cpgs)
    _ = next(cpg_table)
    cpg_hash = {}
    for line in cpg_table:
        _, chr, pos1, strand, _, _, _ = line.split('\t')
        pos0 = int(pos1) - 1
        if chr in cpg_hash:
            if strand in cpg_hash[chr]:
                cpg_hash[chr][strand].append(pos0)
            else:
                cpg_hash[chr][strand] = [pos0]
        else:
            cpg_hash[chr] = {strand: [pos0]}
    cpg_table.close()
    if output_cpgs:
        with open(output_cpgs, 'w') as json_file:
            json.dump(cpg_hash, json_file)
    return cpg_hash


def print_header():
    print('\t'.join([
        'query_name',
        'reference_start',
        'xm_tag',
        'cigartuples',
        'reference_end',
        'is_forward'
    ]))


def print_relevant_read_info(read):
    print('\t'.join([
        str(item)
        for item in [
            read.query_name,
            read.reference_start,
            read.get_tag('XM'),
            read.cigartuples,
            read.reference_end,
            read.is_forward
        ]
    ]))


def get_quartets(read, cpgs):
    pass
#    #print('cpgs: ', cpgs)
#    position = datum['reference_start']
#    #print('position: ', position)
#    bases = list(cpgs.base - position - 1)
#    #print('bases: ', bases, type(bases))
#    meths = []
#    base_index = 0
#    xm_tag_index = 0
#    ct_offset = 0
#    #print('xm tag enumerated:')
#    #for i, x in enumerate(xm_tag):
#    #    print(i, x)
#    dummy = ''
#    for ct in datum['cigartuples']:
#        #print(ct)
#        is_match = ct[0] == 0
#        is_insertion = ct[0] == 1
#        is_deletion = ct[0] == 2
#        if is_match:
#            dummy += xm_tag[ct_offset: ct_offset+ct[1]]
#            ct_offset += ct[1]
#        elif is_insertion:
#            ct_offset += ct[1]
#        elif is_deletion:
#            dummy += 'o'*ct[1]
#        else:
#            pass
#    print([dummy[base] for base in bases])
#    return None


def merge_quartets(all_quartets, quartets):
    pass


def write_quartets(all_quartets, output_bed):
    pass


def run_epipolymorphism(input_bam, input_cpgs, output_bed, loglevel=0):
    bam = pysam.AlignmentFile(input_bam, 'rb')
    cpgs = load_cpgs(input_cpgs)
    all_quartets = {}
    if loglevel > 0:
        print_header()
    for read in bam.fetch():
        if loglevel > 0:
            print_relevant_read_info(read)
        quartets = get_quartets(read, cpgs)
        merge_quartets(all_quartets, quartets)
    write_quartets(all_quartets, output_bed)


def cli():
    if len(sys.argv) == 1:
        print("For more information: seed --help")
        sys.exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b',
        '--bam',
        help='input (BAM)',
        dest='bam',
        required=True
    )
    parser.add_argument(
        '-c',
        '--cpgs',
        help='bismark CpGs (txt)',
        dest='cpgs',
        required=True
    )
    parser.add_argument(
        '-o',
        '--output',
        help='output (BED)',
        dest='output',
        required=True
    )
    parser.add_argument(
        '-l',
        '--log-level',
        help='log level (higher for more information)',
        dest='loglevel',
        type=int,
        default=0
    )
    args = parser.parse_args()
    run_epipolymorphism(args.bam, args.cpgs, args.output, args.loglevel)


if __name__ == '__main__':
    cli()
