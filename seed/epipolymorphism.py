import sys
import json
import argparse
import itertools as it
from bisect import bisect_left, bisect_right

import pysam


def load_cpgs(input_cpgs, output_cpgs=None):
    cpg_table = open(input_cpgs)
    _ = next(cpg_table)
    cpg_hash = {}
    for line in cpg_table:
        _, chromosome, pos1, strand, _, _, _ = line.split('\t')
        pos0 = int(pos1) - 1
        if chromosome in cpg_hash:
            if strand in cpg_hash[chromosome]:
                cpg_hash[chromosome][strand].append(pos0)
            else:
                cpg_hash[chromosome][strand] = [pos0]
        else:
            cpg_hash[chromosome] = {strand: [pos0]}
    cpg_table.close()
    for chromosome in cpg_hash.keys():
        for strand in cpg_hash[chromosome].keys():
            cpg_hash[chromosome][strand].sort()
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
    strand = 'F' if read.is_forward else 'R'
    cpgs_on_chromosome = cpgs[read.reference_name][strand]
    cpg_start = bisect_left(cpgs_on_chromosome, read.reference_start)
    cpg_end = bisect_right(cpgs_on_chromosome, read.reference_end-1)
    bases = [
        cpg - read.reference_start
        for cpg in cpgs_on_chromosome[cpg_start: cpg_end]
    ]
    ct_offset = 0
    matches = ''
    xm_tag = read.get_tag('XM')
    for ct in read.cigartuples:
        is_match = ct[0] == 0
        is_insertion = ct[0] == 1
        is_deletion = ct[0] == 2
        if is_match:
            matches += xm_tag[ct_offset: ct_offset+ct[1]]
            ct_offset += ct[1]
        elif is_insertion:
            ct_offset += ct[1]
        elif is_deletion:
            matches += 'o'*ct[1]
        else:
            sys.stderr.write('no logic for %s!' % read.query_name)
    quartets = []
    indices = range(4)
    for i in range(len(bases) - 3):
        quartet = tuple(
            bases[i + j] + read.reference_start for j in indices
        ) + tuple(
            matches[bases[i + j]] == 'Z' for j in indices
        )
        quartets.append(quartet)
    return quartets


def merge_quartets(all_quartets, quartets):
    for quartet in quartets:
        positions = quartet[:4]
        methylation = quartet[4:]
        if positions in all_quartets:
            if methylation in all_quartets[positions]:
                all_quartets[positions][methylation] += 1
            else:
                all_quartets[positions][methylation] = 1
        else:
            all_quartets[positions] = {
                methylation: 1
            }


def write_quartets(chromosome, all_quartets, output_tsv_file):
    all_methylation = list(
        it.product((False, True), (False, True), (False, True), (False, True))
    )
    for position in all_quartets.keys(): 
        row_start = [chromosome] + [str(i) for i in position]
        row_end = []
        for methylation in all_methylation:
            if methylation in all_quartets[position]:
                row_end.append(str(all_quartets[position][methylation]))
            else:
                row_end.append('0')
        output_tsv_file.write('\t'.join(row_start + row_end) + '\n')


def run_epipolymorphism(input_bam, input_cpgs, output_tsv_path, loglevel=0):
    bam = pysam.AlignmentFile(input_bam, 'rb')
    output_tsv_file = open(output_tsv_path, 'w')
    all_methylation = list(
        it.product((False, True), (False, True), (False, True), (False, True))
    )
    header = '\t'.join([
        'chromosome', 'cpg1_0', 'cpg2_0', 'cpg3_0', 'cpg4_0'
    ] + [
        ''.join([
            'Z' if meth_base else 'z' for meth_base in methylation
        ])
        for methylation in all_methylation
    ]) + '\n'
    output_tsv_file.write(header)
    cpgs = load_cpgs(input_cpgs)
    all_quartets = {}
    if loglevel > 0:
        print_header()
    previous_chromosome = None
    for read in bam.fetch():
        chromosome = read.reference_name
        if loglevel > 0:
            print_relevant_read_info(read)
        quartets = get_quartets(read, cpgs)
        merge_quartets(all_quartets, quartets)
        if previous_chromosome != chromosome:
            write_quartets(chromosome, all_quartets, output_tsv_file)
            all_quartets = {}
        previous_chromosome = chromosome
    write_quartets(chromosome, all_quartets, output_tsv_file)
    output_tsv_file.close()


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
