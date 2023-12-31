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
        'is_forward',
        'mapping_quality'
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
            read.is_forward,
            read.mapping_quality
        ]
    ]))


def get_quartets(read, cpgs, mapping_quality=30):
    if read.mapping_quality < mapping_quality:
        return None
    strand = 'F' if read.is_forward else 'R'
    try:
        cpgs_on_chromosome = cpgs[read.reference_name][strand]
    except KeyError:
        error_message = 'warning: key error with %s' % read.reference_name
        print(error_message, file=sys.stderr)
        return None
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
    if quartets is None:
        return
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


all_methylation = list(
    it.product((False, True), (False, True), (False, True), (False, True))
)


def get_methylation_counts(quartets_at_position, minimum_count=10):
    row = []
    for methylation in all_methylation:
        if methylation in quartets_at_position:
            row.append(quartets_at_position[methylation])
        else:
            row.append(0)
    count = sum(row)
    if count < minimum_count:
        return None
    freqs = [ x / count for x in row ]
    epipolymorphism = 1 - sum([f**2 for f in freqs])
    return row + [count] + freqs + [epipolymorphism]


def stringify(row):
    return [str(i) for i in row]


def write_quartets(chromosome, all_quartets, output_tsv_file):
    positions = sorted(
        all_quartets.keys(),
        key=lambda x: x[0]
    )
    for position in positions:
        row_position = [chromosome] + [str(i) for i in position]
        row_counts = get_methylation_counts(all_quartets[position])
        if not row_counts is None:
            output_tsv_file.write('\t'.join(sum([
                stringify(row)
                for row in [row_position, row_counts]
            ], [])) + '\n')


def run_epipolymorphism(input_bam, input_cpgs, output_tsv_path,
        output_quartets=None, loglevel=0):
    bam = pysam.AlignmentFile(input_bam, 'rb')
    output_tsv_file = open(output_tsv_path, 'w')
    counts_header = [
        ''.join([
            'Z' if meth_base else 'z' for meth_base in methylation
        ])
        for methylation in all_methylation
    ]
    freqs_header = [count_header + '_freq' for count_header in counts_header]
    header = '\t'.join([
        'chromosome', 'cpg1_0', 'cpg2_0', 'cpg3_0', 'cpg4_0'
    ] + counts_header + ['count'] + freqs_header + ['epipolymorphism']) + '\n'
    output_tsv_file.write(header)
    cpgs = load_cpgs(input_cpgs)
    all_quartets = {}
    if loglevel > 0:
        print_header()
    previous_chromosome = None
    for read in bam.fetch():
        chromosome = read.reference_name
        if previous_chromosome == None:
            previous_chromosome = chromosome
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
