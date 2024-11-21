#!/usr/bin/env python3

import sys

def merge_paf_intervals(input_file, output_file, tolerance_percent, min_intervals):
    """
    Merges PAF alignments by chaining intervals into larger ones based on the same order
    and similar distances in both genome assemblies within a specified tolerance.
    PAF file is sorted based on qname and qstart before processing.
    Coordinates are handled according to the PAF specification.
    Chains with fewer intervals than 'min_intervals' are discarded.

    Parameters:
    - input_file: Path to the input PAF file.
    - output_file: Path to the output PAF file with merged alignments.
    - tolerance_percent: Allowed percentage difference between distances to consider for merging.
    - min_intervals: Minimum number of merged intervals required to keep a chain.
    """
    tolerance_fraction = tolerance_percent / 100.0

    # Read PAF records into a list
    records = []
    with open(input_file, 'r') as infile:
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            fields = line.strip().split('\t')
            if len(fields) < 12:
                sys.stderr.write(f"Warning: Line {line_num} does not have at least 12 fields. Skipping.\n")
                continue

            # Extract required fields
            try:
                qname = fields[0]
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                strand = fields[4]
                tname = fields[5]
                tlen = int(fields[6])
                tstart = int(fields[7])
                tend = int(fields[8])
                nmatch = int(fields[9])
                alen = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                sys.stderr.write(f"Warning: Line {line_num} has invalid numerical fields. Skipping.\n")
                continue

            # Ensure qstart ≤ qend and tstart ≤ tend (coordinates are always increasing)
            qstart, qend = min(qstart, qend), max(qstart, qend)
            tstart, tend = min(tstart, tend), max(tstart, tend)

            record = {
                'qname': qname,
                'qlen': qlen,
                'qstart': qstart,
                'qend': qend,
                'strand': strand,
                'tname': tname,
                'tlen': tlen,
                'tstart': tstart,
                'tend': tend,
                'nmatch': nmatch,
                'alen': alen,
                'mapq': mapq,
                'line_num': line_num,
                'line': line,  # Keep the original line in case needed
                'intervals': 1  # Each record starts as one interval
            }
            records.append(record)

    # Sort records based on qname and qstart
    records.sort(key=lambda r: (r['qname'], r['qstart']))

    with open(output_file, 'w') as outfile:
        current_chain = None  # Stores the current chain alignment
        chain_intervals = 1    # Counts the number of intervals in the current chain

        for i, record in enumerate(records):
            qname = record['qname']
            qlen = record['qlen']
            qstart = record['qstart']
            qend = record['qend']
            strand = record['strand']
            tname = record['tname']
            tlen = record['tlen']
            tstart = record['tstart']
            tend = record['tend']
            nmatch = record['nmatch']
            alen = record['alen']
            mapq = record['mapq']
            line_num = record['line_num']

            # For the first alignment
            if current_chain is None:
                current_chain = {
                    'qname': qname,
                    'qlen': qlen,
                    'qstart': qstart,
                    'qend': qend,
                    'strand': strand,
                    'tname': tname,
                    'tlen': tlen,
                    'tstart': tstart,
                    'tend': tend,
                    'nmatch': nmatch,
                    'alen': alen,
                    'mapq': mapq,
                    'prev_qstart': qstart,
                    'prev_tstart': tstart,
                    'prev_tend': tend,
                    'intervals': 1  # Start with one interval
                }
                chain_intervals = 1
                continue

            # Check if we can merge current alignment with the current chain
            can_merge = False

            if (qname == current_chain['qname'] and
                tname == current_chain['tname'] and
                strand == current_chain['strand']):

                delta_query = qstart - current_chain['prev_qstart']

                if strand == '+':
                    delta_target = tstart - current_chain['prev_tstart']
                    sign_consistent = delta_query * delta_target > 0
                elif strand == '-':
                    delta_target = current_chain['prev_tend'] - tend  # As query increases, target decreases
                    sign_consistent = delta_query * delta_target > 0
                else:
                    sys.stderr.write(f"Warning: Line {line_num} has unrecognized strand '{strand}'. Skipping.\n")
                    continue

                if sign_consistent:
                    # Compute mean delta and difference
                    mean_delta = (abs(delta_query) + abs(delta_target)) / 2.0 if (delta_query != 0 or delta_target != 0) else 0
                    difference = abs(abs(delta_query) - abs(delta_target))
                    allowed_difference = mean_delta * tolerance_fraction

                    if mean_delta == 0 or difference <= allowed_difference:
                        # Intervals can be merged
                        can_merge = True

            if can_merge:
                # Update current chain
                current_chain['qstart'] = min(current_chain['qstart'], qstart)
                current_chain['qend'] = max(current_chain['qend'], qend)
                current_chain['nmatch'] += nmatch
                current_chain['alen'] += alen
                current_chain['mapq'] = min(current_chain['mapq'], mapq)
                current_chain['prev_qstart'] = qstart

                current_chain['tstart'] = min(current_chain['tstart'], tstart)
                current_chain['tend'] = max(current_chain['tend'], tend)
                current_chain['prev_tstart'] = tstart
                current_chain['prev_tend'] = tend

                current_chain['intervals'] += 1  # Increment interval count
                chain_intervals += 1
            else:
                # Output current chain if it meets the minimum intervals requirement
                if current_chain['intervals'] >= min_intervals:
                    output_fields = [
                        current_chain['qname'],
                        str(current_chain['qlen']),
                        str(current_chain['qstart']),
                        str(current_chain['qend']),
                        current_chain['strand'],
                        current_chain['tname'],
                        str(current_chain['tlen']),
                        str(current_chain['tstart']),
                        str(current_chain['tend']),
                        str(current_chain['nmatch']),
                        str(current_chain['alen']),
                        str(current_chain['mapq'])
                    ]
                    outfile.write('\t'.join(output_fields) + '\n')
                else:
                    # Optionally, you can print a message or count discarded chains
                    pass

                # Start a new chain
                current_chain = {
                    'qname': qname,
                    'qlen': qlen,
                    'qstart': qstart,
                    'qend': qend,
                    'strand': strand,
                    'tname': tname,
                    'tlen': tlen,
                    'tstart': tstart,
                    'tend': tend,
                    'nmatch': nmatch,
                    'alen': alen,
                    'mapq': mapq,
                    'prev_qstart': qstart,
                    'prev_tstart': tstart,
                    'prev_tend': tend,
                    'intervals': 1  # Start with one interval
                }
                chain_intervals = 1

        # Output the last chain if it meets the minimum intervals requirement
        if current_chain is not None and current_chain['intervals'] >= min_intervals:
            output_fields = [
                current_chain['qname'],
                str(current_chain['qlen']),
                str(current_chain['qstart']),
                str(current_chain['qend']),
                current_chain['strand'],
                current_chain['tname'],
                str(current_chain['tlen']),
                str(current_chain['tstart']),
                str(current_chain['tend']),
                str(current_chain['nmatch']),
                str(current_chain['alen']),
                str(current_chain['mapq'])
            ]
            outfile.write('\t'.join(output_fields) + '\n')
        else:
            # Optionally, you can print a message or count discarded chains
            pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Merge PAF alignments by chaining intervals based on similar distances within a tolerance.')
    parser.add_argument('input_paf', help='Input PAF file')
    parser.add_argument('output_paf', help='Output PAF file with merged alignments')
    parser.add_argument('-t', '--tolerance', type=float, default=0.0,
                        help='Tolerance percentage for distance differences (default: 0.0)')
    parser.add_argument('-m', '--min-intervals', type=int, default=1,
                        help='Minimum number of merged intervals required to keep a chain (default: 1)')

    args = parser.parse_args()

    merge_paf_intervals(args.input_paf, args.output_paf, args.tolerance, args.min_intervals)
