#!/usr/bin/env python
"""
convert bedpe to bed file for jbrowse arcs vidualization
"""
from math import log10

import argparse
import csv

argparser = argparse.ArgumentParser(
    description="convert bedpe to bed file for jbrowse arcs vidualization"
    )
argparser.add_argument(
    "-i", "--input", help="input bedpe file", required=True
    )
argparser.add_argument(
    "-o", "--output", help="output bed file", required=True
    )

args = argparser.parse_args()

with open(args.input, "r") as bedpe_file, open(args.output, "w") as bed_file:
    bedpe_reader = csv.reader(bedpe_file, delimiter="\t")
    bed_writer = csv.writer(bed_file, delimiter="\t")
    # iterate over bedpe file, get also i as counter

    for i, row in enumerate(bedpe_reader):
        # each line of bedpe is converted to two lines of bed
        # verify that start is smaler than end
        # must be on the same chromosome!!
        if row[0] != row[3]:
            print(row)
            print("chromosomes in bedpe line are not the same")
            print("skipping line")
            continue
        coord1 = [int(row[1]), int(row[2]), int(row[4]), int(row[5])]
        coord1.sort()
        name1 = "{}_{}".format(i, row[6])
        # output bed line is chr start end name
        bed_writer.writerow(
            [row[0], coord1[0], coord1[2], name1, -log10(
                float(
                    row[7]
                    )
                )]
            )
        bed_writer.writerow(
            [row[0], coord1[1], coord1[3], name1, -log10(
                float(
                    row[7]
                    )
                )]
            )

# Path: bedp
