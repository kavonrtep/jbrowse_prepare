#!/usr/bin/env python
import hashlib

import os

import shutil

import argparse
import csv
import json
from typing import Dict, List


def read_table(table_file) -> List[Dict[str, str]]:
    """
    read table file and return list of dictionaries
    columns labels are keys of dictionaries
    :param table_file: table file
    :return: list of dictionaries
    """
    table = []
    with open(table_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row: Dict[str, str]
        for row in reader:
            table.append(row)
    return table

def get_synteny_track_from_paf(genome1, genome2, paf):
    track_id = hashlib.md5((genome1 + genome2 + paf).encode()).hexdigest()  # to get unique track id
    track_id = F"{genome1}_{genome2}_{track_id}"
    json_template =  {
            "type": "SyntenyTrack",
            "trackId": track_id,
            "name": paf,
            "assemblyNames": [
                genome1,
                genome2
            ],
            "adapter": {
                "type": "PAFAdapter",
                "targetAssembly": genome2,
                "queryAssembly": genome1,
                "pafLocation": {
                    "locationType": "UriLocation",
                    "uri": F"./{paf}"
                }
            },
            "category": ["Synteny"],
        }
    return json_template

def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description='script for adding synteny track to jbrowse config file'
            )
    parser.add_argument(
            '-c', '--config', required=True,
            help='json config file which should be updated'
            )
    parser.add_argument(
            '-t', '--table', required=True,
            help=('tab delimited table with synteny information.'
                  ' This file mu contain header with three values:'
                  ' GENOME1, GENOME2, PAF.')
            )
    parser.add_argument(
            '-o', '--output', required=True,
            help='output configuration json file'
            )
    parser.add_argument(
            '-d', '--dir', default=None, help="jbrowse directory"
            )

    args = parser.parse_args()

    track_info = read_table(args.table)
    print(track_info)

    with open(args.config) as f:
        config = json.load(f)

    # iterate over PAF in table and create json config for each
    for track in track_info:
        # assume all tracks are already in args.dir
        # check if track PAF file exists
        # only check if arg.dir is provided
        if args.dir:
            if not os.path.exists(os.path.join(args.dir, track["PAF"])):
                raise FileNotFoundError(F"File {track['PAF']} not found in {args.dir}")

        json_track = get_synteny_track_from_paf(track["GENOME1"],
                                                track["GENOME2"],
                                                track["PAF"])
        config["tracks"].append(json_track)
    # export to json
    with open(args.output, "w") as f:
        json.dump(config, f, indent=4)


if __name__ == "__main__":
    main()