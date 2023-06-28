#!/usr/bin/env python
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


def get_synteny_track_from_anchors(
        anchors, bed1, bed2, genome1, genome2, name,
        track_id=""
        ):
    """
    create synteny track json from anchors and bed files

    """
    if track_id == "":
        # calculate unique track id as hash of arguments
        track_id = hash(
            (anchors, bed1, bed2, genome1, genome2, name)
            )
        track_id = F"{name}_{track_id}"

    json_template = {
        "type": "SyntenyTrack",
        "trackId": track_id,
        "name": name,
        "assemblyNames": [
            genome1,
            genome2
            ],
        "adapter": {
            "type": "MCScanAnchorsAdapter",
            "mcscanAnchorsLocation": {
                "locationType": "UriLocation",
                "name": name,
                "uri": anchors
                },
            "bed1Location": {
                "locationType": "UriLocation",
                "name": bed1,
                "uri": bed1
                },
            "bed2Location": {
                "locationType": "UriLocation",
                "name": bed2,
                "uri": bed2
                },
            "assemblyNames": [
                genome1,
                genome2
                ]
            }
        }
    return json_template


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description='script fo adding synteny track to jbrowse config file'
            )
    parser.add_argument(
            '-c', '--config', required=True,
            help='json config file which should be updated'
            )
    parser.add_argument(
            '-t', '--table', required=True,
            help='tab delimited table with synteny information'
            )
    parser.add_argument(
            '-o', '--output', required=True,
            help='output configuration json file'
            )
    parser.add_argument(
            '-d', '--dir', default=".", help="jbrowse directory"
            )
    parser.add_argument(
            '-w', '--overwrite', action="store_true", default=False,
            help="overwrite existing files"
            )

    args = parser.parse_args()

    track_info = read_table(args.table)
    print(track_info)

    with open(args.config) as f:
        config = json.load(f)

    # iterate over track info and add renamed data to args.dir
    for track in track_info:
        new_track = {}
        new_track['filename'] = F'{track["name"]}.{track["type"]}'
        # copy to dir with new name
        new_fn = args.dir + "/" + new_track['filename']
        if args.overwrite or not os.path.exists(new_fn):
            shutil.copy(track['filename'], new_fn)
        if track['bed1'] != "":
            new_track['bed1'] = F'{track["assembly1"]}_{track["name"]}.bed'
            new_fn = args.dir + "/" + new_track['bed1']
            if args.overwrite or not os.path.exists(new_fn):
                shutil.copy(track['bed1'], new_fn)
        if track['bed2'] != "":
            new_track['bed2'] = F'{track["assembly2"]}_{track["name"]}.bed'
            new_fn = args.dir + "/" + new_track['bed2']
            if args.overwrite or not os.path.exists(new_fn):
                shutil.copy(track['bed2'], new_fn)

        if track['type'] == "anchors":
            json_track = get_synteny_track_from_anchors(
                anchors=new_track['filename'],
                bed1=new_track['bed1'],
                bed2=new_track['bed2'],
                genome1=track['assembly1'],
                genome2=track['assembly2'],
                name=track['name']
                )
            config['tracks'].append(json_track)

    # export new json config
    with open(args.output, 'w') as f:
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    main()
