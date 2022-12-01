#!/usr/bin/env python
''' parse table with tracks for jbrowse2 and create config.json file '''
import shutil

import json

import subprocess

import argparse

import os
import sys
import csv
from shlex import quote



def read_config(config_file, directory):
    """
    read csv table, return list of tracks
    """
    input_tracks = []
    print('reading config file:')
    with open(config_file, 'r') as f:
        reader = csv.DictReader(f, delimiter="\t")
        labels = set()
        first = True
        for record in reader:
            if first:
                assembly_name = record['label']
                first = False
            if record['label'] in labels:
                print('duplicated label:', record['label'])
                sys.exit('labels must be unique, exiting...')
            labels.add(record['label'])
            record['new_filename'] = os.path.join(
                directory, F"{assembly_name}_{record['label']}.{record['format']}"
                )
            record['src_filename'] = os.path.join(
                record['dirname'], record['filename']
                )
            input_tracks.append(record)
    print('')
    return input_tracks



def create_reference(ref_path, outdir, assembly_name):
    """
    create conf.json file using jbrowse add-assembly command
    """
    print('creating reference:')
    # check if index exists
    indexfile = F"{ref_path}.fai"
    if not os.path.exists(indexfile):
        # create index of fasta reference
        cmd = F"samtools faidx {quote(ref_path)}"
        subprocess.check_call(cmd, shell=True)
    cmd = (
        F"jbrowse add-assembly {quote(ref_path)} --out {outdir} --load inPlace --name "
        F"{assembly_name} --overwrite"
    )
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    print('')



def get_config_string(track):
    """ use additional information from track to
    to create config string to set color, height, etc
    setting is specific base on the track format"""
    config = ''
    color = track['color'] if track['color'] else 'goldenrod'


    if track['format'] == 'bigwig':
        config = {"displays": [{"displayId": track['label'], "type":
            "LinearWiggleDisplay",
                                "renderers": {"XYPlotRenderer": {"color": color,
                                    "type": "XYPlotRenderer"}}}]}
    if track['format'] == 'bed' and track['type'] != 'arcs':
        config = {"displays": [{"displayId": track['label'], "type": "LinearBasicDisplay",
                                "renderer": {"color1": color,
                                             "type": "SvgFeatureRenderer"}}]}
    if track['format'] == 'gff3':
        config = {"displays": [{"displayId": track['label'], "type": "LinearBasicDisplay",
                                "renderer": {"color1": color,
                                             "type": "SvgFeatureRenderer"}}]}
    if track['format'] == 'bed' and track['type'] == 'arcs':
        config = {"displays": [{"displayId": track['label'], "type": "LinearArcDisplay",
                                "renderer": {"color": color,
                                             "type": "ArcRenderer",
                                             "label": ""}}]}

    if config!= '':
        # escape quotes in config string
        config_str = json.dumps(config).replace('"', '\\"')
        cfg_cmd = "--config \"" + config_str  + "\""
    else:
        cfg_cmd = ''

    return cfg_cmd



def add_track(track, outdir, assembly_name):
    """
    add track to jbrowse using jbrowse add-track command
    """
    print('adding track:')
    # check is csi index exists, is so use it in add-track command
    indexfile = F"{track['new_filename']}.csi"
    cfg_cmd = get_config_string(track)
    if os.path.exists(indexfile):
        # use index:
        cmd = (F"jbrowse add-track {quote(track['new_filename'])} --out {outdir} "
               F"--load inPlace --trackId {quote(track['label'])} --indexFile"
               F" {quote(indexfile)} "
               F" --category {quote(track['category'])}"
               F" {cfg_cmd} --assemblyNames {quote(assembly_name)} --overwrite")
        print(cmd)
    else:
        cmd = (F"jbrowse add-track {quote(track['new_filename'])} --out {outdir} "
               F"--load inPlace --trackId {quote(track['label'])} "
               F" --category {quote(track['category'])} --overwrite"
               F" {cfg_cmd} --assemblyNames {quote(assembly_name)}")
        print(cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        print(F"Error adding track {track['src_filename']}, ")
        print("\n\n")

    print('')

def sort_and_index(track):
    """
    sort and index track, handles bed,  gff and bam files, bigwig files are already sorted
    use always CSI index for every file
    """
    # check if sortet.gz and csi index exists, in not create them

    sortfilegz = F"{track['new_filename']}.sorted.{track['format']}.gz"
    sortfile = F"{track['new_filename']}.sorted.{track['format']}"
    indexfile = F"{track['new_filename']}.sorted.{track['format']}.csi"
    indexfilegz = F"{track['new_filename']}.sorted.{track['format']}.gz.csi"
    # BED
    if track['format'] == 'bed':
        if os.path.exists(sortfilegz) and os.path.exists(indexfilegz):
            print(F"sorted and indexed file {sortfilegz} exists, skipping...")
        else:
            cmd = F"sort -k1,1 -k2,2n {quote(track['new_filename'])} > {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"bgzip {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"tabix -C -p bed {quote(sortfilegz)}"
            subprocess.check_call(cmd, shell=True)
        return sortfilegz
    # GFF
    if track['format'] == 'gff' or track['format'] == 'gff3':
        if os.path.exists(sortfilegz) and os.path.exists(indexfilegz):
            print(F"sorted and indexed file {sortfilegz} exists, skipping...")
        else:
            cmd = F"sort -k1,1 -k4,4n {quote(track['new_filename'])} > {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"bgzip {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"tabix -C -p gff {quote(sortfilegz)}"
            subprocess.check_call(cmd, shell=True)
        return sortfilegz
    # BAM
    if track['format'] == 'bam':
        if os.path.exists(sortfile) and os.path.exists(indexfile):
            print(F"sorted and indexed file {sortfile} exists, skipping...")
        else:
            cmd = F"samtools sort {quote(track['new_filename'])} -o {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"samtools index -c {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
        return sortfile
    if track['format'] == 'bigwig':
        return track['new_filename']
    print('format not recognized, skipping sorting and indexing')
    print(track['format'])
    print('-------------------')
    return track['new_filename']


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--output_dir', type=str, default="data")
    parser.add_argument("-c", "--config_table", required=True)
    parser.add_argument("-C", "--config_path", required=False, default="config.yaml")
    parser.add_argument("-u", "--update", action="store_true", default=False)
    args = parser.parse_args()

    track_list = read_config(args.config_table, args.output_dir)
    # first track is always genome reference
    ref_path = track_list[0]['new_filename']
    print('reference path:', ref_path)

    # first remove old config.json file
    if os.path.exists(os.path.join(args.output_dir, 'config.json')) and not args.update:
        print('removing old config.json file')
        os.remove(os.path.join(args.output_dir, 'config.json'))

    # copy all files to output directory but only if they are not already there
    assembly_name = track_list[0]['label']





    for track in track_list:
        if not os.path.exists(track['new_filename']):
            print('copying file:', track['src_filename'])
            print('to:', track['new_filename'])
            # use shutil.copy2 to keep file metadata
            shutil.copy(track['src_filename'], track['new_filename'])

        else:
            print('file already exists:', track['new_filename'])

    print('creating reference files:')

    create_reference(ref_path, args.output_dir, assembly_name)
    for track in track_list[1:]:
        # sort and make index:
        track['new_filename'] = sort_and_index(track)
        print('adding track:', track['label'])
        add_track(track, args.output_dir, assembly_name)
        print("----------------------------------------")
        print("----------------------------------------")
        print("----------------------------------------")

