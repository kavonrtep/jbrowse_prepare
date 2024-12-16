#!/usr/bin/env python3
''' parse table with tracks for jbrowse2 and create config.json file '''
import hashlib

import shutil

import json

import subprocess

import argparse

import os
import sys
import csv
from shlex import quote


def bedpe2bed(bedpe_file, bed_file):
    # new bed wil be in format chrom [(S1+E1)/2]   [(S2+E2)/2]
    with open(bedpe_file, "r") as fin, open(bed_file, "w") as fout:
        for line in fin:
            items = line.strip().split("\t")
            if items[0] != items[3]:
                print("chromosomes in bedpe line are not the same")
                print("skipping line")
                continue
            coord1 = int(int(items[1]) + int(items[2]))/2
            coord2 = int(int(items[4]) + int(items[5]))/2
            if coord1 > coord2:
                coord1, coord2 = coord2, coord1
            fout.write(F"{items[0]}\t{coord1}\t{coord2}\n")




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
            if record['format'] == 'bedpe':
                # this is special case, bedpe is not correctly implemented in jbrowse
                record['new_filename'] = os.path.join(
                        directory, F"{assembly_name}_{record['label']}.bed"
                        )
            else:
                record['new_filename'] = os.path.join(
                        directory, F"{assembly_name}_{record['label']}.{record['format']}"
                        )
            print(record)
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
    # create chromsize file from indexfile
    chromsize_file = F"{ref_path}.chromsizes"
    if not os.path.exists(chromsize_file):
        with open(indexfile, 'r') as fin, open(chromsize_file, 'w') as fout:
            for line in fin:
                items = line.strip().split("\t")
                fout.write(F"{items[0]}\t{items[1]}\n")
    print('')



def get_config_string(track):
    """ use additional information from track to
    to create config string to set color, height, etc
    setting is specific base on the track format"""
    config = ''
    color = track['color'] if track['color'] else 'goldenrod'
    if "displayMode" in track:
        displayMode = track['displayMode'] if track['displayMode'] else 'normal'
    else:
        displayMode = 'normal'

    md5sum = hashlib.md5(track['new_filename'].encode('utf-8')).hexdigest() + "_d"
    displayId = F"{track['label']}_{md5sum}"
    showLabels = True
    if "showLabels" in track:
        if track['showLabels'].lower() == 'false':
            showLabels = False

    if track['format'] == 'bigwig':
        config = {"displays": [{"displayId": displayId, "autoscale": "global", "type":
            "LinearWiggleDisplay",
                                "renderers": {"XYPlotRenderer": {"color": color,
                                    "type": "XYPlotRenderer"}}}]}
    if track['format'] == 'bed' and track['type'] != 'arcs':
        config = {"displays": [{"displayId": displayId, "type": "LinearBasicDisplay",
                                "renderer": {"color1": color,
                                             "displayMode": displayMode,
                                             "type": "SvgFeatureRenderer"}}]}
    if track['format'] == 'gff3':
        config = {"displays": [{"displayId": displayId, "type": "LinearBasicDisplay",
                                "renderer": {"color1": color,
                                             "type": "SvgFeatureRenderer",
                                             "displayMode": displayMode,
                                             "showLabels": showLabels}}]}
    # arcs do not have displayMode!
    if track['format'] == 'bed' and track['type'] == 'arcs':
        config = {"displays": [{"displayId": displayId, "type": "LinearArcDisplay",
                                "renderer": {"color": color,
                                             "type": "ArcRenderer",
                                             "height": "jexl:log10(get(feature,"
                                                       "'end')-get(feature,'start'))*20",
                                             "label": ""}}]}
    if track['format'] == 'bedpe' and track['type'] == 'arcs':
        config = {"displays": [{"displayId": displayId, "type": "LinearArcDisplay",
                                "renderer": {"color": color,
                                             "type": "ArcRenderer",
                                             "label": ""}}]}


    if config!= '':
        # escape quotes in config string
        config_str = json.dumps(config).replace('"', '\\"')
        print('config string:', config_str)
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
    md5sum = hashlib.md5(track['new_filename'].encode('utf-8')).hexdigest()
    trackId = F"{track['label']}_{md5sum}"
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(F"{track['label']}_{track['new_filename']}")
    print('trackId:', trackId)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
    if os.path.exists(indexfile):
        # use index:
        cmd = (F"jbrowse add-track {quote(track['new_filename'])} --out {outdir} "
               F"--load inPlace --trackId {quote(trackId)} --indexFile"
               F" {quote(indexfile)} "
               F" --category {quote(track['category'])}"
               F" --name {quote(track['label'])}"
               F" {cfg_cmd} --assemblyNames {quote(assembly_name)} --overwrite")
        print(cmd)
    else:
        cmd = (F"jbrowse add-track {quote(track['new_filename'])} --out {outdir} "
               F"--load inPlace --trackId {quote(trackId)} "
               F" --category {quote(track['category'])} --overwrite"
               F" --name {quote(track['label'])}"
               F" {cfg_cmd} --assemblyNames {quote(assembly_name)}")
        print(cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        print(F"Error adding track {track['src_filename']}, ")
        print("\n\n")

    print('')

def guess_file_format(filepath):
    with open(filepath, 'rb') as f:
        # Read the first few bytes to check for the BigWig magic number
        header = f.read(4)
        if header in [b'\x88\x8F\xFC\x26', b'\x26\xFC\x8F\x88']:
            return "BigWig"
    # Check for BedGraph format
    with open(filepath, 'r', errors='ignore') as f:
        lines_checked = 0
        for line in f:
            columns = line.strip().split()
            if len(columns) < 4:
                continue
            try:
                int(columns[1])  # Test if the 2nd column is an integer
                int(columns[2])  # Test if the 3rd column is an integer
                float(columns[3])  # Test if the 4th column is a float
                return "BedGraph"
            except ValueError:
                # If casting fails, move on to the next line
                pass
            lines_checked += 1
            if lines_checked > 10:  # You can adjust the number of lines you want to check
                break
    return "Unknown"



def sort_and_index(track, assembly_name, force=False):
    """
    sort and index track, handles bed,  gff and bam files, bigwig files are already sorted
    use always CSI index for every file
    """
    # check if sorted.gz and csi index exists, if not create them

    sortfilegz = F"{track['new_filename']}.sorted.{track['format']}.gz"
    sortfile = F"{track['new_filename']}.sorted.{track['format']}"
    indexfile = F"{track['new_filename']}.sorted.{track['format']}.csi"
    indexfilegz = F"{track['new_filename']}.sorted.{track['format']}.gz.csi"
    # BED

    if track['format'] == 'bed':
        if os.path.exists(sortfilegz) and os.path.exists(indexfilegz) and not force:
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
        if os.path.exists(sortfilegz) and os.path.exists(indexfilegz) and not force:
            print(F"sorted and indexed file {sortfilegz} exists, skipping...")
        else:
            cmd = F"sort -k1,1 -k4,4n {quote(track['new_filename'])} > {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"bgzip {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"tabix -C -p gff {quote(sortfilegz)}"
            subprocess.check_call(cmd, shell=True)
        return sortfilegz
    # BEDPE
    if track['format'] == 'bedpe':
        # this should not be used, only if jbrowse implementation for bedpe is fixed
        print("WARNING: BEDPE tracks are not supported in jbrowse, use bed + arcs "
              "instead")
        if os.path.exists(sortfilegz) and os.path.exists(indexfilegz) and not force:
            print(F"sorted and indexed file {sortfilegz} exists, skipping...")
        else:
            cmd = F"sort -k1,1 -k2,2n -k4,4n {quote(track['new_filename'])} > {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"bgzip {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"tabix -C -p bed {quote(sortfilegz)}"
            subprocess.check_call(cmd, shell=True)
        return sortfilegz
    # BAM
    if track['format'] == 'bam':
        if os.path.exists(sortfile) and os.path.exists(indexfile) and not force:
            print(F"sorted and indexed file {sortfile} exists, skipping...")
        else:
            cmd = F"samtools sort {quote(track['new_filename'])} -o {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
            cmd = F"samtools index -c {quote(sortfile)}"
            subprocess.check_call(cmd, shell=True)
        return sortfile
    if track['format'] == 'bigwig':
        # check if it is really bigwig and not bedgraph
        if guess_file_format(track['new_filename']) == 'BigWig':
            print(F"file {track['new_filename']} - format OK")
        else:
            print(F"file {track['new_filename']} is not bigwig but bedgraph, converting "
                  F"to bigwig...")
            chrsizes =  F"{assembly_name}_{assembly_name}.fasta.chromsizes"
            script_path = os.path.dirname(os.path.realpath(__file__))
            bgr_file = F"{track['new_filename']}"
            bgr_file_sorted = F"{track['new_filename']}.sorted"
            bw_file = F"{track['new_filename']}.bigwig"
            # sort bedgraph
            cmd = F"sort -k1,1 -k2,2n {quote(bgr_file)} | cut -f 1-4 >" \
                  F" {quote(bgr_file_sorted)}"
            subprocess.check_call(cmd, shell=True)
            # verify that all intervals are within the chromosome sizes, if not remove them
            # for example epyc2 sometimes reports intervals that are outside of the
            # chromosome!
            def filter_intervals_by_chromsizes(bedgraph_file, chromsizes_file):
                with open(chromsizes_file, 'r') as f:
                    chromsizes = {}
                    for line in f:
                        chrom, size = line.strip().split()
                        chromsizes[chrom] = int(size)
                with open(bedgraph_file, 'r') as f, open(F"{bedgraph_file}.tmp", 'w') as fout:
                    for line in f:
                        chrom, start, end, score = line.strip().split()
                        if chrom in chromsizes and int(end) <= chromsizes[chrom]:
                             fout.write(line)
                        else:
                            print(F"removing interval {chrom}:{start}-{end}")
                shutil.move(F"{bedgraph_file}.tmp", bedgraph_file)

            filter_intervals_by_chromsizes(bgr_file_sorted, chrsizes)

            # replace bgr_file with bw_file
            cmd = (F"{script_path}/bedGraphToBigWig {quote(bgr_file_sorted)} "
                   F"{quote(chrsizes)} {quote(bw_file)}")
            subprocess.check_call(cmd, shell=True)
            shutil.move(bw_file, bgr_file)
        print('-------x-x-x-x-x bigwig_file:')
        print(track['new_filename'])
        return track['new_filename']
    print('format not recognized, skipping sorting and indexing')
    print(track['format'])
    print('-------------------')
    return track['new_filename']


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--output_dir', type=str, default="data")
    parser.add_argument("-c", "--config_table", required=True)
    parser.add_argument("-u", "--update", action="store_true", default=False)
    parser.add_argument("-f", "--force", action="store_true", default=False,
                        help="force overwrite of existing files")
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
        print("track:", track['label'])
        print("format:", track['format'])
        print("//////////////////////////////////////////")
        if not os.path.exists(track['new_filename']) or args.force:
            if track['format'] == 'bedpe':
                # this is not implemented correctly in jbrowse
                # convert bedpe to bed.
                bedpe2bed(track['src_filename'], track['new_filename'])
                track['format'] = 'bed'
                track['type'] = 'arcs'
                print('converting bedpe to bed:', track['src_filename'])
                print("========================================")
            else:
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
        track['new_filename'] = sort_and_index(track, assembly_name, args.force)
        print('adding track:', track['label'])
        add_track(track, args.output_dir, assembly_name)
        print("----------------------------------------")
        print("----------------------------------------")

    # add configuration to json file
    json_file = os.path.join(args.output_dir, 'config.json')
    # load json file
    with open(json_file, 'r') as f:
        data = json.load(f)
    # add configuration to json file
    # "configuration": {"logoPath": {"uri": "./elixir_150x48.svg"}},
    data["configuration"] = {"logoPath": {"uri": "./elixir_150x48.svg"}}
    # write to json file
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)

