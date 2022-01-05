#!/usr/bin/env python3
"""
this script is using csv table as input to create data fro jbrowse, see tamplate for
csv formatting example
"""
import argparse
import csv
import os
import sys
import subprocess
import glob


# environment created using
#
# mamba  create -n apollo_utils -c conda-forge -c bioconda blat blast
# bioconductor-rtracklayer  jbrowse "python>=3.8" ucsc-fatotwobit
#


def read_config(config_file, directory):
    """
    read csv table, return list of tracks
    """
    input_tracks = []
    print('reading config file:')
    with open(config_file, 'r') as f:
        reader = csv.DictReader(f, delimiter="\t")
        labels = set()
        for record in reader:
            if record['label'] in labels:
                print('duplicated label:', record['label'])
                sys.exit('labels must be unique, exiting...')
            labels.add(record['label'])
            record['new_filename'] = os.path.join(
                directory, F"{record['label']}.{record['format']}"
                )
            record['src_filename'] = os.path.join(
                record['dirname'], record['filename']
                )
            input_tracks.append(record)
    print('')
    return input_tracks


def makesymlink(src, dest):
    """ creates or update symlinks for file, if file exist it raises expetion
    """
    if os.path.exists(dest):
        if os.path.islink(dest):
            if os.readlink(dest) == src:
                print(dest, " already point to ", src)
            else:
                print("replacing link to", os.readlink(dest), ' with ', src)
                os.unlink(dest)
                os.symlink(src, dest)
        else:
            print(dest, " already exists and it is not link")
            sys.exit()
    else:
        os.symlink(src, dest)


def create_reference(fasta, directory):
    """ format fasta file wo 2bit and fai and create blast databases
    """
    if not os.path.exists(directory):
        os.mkdir(directory)
    # make symlinks to data dir
    reference = os.path.join(directory, "reference.fasta")
    if os.path.exists(reference):
        print('reference file  exists')
    else:
        print('adding link to reference ', )
        os.symlink(fasta, reference)
    if os.path.exists(F"{reference}.2bit"):
        print('2bit reference already exist')
    else:
        subprocess.check_call(["faToTwoBit", reference, F"{reference}.2bit"])
    if glob.glob(F"{reference}.2bit.*nsq"):
        print('blast database already exists')
    else:
        print("creating database")
        subprocess.check_call(
            ["makeblastdb", '-in', reference, '-out', F"{reference}.2bit", "-dbtype",
             "nucl"]
            )
    if os.path.exists(F"{reference}.fai"):
        print("faidx already exist")
    else:
        print("creating faidx database")
        subprocess.check_call(['samtools', 'faidx', reference])

    subprocess.check_call(
        ["prepare-refseqs.pl", "--indexed_fasta", reference, "--out", directory]
        )


def gff_to_tabix(record):
    """sort, compress and create tabix index"""
    gff = record['new_filename']
    output_file = F"{gff}.sorted.gff.gz"
    if os.path.exists(output_file):
        print(output_file, "already exists, skipping..")
    else:
        makesymlink(record['src_filename'], gff)
        os.system(F"sort -k1,1 -k 4,4n {gff} > {gff}.sorted.gff")
        os.system(F"bgzip -f {gff}.sorted.gff")
        # make csi index
        os.system(F"tabix -C -p gff {output_file}")
        # make tbi index
        os.system(F"tabix -p gff {output_file}")
    return os.path.basename(output_file)


def add_bam(record, f):
    """ create bam index if does not exists and return bam track configuration"""
    print('adding bam file')
    makesymlink(record['src_filename'], record['new_filename'])
    if not os.path.exists(record['src_filename'] + ".bai"):
        # make index on source location
        subprocess.check_call(['samtools', 'index', record['src_filename']])
    makesymlink(record['src_filename'] + ".bai", record['new_filename'] + ".bai")

    bam = os.path.basename(record['new_filename'])
    conf_str = "\n".join(
        [F"[tracks.{record['label']}]", F"urlTemplate={bam}",
         "storeClass=JBrowse/Store/SeqFeature/BAM", F"type={record['type']}",
         F"label={record['label']}", F"category={record['category']}", "\n"]
        )
    f.write(conf_str)


def add_gff3(record, f):
    """ convert gff3 to tabix format and return gff3 track configuration"""
    gffformated = gff_to_tabix(record)
    conf_str = "\n".join(
        [F"[tracks.{record['label']}]", F"urlTemplate={gffformated}",
         "storeClass=JBrowse/Store/SeqFeature/GFF3Tabix", F"type={record['type']}",
         F"label={record['label']}", F"category={record['category']}", "\n"]
        )
    f.write(conf_str)


# def add_track(record, f):


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--output_dir', type=str, default="data")
    parser.add_argument("-c", "--config_table", required=True)

    args = parser.parse_args()

    track_list = read_config(args.config_table, args.output_dir)

    # first track is always genome reference
    ref_path = os.path.join(
        track_list[0]['dirname'], track_list[0]['filename']
        )
    print('creating reference files:')
    create_reference(ref_path, args.output_dir)

    add_track = {'gff3': add_gff3, 'bam': add_bam}

    # add tracks
    with open(os.path.join(args.output_dir, "tracks.conf"), 'w') as cfile:
        for i in range(1, len(track_list)):
            add_track[track_list[i]['format']](track_list[i], cfile)

    # TODO - add force option to overwrite everything  # TODO a add option to use csi
    #  (csiUrlTemplate) for large scaffolds
