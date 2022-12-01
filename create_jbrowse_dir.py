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
    print('adding link to ', src)
    print('creating symlink to ', dest)
    print(os.path.exists(dest))
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
        try:
            os.symlink(src, dest)
        except FileExistsError:
            if os.path.islink(dest):
                print(dest, " already exists")
                print("replacing link to", os.readlink(dest), ' with ', src)
                os.unlink(dest)
                os.symlink(src, dest)
            else:
                print(dest, " already exists and it is not link")
                sys.exit()


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

def bed_to_tabix(record):
    """sort, compress and create tabix index"""
    bed = record['new_filename']
    output_file = F"{bed}.sorted.bed.gz"
    if os.path.exists(output_file):
        print(output_file, "already exists, skipping..")
    else:
        makesymlink(record['src_filename'], bed)
        os.system(F"sort -k1,1 -k 2,2n {bed} > {bed}.sorted.bed")
        # remove lines wich started with "track name=" and "browser position"
        os.system(F"grep -v '^track name=' {bed}.sorted.bed | grep -v '^browser position' > {bed}.sorted.bed.tmp")
        os.system(F"mv {bed}.sorted.bed.tmp {bed}.sorted.bed")
        os.system(F"bgzip -f {bed}.sorted.bed")
        # make csi index
        try:
            subprocess.check_call(['tabix', '-p', 'bed', output_file])
            print('tabix index created')
        except subprocess.CalledProcessError:
            print("tabix tbi failed, trying with csi index")
            subprocess.check_call(['tabix', '-C', '-p', 'bed', output_file])
    return os.path.basename(output_file)


def gff_to_tabix(record):
    """sort, compress and create tabix index"""
    gff = record['new_filename']
    print("gff", gff)
    output_file = F"{gff}.sorted.gff.gz"
    output_bigwig = F"{gff}.sorted.gff.gz.bw"
    print("output_file", output_file)
    print("------------------------")
    # check is corresponding bigwig file exist - same name but with bw extension
    print("checking if bigwig file exist")
    print(record['src_filename']+".bw")
    if os.path.exists(record['src_filename']+".bw"):
        print('bigwig file exist, adding it to gff record')
        # make symlink to bigwig file, owerride if exist
        makesymlink(record['src_filename'] + '.bw', output_bigwig)

    if os.path.exists(output_file):
        print(output_file, "already exists, skipping..")
    else:
        makesymlink(record['src_filename'], gff)
        os.system(F"sort -k1,1 -k 4,4n {gff} > {gff}.sorted.gff")
        os.system(F"bgzip -f {gff}.sorted.gff")
        # make csi index
        try:
            subprocess.check_call(['tabix', '-p', 'gff', output_file])
        except subprocess.CalledProcessError:
            print("tabix tbi failed, trying with csi index")
            subprocess.check_call(['tabix', '-C', '-p', 'gff', output_file])
    return os.path.basename(output_file)


def add_bam(record, f):
    """ create bam index if does not exists and return bam track configuration"""
    print('adding bam file')
    makesymlink(record['src_filename'], record['new_filename'])
    if not os.path.exists(record['src_filename'] + ".csi"):
        # make index on source location
        subprocess.check_call(['samtools', 'index','-c', record['src_filename']])
    makesymlink(record['src_filename'] + ".csi", record['new_filename'] + ".csi")

    bam = os.path.basename(record['new_filename'])
    conf_str = "\n".join(
        [F"[tracks.{record['label']}]", F"urlTemplate={bam}",
         "storeClass=JBrowse/Store/SeqFeature/BAM", F"type={record['type']}",
         F"csiUrlTemplate={bam}.csi",
         F"label={record['label']}", F"category={record['category']}",
         "\n"]
        )
    f.write(conf_str)

def add_general_track_settings2(record):
    """ add general track settings"""
    conf_str = "\n".join(
        [F"chunkSizeLimit={record['chunkSizeLimit']}",
         F"style.color={record['color']}",
         F"histograms.color={record['color']}",
         F"histograms.height={record['histograms.height']}",
         F"displayMode={record['displayMode']}",
         "\n"]
        )
    return conf_str


def add_general_track_settings(record):
    """check for additional settings available in record
    """
    excluded_keys = ['type', 'src_filename', 'new_filename', 'label', 'category',
                     'color', 'new_src_filename', 'dirname', 'filename', 'format']
    conf_str = ""
    for key, value in record.items():
        if key not in excluded_keys:
            if value != '':
                conf_str += F"{key}={value}\n"
    conf_str += "\n"
    return conf_str


def add_gff3(record, f):
    """ convert gff3 to tabix format and return gff3 track configuration"""
    gffformated = gff_to_tabix(record)
    gffformated_full_path = os.path.join(args.output_dir, gffformated)
    # check if tbi index is available
    print("checking if tbi index is available")
    print(gffformated + ".tbi")
    if os.path.exists(gffformated_full_path + ".tbi"):
        print("tbi index found")
        conf_str = "\n".join(
            [F"[tracks.{record['label']}]", F"urlTemplate={gffformated}",
             "storeClass=JBrowse/Store/SeqFeature/GFF3Tabix", F"type={record['type']}",
             F"label={record['label']}", F"category={record['category']}\n"
             ]
            )
    else:
        print("csi index found")
        conf_str = "\n".join(
            [F"[tracks.{record['label']}]", F"urlTemplate={gffformated}",
             F"csiUrlTemplate={gffformated}.csi",
             "storeClass=JBrowse/Store/SeqFeature/GFF3Tabix", F"type={record['type']}",
             F"label={record['label']}", F"category={record['category']}\n"
             ]
            )

    # add link to bigwig file if exist
    print("--checking if bigwig file exist--")
    print(gffformated_full_path + ".bw")
    if os.path.exists(gffformated_full_path + ".bw"):
        print('bigwig file exist, adding to track configuration')
        conf_str += "histograms.storeClass=JBrowse/Store/SeqFeature/BigWig\n"
        conf_str += F"histograms.urlTemplate={gffformated}.bw\n"
    conf_str += add_general_track_settings(record)
    
    f.write(conf_str)

def add_bed(record, f):
    """ convert bed to tabix format and return bed track configuration"""
    bedformated = bed_to_tabix(record)
    bedformated_full_path = os.path.join(args.output_dir, bedformated)
    # check if tbi index is available
    print("checking if tbi index is available: " + bedformated + ".tbi")
    if os.path.exists(bedformated_full_path + ".tbi"):
        print("tbi index found")
        conf_str = "\n".join(
            [F"[tracks.{record['label']}]", F"urlTemplate={bedformated}",
             "storeClass=JBrowse/Store/SeqFeature/BEDTabix", F"type={record['type']}",
             F"label={record['label']}", F"category={record['category']}",
             "\n"]
            )
    else:
        print("tbi index not found, trying with csi index")
        conf_str = "\n".join(
            [F"[tracks.{record['label']}]", F"urlTemplate={bedformated}",
             F"csiUrlTemplate={bedformated}.csi",
             "storeClass=JBrowse/Store/SeqFeature/BEDTabix", F"type={record['type']}",
             F"label={record['label']}", F"category={record['category']}",
             "\n"]
            )
    conf_str += add_general_track_settings(record)
    f.write(conf_str)



def add_bed12(record, f):
    """ first convert bed12 to gff3 and returs track configuration"""
    record['new_src_filename'] = record['new_filename'] + ".gff3"


    print("new_src_filename", record['new_src_filename'])
    print("new_file", record['new_filename'])
    bed12_to_gff3(record['src_filename'], record['new_src_filename'])
    # remove leading directory name from new_filename
    record['new_src_filename'] = os.path.basename(record['new_src_filename'])
    record['src_filename'] = record['new_src_filename']
    print('src_filename', record['src_filename'], "updated")
    print("=====================================")
    add_gff3(record, f)


def bed12_to_gff3(bedfile, gff3out):
    """convert bed12 to gff3"""
    with open(bedfile) as f, open(gff3out, 'w') as gff3:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            line = line.split("\t")
            # check if bed12
            if len(line) == 12:
                chrom = line[0]
                start = line[1]
                end = line[2]
                name = line[3]
                score = line[4]
                strand = line[5]
                thickStart = line[6]
                thickEnd = line[7]
                itemRgb = line[8]
                blockCount = line[9]
                blockSizes = line[10]
                blockStarts = line[11]
                # get block sizes
                blockSizes = blockSizes.split(",")
                blockSizes = [int(x) for x in blockSizes if x != ""]
                # get block starts
                blockStarts = blockStarts.split(",")
                blockStarts = [int(x) for x in blockStarts if x != ""]
                # get exon number
                exon_number = len(blockSizes)
                # get exon start and end
                exon_starts = []
                exon_ends = []
                for i in range(exon_number):
                    exon_starts.append(int(start) + blockStarts[i])
                    exon_ends.append(int(start) + blockStarts[i] + blockSizes[i])
                # write gff3

                for i in range(exon_number):
                    gff3.write(F"{chrom}\tbed12\tbiological_region\t{exon_starts[i]}\t"
                               F"{exon_ends[i]}\t{score}\t{strand}\t.\t"
                               F"Parent={name}\n")
                gff3.write(F"{chrom}\tbed12\tsequence_feature\t{start}\t{end}\t{score}\t"
                           F"{strand}\t.\tID={name};Name={name}\n")
            else:
                print("line not bed12 format, skipping..")


# def add_track(record, f):
def add_bigwig(record, f):
    """ takes bigwig file and return bigwig track configuration"""
    print('adding bigwig file')
    makesymlink(record['src_filename'], record['new_filename'])
    bigwig = os.path.basename(record['new_filename'])
    conf_str = "\n".join(
        [F"[tracks.{record['label']}]", F"urlTemplate={bigwig}",
         "storeClass=JBrowse/Store/SeqFeature/BigWig", F"type={record['type']}",
         F"label={record['label']}", F"category={record['category']}",
         "\n"]
        )
    conf_str += add_general_track_settings(record)
    f.write(conf_str)


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

    add_track = {'gff3': add_gff3, 'bam': add_bam, 'bigwig': add_bigwig,
                 'bed': add_bed, 'bed12': add_bed12}

    # add tracks
    with open(os.path.join(args.output_dir, "tracks.conf"), 'w') as cfile:
        for i in range(1, len(track_list)):
            add_track[track_list[i]['format']](track_list[i], cfile)

    # TODO - add force option to overwrite everything  # TODO a add option to use csi
    #  (csiUrlTemplate) for large scaffolds
