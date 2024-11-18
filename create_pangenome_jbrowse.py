#!/usr/bin/env python3

import os
import pandas as pd
import subprocess

# Constants
GENOMES_CSV = 'genomes.csv'       # Contains headers: genome, abs_path
REPEATS_CSV = 'repeats.csv'       # Contains headers: genome, abs_path
CHIP_SEQ_CSV = 'chip_seq.csv'     # Contains headers: genome, abs_path
OLIGOS_CSV = 'oligos.csv'         # Contains headers: genome, abs_path
TRACK_CONFIG_CSV = 'track_config.csv'  # Contains headers including data_source, relative_path, and label
CONFIG_CSV_DIR = './config_csv_files'
JBROWSE_CONFIG_SCRIPT = '/mnt/raid/users/petr/workspace/jbrowse_prepare/create_jbrowse2_config.py'
MAKE_PAF_SCRIPT = '/mnt/raid/users/petr/workspace/granges_tools/make_paf_from_probes.R'
ADD_SYNTENY_SCRIPT = '/mnt/raid/users/petr/workspace/jbrowse_prepare/add_synteny_tracks_from_paf.py'
# run it from where data will be stored
CURRENT_DIR = "."

# Ensure the output directory exists
if not os.path.exists(CONFIG_CSV_DIR):
    os.makedirs(CONFIG_CSV_DIR)

# Step 1: Read genomes.csv
genomes_df = pd.read_csv(GENOMES_CSV, sep='\t')
print(genomes_df)
genomes_df['genome'] = genomes_df['genome'].str.strip()
genomes_df['abs_path'] = genomes_df['abs_path'].str.strip()
genome_list = genomes_df['genome'].tolist()

# Step 2: Read track_config.csv
track_config_df = pd.read_csv(TRACK_CONFIG_CSV, sep='\t')
track_config_df['relative_path'] = track_config_df['relative_path'].str.strip()
track_config_df['data_source'] = track_config_df['data_source'].str.strip()
track_config_df['label'] = track_config_df['label'].str.strip()

# Separate track configurations by data_source
repeats_track_config_df = track_config_df[track_config_df['data_source'] == 'repeats.csv']
chip_seq_track_config_df = track_config_df[track_config_df['data_source'] == 'chip_seq.csv']

# Step 3: Read repeats.csv
if os.path.exists(REPEATS_CSV):
    repeats_df = pd.read_csv(REPEATS_CSV, sep='\t')
    repeats_df['genome'] = repeats_df['genome'].str.strip()
    repeats_df['abs_path'] = repeats_df['abs_path'].str.strip()
else:
    repeats_df = pd.DataFrame(columns=['genome', 'abs_path'])

# Step 4: Read chip_seq.csv
if os.path.exists(CHIP_SEQ_CSV):
    chip_seq_df = pd.read_csv(CHIP_SEQ_CSV, sep='\t')
    chip_seq_df['genome'] = chip_seq_df['genome'].str.strip()
    chip_seq_df['abs_path'] = chip_seq_df['abs_path'].str.strip()
else:
    chip_seq_df = pd.DataFrame(columns=['genome', 'abs_path'])

# Configuration CSV columns
config_columns = [
    'label', 'format', 'category', 'type', 'dirname', 'filename',
    'style.color', 'color', 'displayMode', 'showLabels'
]

# Initialize a dictionary to store configuration DataFrames for each genome
config_dfs = {}

# Step 5: For each genome in genomes.csv, create the configuration DataFrame
for idx, row in genomes_df.iterrows():
    genome_name = row['genome']
    genome_abs_path = row['abs_path']

    # Initialize a DataFrame for the genome's configuration
    config_df = pd.DataFrame(columns=config_columns)

    # Add the reference track (genome fasta)
    if not os.path.exists(genome_abs_path):
        print(f"Warning: genome.fasta not found for genome {genome_name} in {genome_abs_path}")
        continue

    reference_track = {
        'label': genome_name,
        'format': 'fasta',
        'category': 'reference',
        'type': 'reference',
        'dirname': os.path.dirname(genome_abs_path),
        'filename': os.path.basename(genome_abs_path),
        'style.color': '',
        'color': '',
        'displayMode': '',
        'showLabels': ''
    }
    config_df = pd.concat([config_df, pd.DataFrame([reference_track])], ignore_index=True)


    # Add repeats tracks if available
    repeats_genome_row = repeats_df[repeats_df['genome'] == genome_name]
    if not repeats_genome_row.empty:
        repeats_abs_path = repeats_genome_row.iloc[0]['abs_path']
        genome_repeats_config = repeats_track_config_df.copy()
        for _, track_row in genome_repeats_config.iterrows():
            relative_path = track_row['relative_path']
            full_repeats_path = os.path.join(repeats_abs_path, relative_path)
            if os.path.exists(full_repeats_path):
                # Get configuration values
                label = track_row['label']
                format = track_row['format']
                category = track_row['category']
                style_color = track_row['style.color']
                color = track_row['color']
                displayMode = track_row['displayMode']
                showLabels = track_row['showLabels']

                # Define the type based on format
                if format in ['gff3', 'bed']:
                    track_type = 'CanvasFeature'
                elif format == 'bigwig':
                    track_type = 'JBrowse/View/Track/Wiggle/XYPlot'
                else:
                    track_type = ''

                repeats_track = {
                    'label': label,
                    'format': format,
                    'category': category,
                    'type': track_type,
                    'dirname': repeats_abs_path,
                    'filename': relative_path,
                    'style.color': style_color,
                    'color': color,
                    'displayMode': displayMode,
                    'showLabels': showLabels
                }
                config_df = pd.concat([config_df, pd.DataFrame([repeats_track])], ignore_index=True)

            else:
                print(f"Warning: {relative_path} not found for genome {genome_name} in {repeats_abs_path}")

    # Add chip-seq tracks if available
    chipseq_genome_row = chip_seq_df[chip_seq_df['genome'] == genome_name]
    if not chipseq_genome_row.empty:
        chipseq_abs_path = chipseq_genome_row.iloc[0]['abs_path']
        genome_chipseq_config = chip_seq_track_config_df.copy()
        for _, track_row in genome_chipseq_config.iterrows():
            relative_path = track_row['relative_path']
            full_chipseq_path = os.path.join(chipseq_abs_path, relative_path)
            if os.path.exists(full_chipseq_path):
                # Get configuration values
                label = track_row['label']
                format = track_row['format']
                category = track_row['category']
                style_color = track_row['style.color']
                color = track_row['color']
                displayMode = track_row['displayMode']
                showLabels = track_row['showLabels']

                # Define the type based on format
                if format == 'bigwig':
                    track_type = 'JBrowse/View/Track/Wiggle/XYPlot'
                else:
                    track_type = ''

                chip_seq_track = {
                    'label': label,
                    'format': format,
                    'category': category,
                    'type': track_type,
                    'dirname': chipseq_abs_path,
                    'filename': relative_path,
                    'style.color': style_color,
                    'color': color,
                    'displayMode': displayMode,
                    'showLabels': showLabels
                }
                config_df = pd.concat([config_df, pd.DataFrame([chip_seq_track])], ignore_index=True)
            else:
                print(f"Warning: {relative_path} not found for genome {genome_name} in {chipseq_abs_path}")

    # Save the configuration DataFrame for this genome
    config_dfs[genome_name] = config_df

# Step 6: Write the configuration CSV files and run create_jbrowse2_config.py
first_genome = True
for genome_name, config_df in config_dfs.items():
    config_csv_file = os.path.join(CONFIG_CSV_DIR, f'genome_{genome_name}.csv')
    config_df.to_csv(config_csv_file, sep='\t', index=False)

    # Run create_jbrowse2_config.py
    cmd = [
        JBROWSE_CONFIG_SCRIPT,
        '-d', CURRENT_DIR,
        '-c', config_csv_file,
    ]
    if not first_genome:
        # For genomes after the first one, use the -u option to update
        cmd.append('-u')
    try:
        subprocess.run(cmd, check=True)
        print(f"Created jbrowse configuration for genome {genome_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running create_jbrowse2_config.py for genome {genome_name}: {e}")
    first_genome = False

# Step 7: Process synteny tracks

# Read oligos.csv
if os.path.exists(OLIGOS_CSV):
    oligos_df = pd.read_csv(OLIGOS_CSV, sep='\t')
    oligos_df['genome'] = oligos_df['genome'].str.strip()
    oligos_df['abs_path'] = oligos_df['abs_path'].str.strip()
else:
    oligos_df = pd.DataFrame(columns=['genome', 'abs_path'])


# Create symbolic links to oligo.bed files
genome_oligo_files = {}
for idx, row in oligos_df.iterrows():
    genome_name = row['genome']
    abs_path = row['abs_path']  # Direct path to the oligo BED file
    if os.path.exists(abs_path):
        # Create symbolic link in the current directory
        symlink_name = f"{genome_name}_oligos.bed"
        symlink_path = os.path.join(CURRENT_DIR, symlink_name)
        if not os.path.exists(symlink_path):
            os.symlink(abs_path, symlink_name)
            print(f"Created symlink {symlink_name} -> {abs_path}")
        else:
            print(f"Symlink {symlink_name} already exists")
        genome_oligo_files[genome_name] = symlink_name
    else:
        print(f"Oligo BED file {abs_path} does not exist for genome {genome_name}")

# Step 8: For all pairs of genomes, generate PAF files
# Create synteny_paf_tracks.csv
synteny_paf_tracks_csv = 'synteny_paf_tracks.csv'
with open(synteny_paf_tracks_csv, 'w') as paf_csv:
    paf_csv.write("GENOME1\tGENOME2\tPAF\n")
    genomes = list(genome_oligo_files.keys())
    length = len(genomes)
    for i in range(length):
        for j in range(i+1, length):
            genome1 = genomes[i]
            genome2 = genomes[j]
            g1_oligo = genome_oligo_files[genome1]
            g2_oligo = genome_oligo_files[genome2]
            g1_oligo_sample = f"{genome1}_oligos_sample.bed"
            g2_oligo_sample = f"{genome2}_oligos_sample.bed"

            # Create sample files using awk command
            # Here we will use subprocess to run the awk commands
            seed = '123'
            if not os.path.exists(g1_oligo_sample):
                awk_cmd_g1 = f"awk -v seed=\"{seed}\" 'BEGIN {{srand(seed)}} rand() < 0.20' {g1_oligo} > {g1_oligo_sample}"
                subprocess.run(awk_cmd_g1, shell=True)
            if not os.path.exists(g2_oligo_sample):
                awk_cmd_g2 = f"awk -v seed=\"{seed}\" 'BEGIN {{srand(seed)}} rand() < 0.20' {g2_oligo} > {g2_oligo_sample}"
                subprocess.run(awk_cmd_g2, shell=True)

            # Get the fasta index files
            g1_fai = f"{genome1}_{genome1}.fasta.fai"
            g2_fai = f"{genome2}_{genome2}.fasta.fai"

            # The fai files should be generated from the fasta files
            # Use genomes_df to get the fasta paths
            genome1_row = genomes_df[genomes_df['genome'] == genome1]
            genome2_row = genomes_df[genomes_df['genome'] == genome2]

            if not genome1_row.empty and not genome2_row.empty:
                genome1_fasta_path = os.path.join(genome1_row['abs_path'].values[0], 'genome.fasta')
                genome2_fasta_path = os.path.join(genome2_row['abs_path'].values[0], 'genome.fasta')

                # Create the fai files if they don't exist
                if not os.path.exists(g1_fai):
                    cmd_index_g1 = ['samtools', 'faidx', genome1_fasta_path]
                    subprocess.run(cmd_index_g1)
                    # Copy the fai file to the current directory with the appropriate name
                    os.rename(f"{genome1_fasta_path}.fai", g1_fai)
                if not os.path.exists(g2_fai):
                    cmd_index_g2 = ['samtools', 'faidx', genome2_fasta_path]
                    subprocess.run(cmd_index_g2)
                    os.rename(f"{genome2_fasta_path}.fai", g2_fai)
            else:
                print(f"Genome fasta paths not found for {genome1} or {genome2}")
                continue

            # Define output PAF files
            outpaf = f"{genome1}_{genome2}.oligo.paf"
            outpaf_sample = f"{genome1}_{genome2}.oligo_sample.paf"

            # Write to synteny_paf_tracks.csv
            paf_csv.write(f"{genome1}\t{genome2}\t{outpaf}\n")
            paf_csv.write(f"{genome1}\t{genome2}\t{outpaf_sample}\n")

            # Run make_paf_from_probes.R
            if not os.path.exists(outpaf):
                cmd_paf = [
                    'Rscript', MAKE_PAF_SCRIPT,
                    '-a', g1_oligo,
                    '-b', g2_oligo,
                    '-A', g1_fai,
                    '-B', g2_fai,
                    '-o', outpaf
                ]
                subprocess.run(cmd_paf)
            if not os.path.exists(outpaf_sample):
                cmd_paf_sample = [
                    'Rscript', MAKE_PAF_SCRIPT,
                    '-a', g1_oligo_sample,
                    '-b', g2_oligo_sample,
                    '-A', g1_fai,
                    '-B', g2_fai,
                    '-o', outpaf_sample
                ]
                subprocess.run(cmd_paf_sample)

# Step 9: Run add_synteny_tracks_from_paf.py
config_json = 'config.json'
output_config_json = 'config_with_synteny.json'
cmd_add_synteny = [
    ADD_SYNTENY_SCRIPT,
    '-c', config_json,
    '-t', synteny_paf_tracks_csv,
    '-o', output_config_json
]
try:
    subprocess.run(cmd_add_synteny, check=True)
    print(f"Generated final configuration with synteny tracks in {output_config_json}")
except subprocess.CalledProcessError as e:
    print(f"Error running add_synteny_tracks_from_paf.py: {e}")
