#!/usr/bin/env python3
"""
Create JBrowse2 configuration for multiple genomes from CSV input.

This script takes a CSV table with multiple genomes and creates JBrowse2
directories and configurations. It handles reference genomes, tracks, and
synteny tracks with PAF generation.
"""

import os
import sys
import csv
import json
import argparse
import subprocess
import hashlib
from collections import defaultdict
from typing import Dict, List, Tuple
from shlex import quote

# Import functions from existing scripts
from paf_chain_filter import merge_paf_intervals


def blast_to_bed(blast_file: str, bed_file: str) -> None:
    """
    Convert BLAST tabular output to BED format for synteny analysis.
    
    BLAST format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    BED format: chr start end name score strand
    
    Strand is determined from start/end positions: + if start < end, - otherwise
    """
    print(f"Converting BLAST file {blast_file} to BED format {bed_file}")
    
    with open(blast_file, 'r') as blast_in, open(bed_file, 'w') as bed_out:
        for line in blast_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
                
            qseqid = fields[0]  # This is the marker ID
            sseqid = fields[1]  # This is the chromosome
            sstart = int(fields[8])
            send = int(fields[9])
            
            # Determine strand and coordinates
            if sstart < send:
                start = sstart
                end = send
                strand = '+'
            else:
                start = send
                end = sstart
                strand = '-'
            
            # Convert to 0-based coordinates (BED format)
            start = start - 1
            
            # Write BED line: chr start end name score strand
            bed_out.write(f"{sseqid}\t{start}\t{end}\t{qseqid}\t1\t{strand}\n")


def read_config_csv(csv_file: str) -> Dict[str, List[Dict]]:
    """
    Read CSV configuration file and group tracks by genome.
    
    Returns dictionary with genome names as keys and list of track records as values.
    """
    genomes = defaultdict(list)
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Skip empty rows or rows where Genome column is empty/whitespace
            genome = row.get('Genome', '').strip()
            if not genome:
                continue
            genomes[genome].append(row)
    
    return dict(genomes)


def check_files_exist(genomes: Dict[str, List[Dict]]) -> bool:
    """
    Check if all files specified in the CSV configuration exist.
    
    Returns True if all files exist, False otherwise.
    """
    print("Checking if all files exist...")
    all_files_exist = True
    missing_files = []
    
    for genome_name, tracks in genomes.items():
        print(f"\nChecking files for genome: {genome_name}")
        
        for track in tracks:
            file_path = os.path.join(track['Dirname'], track['Filename'])
            if os.path.exists(file_path):
                print(f"  ✓ {track['Label']}: {file_path}")
            else:
                print(f"  ✗ {track['Label']}: {file_path} (MISSING)")
                missing_files.append(file_path)
                all_files_exist = False
        
        # Check for optional rDNA files
        ref_dirname = None
        for track in tracks:
            if track['type'] == 'reference':
                ref_dirname = track['Dirname']
                break
        
        if ref_dirname:
            rdna_bed_path = os.path.join(ref_dirname, 'rdna.bed')
            if os.path.exists(rdna_bed_path):
                print(f"  ✓ rDNA (optional): {rdna_bed_path}")
            else:
                print(f"  - rDNA (optional): {rdna_bed_path} (not found, will skip)")
    
    if all_files_exist:
        print(f"\n✓ All {sum(len(tracks) for tracks in genomes.values())} required files exist!")
        return True
    else:
        print(f"\n✗ {len(missing_files)} file(s) are missing:")
        for missing_file in missing_files:
            print(f"  - {missing_file}")
        return False


def create_reference(genome_name: str, ref_record: Dict, output_dir: str) -> str:
    """
    Create reference assembly for JBrowse2 using jbrowse add-assembly command.
    
    Returns the path to the reference file in the output directory.
    """
    src_path = os.path.join(ref_record['Dirname'], ref_record['Filename'])
    ref_filename = f"{genome_name}.fasta"
    ref_path = os.path.join(output_dir, ref_filename)
    
    print(f"Creating reference for {genome_name}")
    
    # Check if source file exists
    if not os.path.exists(src_path):
        raise FileNotFoundError(f"Reference file not found: {src_path}")
    
    # Create symlink to reference file
    if not os.path.exists(ref_path):
        os.symlink(src_path, ref_path)
        print(f"Created symlink {ref_path} -> {src_path}")
    else:
        print(f"Skipping symlink creation (already exists): {ref_path}")
    
    # Create FAI index if it doesn't exist
    fai_path = f"{ref_path}.fai"
    if not os.path.exists(fai_path):
        cmd = f"samtools faidx {quote(ref_path)}"
        subprocess.check_call(cmd, shell=True)
        print(f"Created FAI index: {fai_path}")
    else:
        print(f"Skipping FAI index creation (already exists): {fai_path}")
    
    # Add assembly to JBrowse2 (run from output directory to avoid path warnings)
    rel_ref_path = os.path.basename(ref_path)  # Use relative path within output dir
    cmd = f"jbrowse add-assembly {quote(rel_ref_path)} --load inPlace --name {quote(genome_name)} --overwrite"
    print(f"Running: {cmd} (in directory: {output_dir})")
    
    # Change to output directory and run command
    original_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        subprocess.check_call(cmd, shell=True)
    finally:
        os.chdir(original_cwd)
    
    return ref_path


def sort_and_index_track(track_record: Dict, genome_name: str, output_dir: str) -> str:
    """
    Sort and index track files for JBrowse2 compatibility.
    
    Returns the path to the processed track file.
    """
    src_path = os.path.join(track_record['Dirname'], track_record['Filename'])
    track_format = track_record['Format']
    label = track_record['Label']
    
    # Check if source file exists
    if not os.path.exists(src_path):
        raise FileNotFoundError(f"Track file not found: {src_path} (Label: {label})")
    
    # Create track filename
    track_filename = f"{genome_name}_{label}.{track_format}"
    track_path = os.path.join(output_dir, track_filename)
    
    # Create symlink to source file
    if not os.path.exists(track_path):
        os.symlink(src_path, track_path)
        print(f"Created symlink {track_path} -> {src_path}")
    else:
        print(f"Skipping symlink creation (already exists): {track_path}")
    
    if track_format in ['gff3', 'gff']:
        # Sort and compress GFF3 files
        sorted_file = f"{track_path}.sorted.{track_format}.gz"
        index_file = f"{sorted_file}.csi"
        if not os.path.exists(sorted_file) or not os.path.exists(index_file):
            temp_sorted = f"{track_path}.sorted.{track_format}"
            cmd = f"sort -k1,1 -k4,4n {quote(track_path)} > {quote(temp_sorted)}"
            subprocess.check_call(cmd, shell=True)
            cmd = f"bgzip {quote(temp_sorted)}"
            subprocess.check_call(cmd, shell=True)
            cmd = f"tabix -C -p gff {quote(sorted_file)}"
            subprocess.check_call(cmd, shell=True)
            print(f"Sorted and indexed GFF3: {sorted_file}")
        else:
            print(f"Skipping GFF3 sorting/indexing (already exists): {sorted_file}")
        return sorted_file
        
    elif track_format == 'bed':
        # Sort and compress BED files
        sorted_file = f"{track_path}.sorted.bed.gz"
        index_file = f"{sorted_file}.csi"
        if not os.path.exists(sorted_file) or not os.path.exists(index_file):
            temp_sorted = f"{track_path}.sorted.bed"
            cmd = f"sort -k1,1 -k2,2n {quote(track_path)} > {quote(temp_sorted)}"
            subprocess.check_call(cmd, shell=True)
            cmd = f"bgzip {quote(temp_sorted)}"
            subprocess.check_call(cmd, shell=True)
            cmd = f"tabix -C -p bed {quote(sorted_file)}"
            subprocess.check_call(cmd, shell=True)
            print(f"Sorted and indexed BED: {sorted_file}")
        else:
            print(f"Skipping BED sorting/indexing (already exists): {sorted_file}")
        return sorted_file
        
    elif track_format == 'bam':
        # Index BAM files
        index_file = f"{track_path}.csi"
        if not os.path.exists(index_file):
            cmd = f"samtools index -c {quote(track_path)}"
            subprocess.check_call(cmd, shell=True)
            print(f"Indexed BAM: {index_file}")
        else:
            print(f"Skipping BAM indexing (already exists): {index_file}")
        return track_path
        
    elif track_format == 'bigwig':
        # BigWig files don't need additional processing
        return track_path
        
    else:
        print(f"Unknown format {track_format}, returning original path")
        return track_path


def get_track_config(track_record: Dict, track_path: str, genome_name: str) -> Dict:
    """
    Generate JBrowse2 track configuration based on track type and format.
    """
    track_format = track_record['Format']
    track_type = track_record['type']
    label = track_record['Label']
    category = track_record['Category']
    
    # Generate unique track ID
    track_id = hashlib.md5(f"{genome_name}_{label}_{track_path}".encode()).hexdigest()
    track_id = f"{label}_{track_id[:8]}"
    
    base_config = {
        "trackId": track_id,
        "name": label,
        "assemblyNames": [genome_name],
        "category": [category]
    }
    
    if track_format in ['gff3', 'gff']:
        config = {
            **base_config,
            "type": "FeatureTrack",
            "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                    "uri": f"./{os.path.basename(track_path)}",
                    "locationType": "UriLocation"
                },
                "index": {
                    "indexType": "CSI",
                    "location": {
                        "locationType": "UriLocation",
                        "uri": f"./{os.path.basename(track_path)}.csi"
                    }
                }
            }
        }
        
    elif track_format == 'bed':
        config = {
            **base_config,
            "type": "FeatureTrack",
            "adapter": {
                "type": "BedTabixAdapter",
                "bedGzLocation": {
                    "uri": f"./{os.path.basename(track_path)}",
                    "locationType": "UriLocation"
                },
                "index": {
                    "indexType": "CSI",
                    "location": {
                        "locationType": "UriLocation",
                        "uri": f"./{os.path.basename(track_path)}.csi"
                    }
                }
            }
        }
        
    elif track_format == 'bam':
        config = {
            **base_config,
            "type": "AlignmentsTrack",
            "adapter": {
                "type": "BamAdapter",
                "bamLocation": {
                    "uri": f"./{os.path.basename(track_path)}",
                    "locationType": "UriLocation"
                },
                "index": {
                    "indexType": "CSI",
                    "location": {
                        "locationType": "UriLocation",
                        "uri": f"./{os.path.basename(track_path)}.csi"
                    }
                }
            }
        }
        
    elif track_format == 'bigwig':
        config = {
            **base_config,
            "type": "QuantitativeTrack",
            "adapter": {
                "type": "BigWigAdapter",
                "bigWigLocation": {
                    "uri": f"./{os.path.basename(track_path)}",
                    "locationType": "UriLocation"
                }
            }
        }
        
    else:
        raise ValueError(f"Unsupported track format: {track_format}")
    
    return config


def create_synteny_tracks(genomes: Dict[str, List[Dict]], output_dir: str, 
                         make_paf_script: str) -> List[Dict]:
    """
    Create synteny tracks by converting BLAST files to BED and generating PAF files.
    
    Returns list of synteny track configurations.
    """
    synteny_tracks = []
    synteny_records = {}
    
    # Extract synteny records for each genome
    for genome_name, tracks in genomes.items():
        for track in tracks:
            if track['type'] == 'SyntenyTrack':
                synteny_records[genome_name] = track
                break
    
    if len(synteny_records) < 2:
        print("Need at least 2 genomes with SyntenyTrack to create synteny tracks")
        return []
    
    # Convert BLAST files to BED format
    bed_files = {}
    for genome_name, synteny_record in synteny_records.items():
        blast_file = os.path.join(synteny_record['Dirname'], synteny_record['Filename'])
        bed_file = os.path.join(output_dir, f"{genome_name}_synteny.bed")
        
        # Check if BLAST file exists
        if not os.path.exists(blast_file):
            print(f"Warning: BLAST file not found for {genome_name}: {blast_file}")
            print(f"Skipping synteny tracks involving {genome_name}")
            continue
        
        if not os.path.exists(bed_file):
            blast_to_bed(blast_file, bed_file)
        else:
            print(f"Skipping BLAST to BED conversion (already exists): {bed_file}")
        
        bed_files[genome_name] = bed_file
    
    # Generate PAF files for all genome pairs (only for genomes with valid BED files)
    genome_names = list(bed_files.keys())  # Only use genomes with valid BED files
    for i in range(len(genome_names)):
        for j in range(i + 1, len(genome_names)):
            genome1 = genome_names[i]
            genome2 = genome_names[j]
            
            # Skip if either genome doesn't have a BED file
            if genome1 not in bed_files or genome2 not in bed_files:
                continue
            
            bed1 = bed_files[genome1]
            bed2 = bed_files[genome2]
            
            # Get reference FAI files
            fai1 = os.path.join(output_dir, f"{genome1}.fasta.fai")
            fai2 = os.path.join(output_dir, f"{genome2}.fasta.fai")
            
            # Generate PAF file
            paf_file = os.path.join(output_dir, f"{genome1}_{genome2}_synteny.paf")
            
            if not os.path.exists(paf_file):
                cmd = [
                    'Rscript', make_paf_script,
                    '-a', bed1,
                    '-b', bed2,
                    '-A', fai1,
                    '-B', fai2,
                    '-o', paf_file
                ]
                print(f"Generating PAF: {' '.join(cmd)}")
                subprocess.run(cmd, check=True)
            else:
                print(f"Skipping PAF generation (already exists): {paf_file}")
            
            # Create chained PAF file
            paf_chained = os.path.join(output_dir, f"{genome1}_{genome2}_synteny_chained.paf")
            if not os.path.exists(paf_chained):
                merge_paf_intervals(paf_file, paf_chained, tolerance_percent=15, min_intervals=10)
                print(f"Created chained PAF: {paf_chained}")
            else:
                print(f"Skipping PAF chaining (already exists): {paf_chained}")
            
            # Create synteny track configuration
            track_id = f"synteny_{genome1}_{genome2}"
            synteny_config = {
                "type": "SyntenyTrack",
                "trackId": track_id,
                "name": f"{genome1} vs {genome2} Synteny",
                "assemblyNames": [genome1, genome2],
                "adapter": {
                    "type": "PAFAdapter",
                    "targetAssembly": genome2,
                    "queryAssembly": genome1,
                    "pafLocation": {
                        "locationType": "UriLocation",
                        "uri": f"./{os.path.basename(paf_chained)}"
                    }
                },
                "category": ["Synteny"]
            }
            
            synteny_tracks.append(synteny_config)
    
    return synteny_tracks


def create_rdna_tracks(genomes: Dict[str, List[Dict]], output_dir: str) -> List[Dict]:
    """
    Create rDNA tracks for genomes where rdna.bed files exist in the 'Dirname' directory.
    
    Returns list of rDNA track configurations.
    """
    rdna_tracks = []
    
    for genome_name, tracks in genomes.items():
        # Get the dirname from the reference track
        ref_dirname = None
        for track in tracks:
            if track['type'] == 'reference':
                ref_dirname = track['Dirname']
                break
        
        if not ref_dirname:
            print(f"Warning: No reference track found for genome {genome_name}, skipping rDNA check")
            continue
        
        # Check if rdna.bed exists in the dirname directory
        rdna_bed_path = os.path.join(ref_dirname, 'rdna.bed')
        if os.path.exists(rdna_bed_path):
            print(f"Found rDNA file for {genome_name}: {rdna_bed_path}")
            
            # Create rDNA track record
            rdna_record = {
                'Dirname': ref_dirname,
                'Filename': 'rdna.bed',
                'Format': 'bed',
                'Label': 'rDNA',
                'Category': 'assembly',
                'type': 'FeatureTrack'
            }
            
            try:
                # Process the rDNA track
                track_path = sort_and_index_track(rdna_record, genome_name, output_dir)
                track_config = get_track_config(rdna_record, track_path, genome_name)
                
                # Add red color and collapsed display mode
                track_config['displays'] = [{
                    'displayId': f"{track_config['trackId']}_d",
                    'type': 'LinearBasicDisplay',
                    'renderer': {
                        'color1': 'red',
                        'type': 'SvgFeatureRenderer',
                        'displayMode': 'collapse',
                        'showLabels': False
                    }
                }]
                
                rdna_tracks.append(track_config)
                print(f"Added rDNA track for {genome_name}")
                
            except Exception as e:
                print(f"Error processing rDNA track for {genome_name}: {e}")
        else:
            print(f"No rDNA file found for {genome_name} in {ref_dirname}")
    
    return rdna_tracks


def create_jbrowse_config(genomes: Dict[str, List[Dict]], output_dir: str, 
                         make_paf_script: str) -> None:
    """
    Create complete JBrowse2 configuration with all genomes and tracks.
    """
    # Remove existing config.json to ensure fresh creation
    config_file = os.path.join(output_dir, 'config.json')
    if os.path.exists(config_file):
        print('Removing old config.json file')
        os.remove(config_file)
    
    all_tracks = []
    
    # Process each genome
    for genome_name, tracks in genomes.items():
        print(f"\nProcessing genome: {genome_name}")
        
        # Find reference track
        ref_record = None
        other_tracks = []
        
        for track in tracks:
            if track['type'] == 'reference':
                ref_record = track
            else:
                other_tracks.append(track)
        
        if not ref_record:
            print(f"Warning: No reference track found for genome {genome_name}")
            continue
        
        # Create reference assembly (this creates/updates config.json)
        ref_path = create_reference(genome_name, ref_record, output_dir)
        
        # Process other tracks
        for track_record in other_tracks:
            if track_record['type'] == 'SyntenyTrack':
                continue  # Handle synteny tracks separately
            
            try:
                track_path = sort_and_index_track(track_record, genome_name, output_dir)
                track_config = get_track_config(track_record, track_path, genome_name)
                all_tracks.append(track_config)
                print(f"Added track: {track_record['Label']}")
            except Exception as e:
                print(f"Error processing track {track_record['Label']}: {e}")
    
    # Create rDNA tracks
    print("\nCreating rDNA tracks...")
    rdna_tracks = create_rdna_tracks(genomes, output_dir)
    all_tracks.extend(rdna_tracks)
    
    # Create synteny tracks
    print("\nCreating synteny tracks...")
    synteny_tracks = create_synteny_tracks(genomes, output_dir, make_paf_script)
    all_tracks.extend(synteny_tracks)
    
    # Read the existing config.json (created by jbrowse add-assembly commands)
    # and add our tracks to it
    config_file = os.path.join(output_dir, 'config.json')
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
    else:
        # Fallback if no config exists
        config = {
            "assemblies": [],
            "tracks": []
        }
    
    # Add our tracks to the existing configuration
    if "tracks" not in config:
        config["tracks"] = []
    config["tracks"].extend(all_tracks)
    
    # Add configuration section
    config["configuration"] = {"logoPath": {"uri": "./elixir_150x48.svg"}}
    
    # Write updated configuration file
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\nJBrowse2 configuration written to: {config_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Create JBrowse2 configuration for multiple genomes"
    )
    parser.add_argument(
        '-c', '--config-csv', 
        required=True,
        help='CSV configuration file with genome and track information'
    )
    parser.add_argument(
        '-o', '--output-dir', 
        required=True,
        help='Output directory for JBrowse2 files'
    )
    parser.add_argument(
        '--make-paf-script',
        default=None,
        help='Path to make_paf_from_probes.R script (default: look in script directory)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Check if all files exist without creating any output (dry run mode)'
    )
    
    args = parser.parse_args()
    
    # Determine path to make_paf_from_probes.R script
    if args.make_paf_script is None:
        # Get the directory where this script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.make_paf_script = os.path.join(script_dir, 'make_paf_from_probes.R')
    
    # Check if make_paf_script exists
    if not os.path.exists(args.make_paf_script):
        print(f"Error: make_paf_from_probes.R script not found at: {args.make_paf_script}")
        print("Please specify the correct path using --make-paf-script option")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read configuration
    print(f"Reading configuration from: {args.config_csv}")
    genomes = read_config_csv(args.config_csv)
    
    print(f"Found {len(genomes)} genomes:")
    for genome_name, tracks in genomes.items():
        print(f"  {genome_name}: {len(tracks)} tracks")
    
    # If dry run, just check files and exit
    if args.dry_run:
        print("\n=== DRY RUN MODE ===")
        files_exist = check_files_exist(genomes)
        if files_exist:
            print("\n✓ Dry run successful - all required files are present")
            sys.exit(0)
        else:
            print("\n✗ Dry run failed - some files are missing")
            sys.exit(1)
    
    # Create JBrowse2 configuration
    create_jbrowse_config(genomes, args.output_dir, args.make_paf_script)
    
    print(f"\nJBrowse2 setup complete in: {args.output_dir}")


if __name__ == "__main__":
    main()