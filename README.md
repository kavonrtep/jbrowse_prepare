# scripts for preparation of jbrowse directories and config

Script uses csv table as input to create and formate reference and track. For coping to remote apollo server use `rsync` command with `-L -r --update --progress` option to ensure that symlink are transformed into referent file/dir.

## example of input csv table:


| **label**                  | **format** | **category**    | **type**       | **dirname**                                                                                     | **filename**                                                    | **color** |
| ---------------------------- | ------------ | ----------------- | ---------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- | ----------- |
| CEPIT                      | fasta      | reference       | reference      | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203                      | asm.bp.hap1+2.p.ctg.fa                                          |           |
| stringtie_transdecoder_met | gff3       | gene_annotation | HTMLFeatures   | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203/stringtie            | annotation_stringtie_conservative_with_CDS_MET_clean_eggnog.gff |           |
| stringtie_transdecoder     | gff3       | gene_annotation | HTMLFeatures   | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203/stringtie            | annotation_stringtie_conservative_with_CDS_clean_eggnog.gff     |           |
| kinetochore_gth            | gff3       | gene_annotation | CanvasFeatures | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203/kinetochore_proteins | kinetochore_proteins_selection_for_annot_01_gth.gff             |           |
| Coleus-F1                  | bam        | expression      | Alignments2    | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203/RNA-seq_mapping      | CEPIT-Coleus-F1Aligned.sortedByCoord.out.bam                    |           |
| Betonica-F2                | bam        | expression      | Alignments2    | /mnt/raid/454_data/cuscuta/Cepithymum_assembly_v4/0_final_hifiasm_20211203/RNA-seq_mapping      | CEPIT-Betonica-F2Aligned.sortedByCoord.out.bam                  |           |

first row is header, second row must be reference genome!

install requirements:

conda install -c conda-forge -c biconda jbrowse \
ucsc-fatotwobit gff3sort tabix samtools blast blat \
bioconductor-rtracklayer bioconductor-biostrings bioconductor-bsgenome \
bzip2
Included script:

- `bedpe2bed.py`convert six column bedpe to three column bed
- `create_jbrowse_dir.py`  create data directory and config for Jbrowse1
- `create_jbrowse2_config.py` create data directory and `json.conf`for Jbrowse 2

## Example of usage for Jbrowse2

```bash

mkdir genome_directory
cd genome_directory
$PATH_TO_JBROWSE_PREPARE/create_jbrowse2_config.py -d . -c $PATH_TO_CONFIG_CSV

``` 
