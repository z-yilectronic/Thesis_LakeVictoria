# MetagenomicLakeVictoria
For my Msc thesis on sedimentary DNA from Lake Victoria. Lake Victoria is a shallow tropical lake in Africa, having experienced desiccation and refilling before the modern lake finally formed over the past 15k years.

The commands are based on the ancient metagenomic analysis pipeline Holi(https://github.com/miwipe/Holi/tree/main/bin)

# Sequencing data
Reference genomes were from microbial genome database such as GTDB and IMG/VR, used for mapping bacteria, virus and archaea. Raw paired-end reads sequencing was performed on Illumina HiSeq 4000 platform. Quality control of raw reads involved trimming, filtering, and merging using Fastp.

Command for reference:
```bash
# Trimming and merging reads with fastp
fastp -i *R1*.fastq.gz -I *R2*.fastq.gz -m --merged_out '.ppm.fq' -V --detect_adapter_for_pe -D --dup_calc_accuracy 5 \
      -g -x -q 30 -e 25 -l 30 -y -c -p -h '.fastp.report.html' -w 1
# Removing duplicates with vsearch
vsearch --fastx_uniques '.ppm.fq' --fastqout '.ppm.vs.fq' --minseqlength 30 --strand both
# Filtering low-complexity reads with SGA
sga preprocess --dust-threshold=1 -m 30 '.ppm.vs.fq' -o '.ppm.vs.d1.fq'
```

# Reads mapping
FASTQ format were mapped to the reference databases using Bowtie2. 
```bash
bowtie2 -x "$db_path" -k 1000 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --np 1 --mp "1,1" \
            --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1" --no-unal -U "${fastq}.ppmR1.fq" -S "$bam_output"
```

BAM file headers then were minimized using the Compressbam function from metaDMG. 
```bash
metaDMG-cpp/misc/compressbam --input "$bam" --output "$(basename "$bam" .bam).comp.bam"
```

Sort and merge all compressed BAM files, and compress them again.
```bash
samtools sort -n -m 4G -o $sorted_bam $bam
samtools merge -n -f '.comp.sam.gz' '.comp.bam.sorted.bam'
metaDMG-cpp/misc/compressbam --input '.comp.sam.gz' --output '.comp.bam'
```

# Alignment filtering
This step involved filtering each BAM to detect reliable aligned reads. Filter tool FilterBAM was used to determine the "best" alignments. 

Firstly reassign the multiple-aligned reads. (If there's no error but no reassign file is generated, consider adjust the ani and count threshold.)
```bash
filterBAM reassign --bam "$bam_output" -i 0 --min-read-ani 92 --min-read-count 3 -M 8G -o "$reassigned_bam"
```

Then filter the reassigned BAM and remember to sort them again.
```bash
filterBAM filter --bam "$reassigned_bam" -m 8G --stats "${filename}_stats.tsv.gz" --stats-filtered "${filename}_stats_filtered.tsv.gz" \
        --min-read-ani 92 --min-read-count 100 --min-normalized-entropy 0.6 --min-normalized-gini 0.4 --min-breadth 0.01 --include-low-detection \
        --bam-filtered "$filtered_bam"
samtools sort -n -@ 12 -m 8G -o "$sorted_bam" "$filtered_bam"
```

# Taxon classification
Taxonomic classification step used the tool metaDMG.

MetaDMG-cpp lca module performs the lowest common ancestor (LCA) algorithm to assign reads to taxonomic groups. It provides a fast and flexible classification of DNA reads aligned against reference databases containing multiple organisms. (Classification was performed against a taxonomy. Taxon names and nodes were extracted from the taxdump files and supplemented with a dictionary mapping accession numbers to taxon IDs.)
```bash
metaDMG-cpp lca --threads 12 --bam "$sorted_bam" --nodes "$taxonomy_path/nodes.dmp" --names "$taxonomy_path/names.dmp" \
        --acc2tax "$taxonomy_path/acc2taxid.map.gz" --weight_type 1 --fix_ncbi 0 \
        --sim_score_low 0.6 --how_many 30 --out_prefix "$lca_prefix"
```
For each read, 30 nucleotide positions from both the 5’ and 3’ ends were included in generated substitution matrices which serve for down
stream damage pattern estimation.

# Damage estimate
The next step involves computing observed DNA deamination patterns (C -> T on forward strand, G -> A on reverse strand) and estimating damage of nodes within LCA tree, using the matrix of misincorporated bases.
```bash
# Damage fit calculation
metaDMG-cpp dfit "$lca_prefix.bdamage.gz" --nodes "$taxonomy_path/nodes.dmp" --names "$taxonomy_path/names.dmp" \
        --lib ds --nopt 10 --doboot 1 --nbootstrap 20 --showfits 2 --seed 25487 \
        --out_prefix "$dfit_prefix"
# Aggregate files of tree nodes, taxa names, lca statistics and damage
metaDMG-cpp aggregate "$lca_prefix.bdamage.gz" --nodes "$taxonomy_path/nodes.dmp" --names "$taxonomy_path/names.dmp" --lcastat "$lca_prefix.stat.gz" --dfit "$dfit_prefix.dfit.gz" \
        --out_prefix "$aggregate_prefix"
```
metaDMG dfit performed numerical optimization of a deamination frequencies model based on the mismatch matrices. It estimates three key parameters of the DNA damage model: amplitude of damage at the first position(A), relative decrease in damage per position(q) and background noise(c, equivalent to sequencing errors). 

# aDNA authenication
MetaDMG analysis outputs a particular Briggs damage model for each sample, with deamination frequencies and associated confidence intervals for 30 nucleotide positions on both the forward and reverse strands. 

To assess ancient state, the concordance correlation coefficient (CCC) was calculated to evaluate the agreement between the damage model and the observed nucleotide misincorporation patterns in each species.
