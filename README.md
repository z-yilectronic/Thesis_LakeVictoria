# Thesis_LakeVictoria
For the master thesis about aDNA from lake Victoria sediment. 

# Sequencing data
Reference genomes were from microbial genome database such as GTDB and IMG/VR, used for mapping bacteria, virus and archaea. Raw paired-end reads sequencing was performed on Illumina HiSeq 4000 platform. Quality control of raw reads involved trimming, filtering, and merging using Fastp.
Command for reference:
```bash
# Trimming and merging reads with fastp
fastp -i *R1*.fastq.gz -I *R2*.fastq.gz -m --merged_out '.ppm.fq' -V --detect_adapter_for_pe -D --dup_calc_accuracy 5 \
      -g -x -q 30 -e 25 -l 30 -y -c -p -h '.fastp.report.html' -w 1
# Removing duplicates with vsearch
vsearch --fastx_uniques '.ppm.fq' --fastqout '.ppm.vs.fq' --minseqlength 30 --strand both
```

# Reads mapping
FASTQ format were mapped to the reference databases using Bowtie2. When mapping with Bowtie2, a maximum of 1000 valid alignments per read
was allowed, sorted by alignment score. The upper limit on the number of seed
extensions was set to 15, with re-seeding allowed twice, and the seed length was
set to 22. No penalty was applied for seed-gap opening, and a penalty of 1
was applied for gap extension for both the reads and the reference. A penalty
of 1 was also applied for unknown or mismatched nucleotides. A lower score
threshold was used for longer reads.
Only aligned reads were written to the output SAM file. The SAM file was
then converted to binary file (BAM) using Samtools[27] for faster indexing and
reading, as well as to save storage space.
The BAM files were very large, up to tens or even hundreds of thousands of
megabytes. In this case, most BAM file headers were too large to permit merging
all BAM files into a single dataset. The headers were minimized by deleting all
entries that did not have associated mappings, using the Compressbam function
from metaDMG. After simplifying the headers, all compressed BAM files for a
11
single sample were merged into one BAM file using Samtools.
Specifically, microbial mapped reads were filtered out prior to eukaryotic
mapping to reduce potential contamination and interference, thereby improving
the accuracy and resolution of eukaryotic read assignment.
2.4 Alignment filtering
The next step involved filtering each BAM file for every sample ID to detect
reliable aligned reads. The first filter, filterBAM[28], was used to determine the
best, or set of best, alignments.
In the reassigning step, the EM algorithm of filterBAM was not used, which
means filtering was performed directly at a threshold of 92% ANI (average
nucleotide identity) and a minimum of 10 reads per aligned reference.
Following the initial reassign, additional filtering thresholds were applied:
minimum read ANI of 92%, minimum read count per reference of 100, minimum
normalized entropy of 0.6 (normalized entropy of the read coverage distribution),
minimum normalized Gini coefficient of 0.4 (normalized Gini coefficient of the
read coverage distribution), and a minimum coverage breadth of 0.01.
The output BAM files were sorted by read names. The filtering step in
f
ilterBAM also produced statistical outputs in TSV format for further analysis.
2.5 Taxon classification
Taxonomic classification step used the tool metaDMG[29].
MetaDMG-cpp lca module performs the lowest common ancestor (LCA)
algorithm to assign reads to taxonomic groups. LCA reassigns the reads by ex
amining all references (after filtering) to which a read mapped and identifies the
lowest shared taxonomic ancestor. It provides a fast and flexible classification of
DNA reads aligned against reference databases containing multiple organisms.
Classification was performed against NCBI taxonomy.
Taxon names and nodes were extracted from the NCBI taxdump files and
supplemented with a dictionary mapping accession numbers to taxon IDs. A
minimum similarity score of 0.95 (i.e., ≥ 95% matching sites) was required for
an alignment to be considered.
Specifically, microbial mapped reads were filtered out prior to eukaryotic
mapping to reduce potential contamination and interference, thereby improving
the accuracy and resolution of eukaryotic read assignment.
For each read, 30 nucleotide positions from both the 5’ and 3’ ends were
included in generated substitution matrices. These matrices serves for down
stream damage pattern estimation, enabling the assessment of characteristic
nucleotide misincorporations associated with ancient DNA damage.
12
2.6 Damage estimate
The next step involves using the matrix of misincorporated bases to extract
observed DNA deamination patterns (C to T substitutions on the forward strand
and G to A substitutions on the reverse strand) and to estimate damage of
nodes within LCA tree. metaDMG dfit performed numerical optimization of a
deamination frequencies model based on the mismatch matrices. It estimates
three key parameters of the DNA damage model: amplitude of damage at the
f
irst position(A), relative decrease in damage per position(q) and background
noise(c, equivalent to sequencing errors). A binomial distribution was used
as the likelihood model, and 20 bootstrap iterations were performed to assess
uncertainty. The number of optimization calls was 10.
Following the calculations, the LCA classification table and corresponding
damage model fit (dfit) statistics were merged.
2.7 aDNA authenication
MetaDMGanalysis outputs a particular damage model for each sample, with
deamination frequencies and associated confidence intervals for 30 nucleotide
positions on both the forward and reverse strands.
To assess ancient state, the concordance correlation coefficient (CCC) was
calculated to evaluate the agreement between the damage model and the ob
served nucleotide misincorporation patterns in each classified taxon. The fit, as
reflected by the coefficient, was used to determine whether the taxa (identified
through the mapping, filtering, and classification steps mentioned above) ex
hibit ancient DNA damage patterns. It separated modern taxa inhabiting the
sediment layers and support the interpretation of ancient DNA sources. Ancient
taxa were identified and authenticated based on DNA damage patterns using
the Briggs damage model[30].
