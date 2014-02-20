======
TIGAR2
======

TIGAR2: sensitive and accurate transcript isoform abundance estimation with longer RNA-Seq reads

Naoki Nariai, Kaname Kojima, Takahiro Mimori, Yukuto Sato, Yosuke Kawai, Yumi Yamaguchi-Kabata, Masao Nagasaki

Submitted.

<pre>
 Example: java -jar Tigar2.jar FASTA SAM OUT --alpha_zero DOUBLE --is_paired --polyA
 FASTA          : reference FASTA file
 SAM            : target SAM/BAM file
 OUT            : output file
 --alpha_zero DOUBLE : tuning parameter alpha_zero
 --is_paired    : paired-end data. default = 0 (false). Please set 1, if sam
                  file was generated from paired-end reads.
 --polyA BOOLEAN : polyA flag. default = 0 (false). Please set 1 if both read
                  and reference sequences contain polyA tails.
</pre>

## Recommended pipeline to run TIGAR2

Prepare cDNA reference sequences in FASTA format.

<pre>
e.g.) human
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz

e.g.) mouse
http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/refMrna.fa.gz
</pre>

Build FM-index for alignment

<pre>
bwa index refMrna.fa
</pre>

Run bowtie2

For single-end data
<pre>
bwa mem -t 8 -P -L 10000 -a refMrna.fa sample.fastq > sample.sam
</pre>

For paired-end data
<pre>
bwa mem -t 8 -P -L 10000 -a refMrna.fa sample_1.fastq sample_2.fastq > sample.sam
</pre>

Please note that sam files are expected to be sorted by read name.
In order to sort sam files by read name:

<pre>
samtools view -bS sample.sam > sample.bam
samtools sort -n sample.bam sample_sorted
</pre>

Run TIGAR

<pre>
java -jar Tigar2.jar refMrna.fa sample_sorted.bam --alpha_zero 0.1 sample_out.txt
</pre>

Output format

<pre>
ID: transcript (mRNA) ID that the program predicted

LENGTH: transcript length

Z: the number of expected fragments that the program assigned to the transcript

FPKM: normalized expression level (Fragments Per Kilobase of exon per Million mapped fragments)

THETA: estimated parameter (transcript abundance), essentially Z divided by total mapped reads.
</pre>

Please note that the current implementation of TIGAR requires large memory size for large sam/bam files.
For large sam/bam files, please specify:
<pre>
java -Xmx16g -Xms16g FASTA SAM OUT --alpha_zero DOUBLE --is_paired INT --polyA INT
</pre>
or, ideally
<pre>
java -Xmx32g -Xms32g FASTA SAM OUT --alpha_zero DOUBLE --is_paired INT --polyA INT
</pre>

This site is maintained by:
Naoki Nariai<br>
<br>
Contact:<br>
nariai [at] megabank.tohoku.ac.jp

Last updated on 2014/02/20

