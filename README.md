======
TIGAR2
======

TIGAR2: sensitive and accurate estimation of transcript isoform expression with longer RNA-Seq reads

Naoki Nariai, Kaname Kojima, Takahiro Mimori, Yukuto Sato, Yosuke Kawai, Yumi Yamaguchi-Kabata and Masao Nagasaki

*Submitted*.

Please download the jar file by clicking <b>Download ZIP</b> on the right panel.

<pre>
Usage: java -jar Tigar2.jar FASTA BAM OUT
 
 FASTA           : reference FASTA file
 BAM             : target SAM/BAM file
 OUT             : output file
 
 Options:
 
 --alpha_zero DOUBLE : tuning parameter alpha_zero
 --is_paired     : paired-end data. default = FALSE.
 --polyA         : polyA flag. default = FALSE.
 --frag_dist_mean DOUBLE: mean of the fragment length distribution. default = estimation from data
 --frag_dist_std DOUBLE:  standard dev of the fragment length distribution. default = estimation from data
</pre>

## Recommended pipeline to run TIGAR2

<b>1. Prepare cDNA reference sequences in FASTA format.</b>

<pre>
e.g.) human
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz

e.g.) mouse
http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/refMrna.fa.gz
</pre>

<b>2. Build FM-index for alignment</b>

<pre>
bwa index refMrna.fa
</pre>

<b>3. Run BMA-MEM</b>

For single-end data
<pre>
bwa mem -t 8 -P -L 10000 -a refMrna.fa sample.fastq > sample.sam
</pre>

For paired-end data
<pre>
bwa mem -t 8 -P -L 10000 -a refMrna.fa sample_1.fastq sample_2.fastq > sample.sam
</pre>


<b>4. Run TIGAR2</b>

For single-end data
<pre>
java -jar Tigar2.jar refMrna.fa sample.bam --alpha_zero 0.1 sample_out.txt
</pre>

For paired-end data
<pre>
java -jar Tigar2.jar refMrna.fa sample.bam --is_paired --alpha_zero 0.1 sample_out.txt
</pre>


<b>Output format</b>

<pre>
ID: transcript (mRNA) ID that the program predicted

LENGTH: transcript length

Z: the number of expected fragments that the program assigned to the transcript

FPKM: normalized expression level (Fragments Per Kilobase of exon per Million mapped fragments)

THETA: estimated parameter (transcript abundance), essentially Z divided by total mapped reads.
</pre>

Please note that the current implementation of TIGAR2 might requir large memory size for large sam/bam files.
In such cases, please specify:
<pre>
e.g.) java -Xmx16g -Xms16g -jar Tigar2.jar FASTA SAM OUT --alpha_zero 0.1
e.g.) java -Xmx32g -Xms32g -jar Tigar2.jar FASTA SAM OUT --alpha_zero 0.1
e.g.) java -Xmx64g -Xms64g -jar Tigar2.jar FASTA SAM OUT --alpha_zero 0.1
</pre>

Please also note that sam files are expected to be sorted by read name.
In order to sort sam files by read name:

<pre>
samtools view -bS sample.sam > sample.bam
samtools sort -n sample.bam sample_sorted
</pre>


This site is maintained by:
Naoki Nariai<br>
<br>
Contact:<br>
nariai [at] megabank.tohoku.ac.jp

Last updated on 2014/03/06

