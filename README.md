======
TIGAR2
======

TIGAR2: sensitive and accurate estimation of transcript isoform expression with longer RNA-Seq reads

Naoki Nariai, Kaname Kojima, Takahiro Mimori, Yukuto Sato, Yosuke Kawai, Yumi Yamaguchi-Kabata and Masao Nagasaki

*BMC Genomics*, 15(Suppl 10):S5 (2014).

<b>Latest news</b><br>
Feb 19, 2015: TIGAR2.1 was released. Multi-threading is now available.<br>
Dec 12, 2014: TIGAR2.0 was released.

Please download the jar file by clicking <b>Download ZIP</b> on the right panel.

<pre>
Usage: java -jar Tigar2_1.jar FASTA SAM OUT
 
 FASTA           : reference FASTA file
 SAM             : target SAM/BAM file
 OUT             : output file
 
 Options:

 --thread_num INT : number of thread
 --alpha_zero DOUBLE : tuning parameter alpha_zero
 --is_paired     : paired-end data. default = FALSE.
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

<b>2. Build bowtie2 index</b>

<pre>
mkdir ref
bowtie2-build refMrna.fa ./ref/refMrna
</pre>

<b>3. Run bowtie2</b>

For single-end data
<pre>
bowtie2 -p 8 -k 100 --very-sensitive -x ./ref/refMrna sample.fastq > sample.sam
</pre>

For paired-end data
<pre>
bowtie2 -p 8 -k 100 --very-sensitive -x ./ref/refMrna -1 sample_1.fastq -2 sample_2.fastq > sample.sam
</pre>

<b>4. Run TIGAR2</b>

For single-end data
<pre>
java -jar Tigar2_1.jar --thread_num 8 refMrna.fa sample.sam --alpha_zero 0.1 sample_out.txt
</pre>

For paired-end data
<pre>
java -jar Tigar2_1.jar --thread_num 8 refMrna.fa sample.sam --is_paired --alpha_zero 0.1 sample_out.txt
</pre>

<b>Output format</b>

<pre>
ID: transcript (mRNA) ID that the program predicted

LENGTH: transcript length

Z: the number of expected fragments that the program assigned to the transcript

FPKM: normalized expression level (Fragments Per Kilobase of exon per Million mapped fragments)

THETA: estimated parameter (transcript abundance), essentially Z divided by total fragments.
</pre>

<b>5. Visualization</b>

You can visualize the optimized alignment by TIGAR2 as follows:

<pre>
samtools sort sample_out.txt.opt.bam sample_opt_sorted
samtools index sample_opt_sorted.bam
</pre>

Please start IGV_2.3.14 or later, and load refMrna.fa as Genome, and sample_opt_sorted.bam.
You can look at the optimized alignment of reads on each transcript isoform.

<hr>

Please note that the current implementation of TIGAR2 might require large memory size for large sam/bam files.
In such cases, please specify:
<pre>
e.g.) java -Xmx16g -Xms16g -jar Tigar2_1.jar --thread_num 8 FASTA SAM OUT --alpha_zero 0.1
e.g.) java -Xmx32g -Xms32g -jar Tigar2_1.jar --thread_num 8 FASTA SAM OUT --alpha_zero 0.1
e.g.) java -Xmx64g -Xms64g -jar Tigar2_1.jar --thread_num 8 FASTA SAM OUT --alpha_zero 0.1

</pre>



You can also choose BWA-MEM as an aligner as follows:

<b>* Build FM-index for BWA-MEM alignment</b>

<pre>
bwa index refMrna.fa
</pre>

<b>* Run BMA-MEM</b>

For single-end data
<pre>
bwa mem -t 8 -L 10000 -a refMrna.fa sample.fastq > sample.sam
</pre>

For paired-end data
<pre>
bwa mem -t 8 -P -L 10000 -a refMrna.fa sample_1.fastq sample_2.fastq > sample.sam
</pre>

<br>
<b>Please DO NOT sort sam files</b>
<br>

Then, you can run TIGAR2 exactly the same as the way described previously.

<hr>
<b>References</b>
<pre>
Naoki Nariai, Osamu Hirose, Kaname Kojima, Masao Nagasaki.
TIGAR: transcript isoform abundance estimation method with gapped alignment of RNA-Seq data by variational Bayesian inference.
Bioinformatics, 29(18):2292-2299 (2013).

Naoki Nariai, Kaname Kojima, Takahiro Mimori, Yukuto Sato, Yosuke Kawai, Yumi Yamaguchi-Kabata and Masao Nagasaki
TIGAR2: sensitive and accurate estimation of transcript isoform expression with longer RNA-Seq reads.
BMC Genomics, 15(Suppl 10):S5 (2014).
</pre>

This site is maintained by:
Naoki Nariai<br>
<br>
Contact:<br>
nariai [at] megabank.tohoku.ac.jp

Last updated on 2016/07/01

