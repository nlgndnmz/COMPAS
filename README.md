## COMPAS - Comparative alternative splicing detection
###### Author: Nilgun Donmez

COMPAS is a bioinformatics tool that is designed to analyze pairs of RNA-Seq datasets in order to identify significant splice variants.

##### Pre-requisites

To use COMPAS, you need Perl ( http://www.perl.org/ ) and R ( http://www.r-project.org/ ) installed on your machine. COMPAS has primarily been tested on 64-bit Linux environments using the following versions (though older or newer versions may also work):

- Perl (5.14.2)
- R (2.15.1)
	
You also need an efficient RNA-Seq mapper (such as, but not limited to, TopHat http://tophat.cbcb.umd.edu/ ) which can report its output in .sam (or .bam) format. The choice of the mapper is left to the user, however we recommend using a software that can both use annotations and find novel junctions. High mapping sensitivity is also recommended.

##### Installation

To run COMPAS, you need to install the "compas" R package (soon to be made available from CRAN, currently available upon request - see compas-manual.pdf for contact info) on your system. The best way to do this is through install.packages() function in R. In addition to this package, COMPAS needs the following scripts:

- processSam.pl
- mergeCounts.pl
- runCOMPAS.R

To call these scripts from anywhere, simply copy them to your ./bin or add the directory containing COMPAS to your path.

##### Running COMPAS

First, you should map the reads for each dataset to the reference genome using an up-to-date gene annotation file (latest Ensembl annotations are recommended). Keep this annotation file (must be in .gtf format) around since you will need it in later steps. Note that if you have reads from multiple lanes and/or paired-end reads in separate fastq files, make sure to map them together into a single .sam file.

For the rest of this section, I will suppose you have two high-throughput RNA-Seq datasets mapped as "tumor.sam" and "benign.sam".

###### Step 1

Run the following commands (can be run in parallel):

```
cat ./tumor.sam | processSam.pl -g ./gene_annotations.gtf -o tumor -r 100 -s 0
cat ./benign.sam | processSam.pl -g ./gene_annotations.gtf -o benign -r 100 -s 0
```

Above, "gene_annotations.gtf" is a GTF formatted annotation file (ideally the same file that was used during the mapping phase). The directive "-r 100" tells the script that the read length is 100bp and "-o" option specifies the prefix for the output files. "-s 0" option tells the script that this RNA-Seq dataset is not strand-specific. This option can be omitted as long as the dataset is not strand-specific. Otherwise, it should be set to either "-s 1" (for forward strand genes) or "-s 2" (for reverse strand genes). The above command will work for .sam files that are coordinate sorted. If your files are not sorted by genomic coordinates, you have to toggle the "-u" option:

` cat ./unsorted_tumor.sam | processSam.pl -g ./gene_annotations.gtf -o tumor -r 100 -s 0 -u `

Note that processSam.pl reads its input from standard input. This is to facilitate pipelines. For example, if you keep your mapping files in compressed .bam format you can do:

` samtools view ./tumor.bam | processSam.pl -g ./gene_annotations.gtf -o tumor -r 100 -s 0 `

This will save significant hard disk space and time since the intermediate .sam file will never be written to disk. In fact, this is highly recommended since it will also be faster than reading the .sam file from disk.

For the first command above, three output files will be written: "tumor.counts", "tumor.hist" and "tumor.skipped". Similarly, for the second command: "benign.counts", "benign.hist" and "benign.skipped" will be written. The .skipped files are for information only, the other files will be used in the subsequent steps.

###### Step 2 

To use COMPAS in (default) comparative mode, you first have to combine the outputs of processSam.pl from different samples. This is done as follows:

mergeCounts.pl -q benign.counts -t tumor.counts > benign_vs_tumor.counts \n\n";

Note that the order of the parameters "-q" and "-t" above does not matter. Note that mergeCounts.pl writes its output to standard output, so directing this into a file (as above) is highly recommended.

###### Step 3

In this step, we run the main R script as below:

` runCOMPAS.R 2 -b 100000000 -i ./benign_vs_tumor.counts -o ./benignVStumor -h ./tumor.hist -l 100 -n 382000000 -h2 ./benign.hist -l2 100 -n2 391000000 > run.log `

Above "benign_vs_tumor.counts" is the file we obtained in the previous step and "benignVStumor" is a prefix for the output files to be written. The arguments "-h", "-l" and "-n" are the .hist file, the read length in base pairs and the total number of reads mapped for the benign dataset (given with -q option in mergeCounts.pl) respectively. "-h2", "-l2" and "-n2" are the same arguments for the second sample. Note that this script writes its progress to the standard output. If you wish to keep this information, direct it to a file as above.

Note that the above command assumes you have already installed the "compas" package on your system. If the run is successful, the following four files will be written:

- benignVStumor.1.gtf
- benignVStumor.2.gtf
- benignVStumor.stat
- benignVStumor.out

The first two files above are GTF formatted files containing all potential transcripts. The .stat file contains some statistics about the optimization run. The .out file contains the predicted major splicing differences between the two samples.  


