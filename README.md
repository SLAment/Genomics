# Genomics
A collection of scripts related to Genomics, mostly associated with my work on the genomics of *Podospora anserina*.

I often use [conda](https://docs.conda.io/projects/conda/en/latest/) to install the dependencies of these scripts. For example, after installation of conda, I make an environment like this:

    $ conda create -n bioinformatics blast=2.9.0 biopython=1.78 gffutils=0.10.1 -c bioconda

(The BLAST scripts below have been tested and work well with BLAST 2.15.0+, biopython 1.81 and gffutils 0.12).

Then I activate my environment:

    $ conda activate bioinformatics

Notice that the scripts (should) all use python 3. Tested only in Unix environments.

## BLAST

The scripts are made to either parse the BLAST output, filter it or modify it, or to BLAST a query sequence from a given genome.

- `BLAST_tabularparser.py` - The input is a BLAST tabular file produced with the typical `blastn` coomand using `-outfmt 6`.

- `BLASTstitcher.py` - Often the output of BLAST hits is fragmented into several lines that could be collapsed into a single one. This script does just that. The input is also the BLAST tabular file produced with `-outfmt 6`.

- `genometblastn.py` - Small script to do tblastn searches using an input query protein fasta and an input genome fasta. (Depricated in favor of `query2haplotype.py`).

- `query2haplotype.py` - Extract haplotypes (in fasta) of an assembly based on an input query fasta. It makes a blastn (or tblastn if requested with `--task`) search to extract the hits sequences from the genome. If the `--haplo` option is used, then it searches for a haplotype instead using the query sequences as flanking regions.

- `query2hitseq.py` - An earlier version of `query2haplotype.py` without the `--haplo` option (so it only gives the best hits in fasta format). (Depricated in favor of `query2haplotype.py`).

## Fasta Manipulation

- `fastaconcat.py` - Script to concatenate any number of fasta files, either by name (`--name`) or by rank (default).

- `gc_calc.py` - Script to calculate GC content from a (multi) fasta file (reporting for each individual sequence)

- `nexus2markers.py` - Script to extract all partitions in a simple nexus file as individual fasta files.

- `psim_calc.py` - Script to calculate Psim sensu [Jorda and Kajava (2009)](https://academic.oup.com/bioinformatics/article/25/20/2632/193638) in the program T-REKS.

- `purgeFasta.py` - Script to subset or remove sequences in an input fasta file using a one-column file with the name of the sequences (to either remove or extract). The output goes to standard output. It depends on biopython.

- `subsetfastaIDbioStdout.py` - Earlier version of `purgeFasta.py`. Script to produce a subset of an input fasta file using one or several strings that identify specific sequences. The output goes to standard output. It depends on biopython.

## Genome Annotation

- `GFFnumerator.py` - Script to re-name the IDs of the genes and features of a GFF3 file. It prints to the standard output. It depends on the library `gffutils`. **OBS** The script CAN'T deal with genes with multiple isoforms (mRNAs). If forced, it might fused exons were it shouldn't. Use with care.

- `GFFSlicer.py` -  Script to extract sections of a GFF file while correcting their coordinates, making them relative to the start of a given coordinate. The default is to assume there is a single contig in the gff, but if not then the contig of interest can be specified (`--contig`). There is also an experimental option to reverse complement a gff (`--invertcoords`) but I haven't implemented a way to fix the frame of CDS features, so those break.

- `GFFSubset.py` - Script to extract features of a gff into a new gff based on name or ID. It relies on `gffutils`. It assumes that the higher level feature is gene.

- `gffutils2fasta.py` - Script to extract fasta sequences out of an input fasta file using a corresponding GFF3 file. This is useful if you want for example the protein sequences of an annotation, or just the exons, etc. The script can extract the following types of features (case sensitive):
	
	* gene
	* CDS - It can also be written as "cds"
	* exon
	* noutrs - this basically means to return a gene (including the introns) with the UTRs removed if present.
	* similarity - manually modified version of the RepeatMasker gtf (might remove in the future)
	* expressed_sequence_match - this is found in the output of MAKER as alignment of other proteins
    * repeat - modified version of the RepeatMasker gtf as produced by `gtfRM2gff.py`

**OBS** The script CAN'T deal with genes with multiple isoforms (mRNAs). If forced, it might fuse exons were it shouldn't. Use with care.

- `gff3TOtbl.py` - A simple script to transform a basic gff3 to a tbl file, similar to the one produced by `MFannot4ncbi.py`. It can't handle multiple mRNA isoforms. This is ideal if you make manual curation into your gff3 (as produced by `tbl2gff.py`, for example) and need to produce a tbl file again. I did this for the output of MFannot, so the script has a particular option to deal with the endonucleases and introns annotated by this program: they get annotated as mobile_element features. Otherwise they are treated like protein-coding genes (without `--mfannot`).

It's ran like so:

	$ ./gff3TOtbl.py NyNewAnnotation.gff3 MyMtContig.fa --mfannot > NyNewAnnotation.tbl

The script can also make features partial, although this require testing. Notice that NCBI won't accept partial features in fully assembled mitochondrial genomes.

	$ ./gff3TOtbl.py NyNewAnnotation.gff3 MyMtContig.fa --mfannot --partial > NyNewAnnotation.tbl

- `gtfRM2gff.py` - Script to transform the output of RepeatMasker (obtained with option `-gff` in RepeatMasker, which has a misleading name because it's a gtf) into a normal gff3. It also appends a color attribute to normal repeats (`-c`) and to simple repeats (`-s`) to be displayed in [The Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/). 

- `invertingBEDs.py` - Script to invert the coordinates of a BED file along a scaffold. It works with other file types as long as there is a scaffold and at least one coordinate column. This is useful when you wish you had done your analysis on the reverse-complement sequence of your reference fasta.

You require two input files: the BED file and a supporting file used to get the total length of each contig. The supporting file can be a two-columns file in the stile of:

	contig1	100000
	contig2	50000

Or a fasta file with the same contigs.

Use with a standard BED file:

	$ ./invertingBEDs.py myfile.bed mycolumnsfile.txt > output.bed

With a fasta file

	$ ./invertingBEDs.py myfile.bed assembly.fasta --fasta > output.bed

With a simple gtf file, you have to tell it what columns contain the contigs (`-c`), the start coordinates (`-s`), and the end coordinates (`-e`).

	$ ./invertingBEDs.py myfile.gtf assembly.fasta --fasta -c 1 -s 4 -e 5 > output.gtf

**OBS** this won't correct the sense of complicated features like protein-coding genes. See `GFFSlicer.py` for a script more dedicated to deal with GFF3 files.

You can also use it for files with a single coordinate column (e.g. VCF files), just give the same column to `-s` and `-e`. (less tested)

- `MFannot4ncbi.py` - The tbl output of [MFannot](https://www.frontiersin.org/articles/10.3389/fpls.2023.1222186/full) is too bare. This scripts attempts to make it closer to the desired tbl for [`table2asn`](https://www.ncbi.nlm.nih.gov/genbank/table2asn/), the NCBI script that produces sqn files for genome submission.

I ran it as

	$ ./MFannot4ncbi.py mfannot_output.fasta.new.tbl MyMtContig.fa -l XXXX -s 10 -n 114110 > MyMtContig.tbl

Where XXXX is the NCBI locus_tag, and 114110 is a number I chose to start the gene codes for all the genes in the mitochondrial contig. This could be 1, or whatever you want, but if you already have nuclear genome contigs with numbers, you can continue those numbers for your mitochondrial contig.

- `tbl2gff.py` - Experimental little script to transform a tbl file, like the one produced by `MFannot4ncbi.py`, into a gff3 file.

- `TideCluster2RM.py` - Script to process the output of TideCluster (`tarean_report.tsv` and `trc_superfamilies.csv`) into a repeat library compatible with RepeatMasker. The superfamilies are optional, and they get added as  `TRCxx#Satellite/sfY` for a xx family that belongs to a Y superfamily.



## Miscellaneous

TODO

## Pylogenetics

- `orthogrs2fasta.py` - Script to parse the *Orthogroups.csv* and *SingleCopyOrthogroups.txt* files produced by [OrthoFinder](https://github.com/davidemms/OrthoFinder) and to extract each orthogroup in a fasta file.

- `orthogrs_parser.py` - Script to parse the *Orthogroups.csv* output file of [OrthoFinder](https://github.com/davidemms/OrthoFinder) and manage it for the *Podospora* project. Useful if you want to find orthogroups where all samples have a given number of orthologs (as defined by `--nugrps`). It can also sample of a given number of orthologs groups randomly (`--sample`). In addition, the script can look for certain type of orthogroups that I call "cool". An orthogroup is cool if:
    
    * All species are represented
    * It excludes in-paralogs (duplications within species) that are only present in one species
    * It includes groups of paralogs, but not all of the species/samples have to be represented in the clade of each paralog

They are cool because, in theory, shared duplications (paralogs) are phylogenetically informative.

## PopGen

Scripts for population genetics, mostly vcf file manipulation.

- `vcfNAreport.py` - Little script to report ammount of missing data in a vcf file per sample. 

For example, you could look for samples that have more than 1% missing sites like:

	$ python vcfNAreport.py Example.vcf.gz -z | awk '$4 > 0.01'
	Sample  Total_sites     Missing_sites   Missing%
	Sample1    35684   34713   0.97279
	Sample2    35684   1346    0.03772
	Sample3    35684   440     0.01233


----

*Disclaimer:* These scripts and example files are provided "as is" and without any express or implied warranties, including, without limitation, the implied warranties of merchantability and fitness for a particular purpose.
