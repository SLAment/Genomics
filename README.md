# Genomics
A collection of scripts related to Genomics, mostly associated with my work on the genomics of *Podospora anserina*.

I often use [conda](https://docs.conda.io/projects/conda/en/latest/) to install the dependencies of these scripts. For example, after installation of conda, I make an environment like this:

    $ conda create -n bioinformatics blast=2.9.0 biopython=1.78 gffutils=0.10.1 -c bioconda

Then I activate my environment:

    $ conda activate bioinformatics

Notice that the scripts (should) all use python 3. Tested only in Unix environments.

## BLAST

The scripts are made to either parse the BLAST output, filter it or modify it, or to BLAST a query sequence from a given genome.

- `BLAST_tabularparser.py` - The input is a BLAST tabular file produced with the typical `blastn` coomand using `-outfmt 6`.
- `BLASTstitcher.py` - Often the output of BLAST hits is fragmented into several lines that could be collapsed into a single one. This script does just that. The input is also the BLAST tabular file produced with `-outfmt 6`.
- `genometblastn.py` - Small script to do tblastn searches using an input query protein fasta and an input genome fasta.
- `query2haplotype.py` - Extract haplotypes (in fasta) of an assembly based on an input query fasta. It makes a blastn search to extract the hits sequences from the genome. If the `--haplo` option is used, then it searches for a haplotype instead using the query sequences as flanking regions.
- `query2hitseq.py` - An earlier version of `query2haplotype.py` without the `--haplo` option (so it only gives the best hits in fasta format), but it has one extra option to get only the top hit (`--tophit`) which I didn't bother implementing in `query2haplotype.py`.

## Fasta Manipulation

- `fastaconcat.py` - Script to concatenate any number of fasta files, either by name (`--name`) or by rank (default).
- `nexus2markers.py` - Script to extract all partitions in a simple nexus file as individual fasta files.
- `subsetfastaIDbioStdout.py` - Script to produce a subset of an input fasta file using one or several strings that identify specific sequences. The output goes to standard output. It depends on biopython.

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

- `gtfRM2gff.py` - Script to transform the output of RepeatMasker (obtained with option `-gff` in RepeatMasker, which has a misleading name because it's a gtf) into a normal gff3. It also appends a color attribute to normal repeats (`-c`) and to simple repeats (`-s`) to be displayed in [The Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/). 
- `totalcovergff.py` - Script to obtained the merged coordinates of all models in gff file (eg. from RepeatMasker either the gtf or gff produced with `gtfRM2gff.py`). Basically it produces a bed file from the gff file with overlapping features merged. If only the gff is given, then it will collapse all the repeats into non-overlapping intervals. If an associated fasta file is also provided (`--fasta`), then it will calculate the total coverage of the contigs within that fasta annotated in the gff. Notice that the "bed" file produced is in the base of the input file (base 1 with gtf or gff files).

Example: I have an SPAdes assembly of the strain PaZp (`PaZp.nice.fa`), with an annotation file of repeated elements called `PaZp.repeatmasker.gff3`. I want to know what percentage of the genome assembly is repeats. When I produced the assembly I aligned the contigs to a reference, so I assigned them to chromosomes or to the mitochondrion if it they were large enough. Hence, I can exclude the mitochondrial contig from this calcualtion by using the substring `_mt`.

	$ python totalcovergff.py PaZp.repeatmasker.gff3 -f PaZp.nice.fa -E _mt | grep 'Total'
	Total	35798243	1592917	4.450

I used `grep 'Total'` because the script will print the values per contig too but I'm not interested in that. So 4.450% of this genome is annotated as repetitive elements. 
 

- `MFannot4ncbi.oy` - The tbl output of [MFannot](https://www.frontiersin.org/articles/10.3389/fpls.2023.1222186/full) is too bare. This scripts attempts to make it closer to the desired tbl for [`table2asn`](https://www.ncbi.nlm.nih.gov/genbank/table2asn/), the NCBI script that produces sqn files for genome submission.

I ran it as

	% ./MFannot4ncbi.py mfannot_output.fasta.new.tbl MyMtContig.fa -l XXXX -s 10 -n 114110 > MyMtContig.tbl

Where XXXX is the NCBI locus_tag, and 114110 is a number I chose to start the gene codes for all the genes in the mitochondrial contig. This could be 1, or whatever you want, but if you already have nuclear genome contigs with numbers, you can continue those numbers for your mitochondrial contig.

- `tbl2gff.py` - Experimental little script to transform a tbl file, like the one produced by `MFannot4ncbi.py`, into a gff3 file.

- `gff3TOtbl.py` - A simple script to transform a basic gff3 to a tbl file, similar to the one produced by `MFannot4ncbi.py`. It can't handle multiple mRNA isoforms.

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
