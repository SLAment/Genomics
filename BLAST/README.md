# BLAST-related scripts

The most useful script is probably `query2haplotype.py`. I use it to retrieve the hits produced by BLAST in fasta format. The script `query2hitseq.py` can also do that (it's an older version), but the real advantage of `query2haplotype.py` becomes apparent when the option `--haplo` is used. Often, the actual BLAST hit gets fragmented in pieces, in the case of, for example, internal repeats. If you have the BLAST table, you can use `BLASTstitcher.py` to fix some of that. However, you probably want the whole thing in a continuous fasta sequence. In that case `--haplo` fuses all the hits that are too close together (maximum distance controlled by `--vicinity`) into a single "haplotype". The BLAST hits used to create the haplotype can be controlled with `--evalue` and `--minsize`. 

Notice that if the query fasta file contains several sequences, these will be considered as possible 5' or 3' end of the haplotype. This is very convenient when you want a particular locus of your favorite genome and you know the flanking genes. You can take the sequence of the left (5') and right (3') genes of your locus, put them in a fasta file and use it as input query.

It's important not to confuse `--vicinity` with `--minhaplo`. The latter simply filters out the output sequences that are smaller than the value given to `--minhaplo` (default 0 bp). That is useful if you only want the big haplotypes.

You can also ask for extra bp to both sides of the haplotype with `--extrabp`. This is particularly useful when you are doing manual curation of transposable elements and you are looking for the edges of the elements.

For additional options see:

    $ python query2haplotype.py -h


*I might make `--haplo` the default in the future.*

---

For filtering of BLAST tables, you can use `BLAST_tabularparser.py`. As input it expects the table produced with the table format 6 (outfmt=6).

---

Although very useful, `query2haplotype.py` and `query2hitseq.py` only work with nucleotide sequences. I made `genometblastn.py` to do tBLASTn searches using a **protein sequence as query**. It does NOT make haplotypes at the moment.

