# BLAST-related scripts

Remember that all scripts are meant to work with python 3 and require the `biopython` library and BLAST installed.

The most useful script is probably `query2haplotype.py`. I use it to retrieve the hits produced by BLAST in fasta format. Example usage:

    $ python query2haplotype.py path/to/assembly.fasta path/to/queryfile.fasta > hits.fasta

Are you BLASTing something big with internal repeats? Better try:

    $ python query2haplotype.py path/to/assembly.fasta path/to/queryfile.fasta --haplo > hits.fasta

The script `query2hitseq.py` can also do that (it's an older version), but the real advantage of `query2haplotype.py` becomes apparent when the option `--haplo` is used. Often, the actual BLAST hit gets fragmented in pieces, in the case of, for example, internal repeats. If you have the BLAST table, you can use `BLASTstitcher.py` to fix some of that. However, you probably want the whole thing in a continuous fasta sequence. In that case `--haplo` fuses all the hits that are too close together (maximum distance controlled by `--vicinity`) into a single "haplotype". The BLAST hits used to create the haplotype can be controlled with `--evalue` and `--minsize`. 

Notice that if the query fasta file contains several sequences, these will be considered as possible 5' or 3' end of the haplotype. This is very convenient when you want a particular locus of your favorite genome and you know the flanking genes. You can take the sequence of the left (5') and right (3') genes of your locus, put them in a fasta file and use it as input query.

It's important not to confuse `--vicinity` with `--minhaplo`. The latter simply filters out the output sequences that are smaller than the value given to `--minhaplo` (default 0 bp). That is useful if you only want the big haplotypes.

You can also ask for extra bp to both sides of the haplotype with `--extrabp`. This is particularly useful when you are doing manual curation of transposable elements and you are looking for the edges of the elements.

    $ python query2haplotype.py path/to/assembly.fasta path/to/queryfile.fasta --haplo --extrabp 2000 > hits.fasta

That will add 2000 bp to the right and left of the retrieved hits.

For additional options see:

    $ python query2haplotype.py -h


*I might make `--haplo` the default in the future.*

---

For filtering of BLAST tables, you can use `BLAST_tabularparser.py`. As input it expects the file produced with the BLAST table format 6 (outfmt=6).

---

Although very useful, `query2haplotype.py` and `query2hitseq.py` only work with nucleotide sequences as input. I made `genometblastn.py` to do tBLASTn searches using a **protein sequence as query**. It does NOT make haplotypes at the moment.

---

## Note on BLAST behavior with missing data

So say that you have a sequence with `?` symbols on it. For example this toy 50bp long sequence:
    
    >toy
    ????GAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGATGAGG

If you would BLAST the sequence to itself, one would expect to have a 100% identity hit that covers the whole sequence, right? Wait, why would you want to do that, you ask? Sometimes you want to BLAST a query to a collection of sequences and see if your query is found there exactly. 

So let's use my script `query2haplotype.py`!

    $ python query2haplotype.py toy.fa toy.fa
    FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 1-4
    >toy_1-46
    ????GAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGAT

Wait, the sequence is shorter, why?! And it complained about positions 1-4. If we look at the output BLAST table:

    $ cat toyVStoy-hits.tab
    toy toy 100.000 46  0   0   1   46  1   46  7.00e-23    84.2

Suddenly the coordinates go from 1 to 46, even tho the original sequence should go from 1 to 50. If the `?` symbols are ignored, then the coordinates should still be from 5 to 50. So probably the sequence gets `?` symbols removed when making the database, and the coordinates no longer match your input fastas! My script `query2haplotype.py`, however, is still working with the original coordinates, and so the output is wrong.

What happens then with `N` symbols? Imagine that the toy sequence is now:

    $ cat toy.fa
    >toy
    NNNNGAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGATGAGG

    $ rm -r toy_db # Otherwise it will use the old sequence
    $ python query2haplotype.py toy.fa toy.fa
    >toy_1-46
    NNNNGAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGAT

Again, the sequenced is chopped. How is the table looking?
    
    $ cat toyVStoy-hits.tab
    toy toy 100.000 46  0   0   5   50  1   46  7.70e-23    84.2

Now, the coordiantes of the query are correct, but **not the coordinates of the subject**!! Hence, when my script is slicing the original subject sequence, the output is again incorrect.

What if the track of `N`s is in the middle?

    $ cat toy.fa
    >toy
    GAGCTGGATGAGCTGGATGANNNNGGAGCTGGATGAGGAGCTGGATGAGG 
    $ rm -r toy_db
    $ python query2haplotype.py toy.fa toy.fa
    >toy_1-50
    GAGCTGGATGAGCTGGATGANNNNGGAGCTGGATGAGGAGCTGGATGAGG
    $ cat toyVStoy-hits.tab
    toy toy 100.000 50  0   0   1   50  1   50  1.26e-20    77.0

That looks much better, at least the `N`s don't seem to be breaking the hit, at least at this small size. Maybe bigger tracks would be an issue.
    
    $ cat toy.fa
    >toy
    GAGCTGGATGAGCTGGATGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGAGCTGGATGAGGAGCTGGATGAGG
    $ rm -r toy_db
    $ /Users/sandralorena/Dropbox/Programming/GitHubRepos/Genomics/BLAST/query2haplotype.py toy.fa toy.fa
    >toy_75-100
    GGAGCTGGATGAGGAGCTGGATGAGG
    >toy_1-20
    GAGCTGGATGAGCTGGATGA
    $ cat toyVStoy-hits.tab
    toy toy 100.000 26  0   0   75  100 75  100 2.67e-11    48.2
    toy toy 100.000 20  0   0   1   20  1   20  4.82e-08    37.4

Indeed! It got broken, but using `-H` that get's fixed.
    
    $ /Users/sandralorena/Dropbox/Programming/GitHubRepos/Genomics/BLAST/query2haplotype.py toy.fa toy.fa -H
    >toy_1-100
    GAGCTGGATGAGCTGGATGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    NNNNNNNNNNNNNNGGAGCTGGATGAGGAGCTGGATGAGG

Finally, what if you are using other IUPAC symbols, instead of just `N`.

    $ cat toy.fa
    >toy
    YRYRGAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGATGAGG

    $ rm -r toy_db
    $ python query2haplotype.py toy.fa toy.fa
    >toy_5-50
    GAGCTGGATGAGCTGGATGAGGAGCTGGATGAGGAGCTGGATGAGG

    $ cat toyVStoy-hits.tab
    toy toy 100.000 46  0   0   5   50  5   50  8.47e-23    84.2

Finally the coordinates match! But I guess here the ambigous bases were not good to be counted as a match. What if the bases are in the middle of the sequence?

    $ cat toy.fa
    >toy
    GAGCTGGATGAGCTGGATGAGGYRYRAGCTGGATGAGGAGCTGGATGAGG
    $ rm -r toy_db
    $ python query2haplotype.py toy.fa toy.fa
    >toy_1-50
    GAGCTGGATGAGCTGGATGAGGYRYRAGCTGGATGAGGAGCTGGATGAGG

This time it got the full sequence! Probably because it's a short small track of ambigous bases in between identical segments, so it still works. However, if I put the ambigous sites at the end, they get chopped (try it for yourself!).

In conclusion, do not use sequences with `?` and `N`. Ambigous bases are also dangerous if they are on the start or the end of the sequences.

(This was tested in BLAST 2.9.0+ and BLAST 2.10.1+)
