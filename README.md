# Linear-time detection of repetitive motifs in DNA

This is a project I did for my MSc in Stochastics and Data Science at the University of Turin. Everything is in one Jupyter notebook: `SimoneMariaGervasoniCodeKMP.ipynb`.

The idea is simple. Some diseases are caused by the same short piece of DNA repeating too many times in a row. The standard example is Huntington's disease, where the trinucleotide CAG repeats inside the huntingtin gene on chromosome 4: up to about 26 repeats is normal, past 35–40 the disease shows up, and in between there is an intermediate range where the person is healthy but the repeat can expand in their children. So the practical problem is: given a chromosome of roughly 190 million letters, find every occurrence of a pattern and measure how long the runs of consecutive repeats are. Doing that with a naive scan is too slow, which is where the algorithms come in.

## What the notebook does

First I implement the two matchers from scratch and check them against each other:

- a brute force matcher, which tries the pattern at every position;
- the Knuth–Morris–Pratt algorithm, split into the failure function (`error_fun`) that preprocesses the pattern and the actual search (`KMP`).

On the same test string they return exactly the same positions, which is the sanity check for everything that follows. There is also a short Manim animation of the KMP scan on a toy example, because the way the two indices move is much easier to see than to read.

Then the notebook uses KMP on real data. A small function (`pattern_detector`) takes the match positions and works out the lengths of the runs of consecutive repeats. Three cases are analysed:

- a plasmid taken from Addgene that carries part of the mutated gene, used in place of a full "ill patient" genome (full genomes of patients are not public, for obvious privacy reasons);
- chromosome 4 from the T2T-CHM13v2.0 assembly, as the intermediate case: the longest run of consecutive CAG found is 28 repeats, which sits in the intermediate range;
- chromosome 4 from the GRCh38 reference (Ensembl release 110), as the healthy case: the longest run is 19 repeats, inside the normal range.

The full chromosome scan takes under a minute.

After that I implement Aho–Corasick, which searches for several patterns in a single pass over the text using a trie with failure links. Running it on the healthy chromosome with the patterns TTAGGG, CAG and CCG gives the same counts reported in Table 4 of the paper *Genomic Repeat Detection Using the Knuth-Morris-Pratt Algorithm on R High-Performance-Computing Package*, which was a second check that the implementations are right.

The last part is a speed comparison. Brute force, KMP and Aho–Corasick are timed on slices of chromosome 4 from 10 to about 16 million base pairs, with a long pattern of 200 characters. With a single pattern KMP is the fastest of the three. The comparison is then repeated with many patterns at once, which is the situation Aho–Corasick is built for: one pass instead of one full scan per pattern.

## Files

- `SimoneMariaGervasoniCodeKMP.ipynb` — the whole project: implementations, animations, the chromosome scans and the benchmarks.
- `SimoneMariaGervasoniReportKMP.pdf` — the written report, with the complexity analysis (KMP searches in linear time in the size of the text plus the pattern; the failure function is linear in the pattern) and the details on the datasets.

## Notes on running it

The notebook downloads the two chromosome files itself (from Ensembl and NCBI) into a `data` folder on the first run; together they are about 190 million base pairs each once extracted, so give it a minute and some disk space. It needs Biopython to parse the FASTA files and Manim only for the animation cells — the rest is plain Python with matplotlib and numpy.
