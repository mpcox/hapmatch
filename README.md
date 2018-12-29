# hapmatch

This C++ program calculates haplotypes with given ingroup and outgroup frequencies in a set of coalescent genealogies.  The code was used to calculate lineage frequencies for a bottlenecked Madagascar ingroup versus its source Indonesian outgroup in:

Cox MP, MG Nelson, MK Tumonggor, FX Ricaut and H Sudoyo. 2012. [A small cohort of Island Southeast Asian women founded Madagascar](https://doi.org/10.1098/rspb.2012.0012). *Proceedings of the Royal Society B* 279: 2761-2768.

*hapmatch* reads input from [Richard Hudson's](http://home.uchicago.edu/~rhudson1/) [*ms*](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html), like the following worked example file, [*test.ms*](test.ms):

```
ms 10 1 –t 4 –I 2 4 6 ...
...
1    111001100000110110000011110000000011101101110001110011 ⎤
2    110001100000110110000011110000000011101101110001110011   Malagasy
3    000001000000100110000011000000010000001010001100111000 
4    000001000000100110000011000000010000001010001100111000 ⎦
5    000001000000100110000011000000010000000000001100111000 ⎤
6    000001000000100110000011000000010000001010001100111000 
7    000000000010100110000011000110010000001110001100000000   Indonesians
8    000000000010100110000011000110010000001110001100000000 
9    000000000010100110000011000110000000001110001000000000 
10   000000000010100110000011000110011000001110011100000000 ⎦
```

Where sequences 1–4 are simulated data representing the Madagascar ingroup and 5–10 represent the Indonesian outgroup.

Running *hapmatch* with an ingroup frequency cutoff of ≥0.2, an outgroup frequency cutoff of ≤0.8 and requiring exact matches:

```
cat test.ms | hapmatch 4 0.2 6 0.8 0
```

returns:

```
S       Sin     Sout    hapmatch
54      31      20      1
```

This output occurs because:

* Ingroup sequence 1 has >2 mismatches with all six outgroup sequences.
* Ingroup sequence 2 has >2 mismatches with all six outgroup sequences.
* Ingroup sequence 3 has 2 mismatches with outgroup sequence 5. This meets the frequency cutoff criteria (1/4 ≥ 0.2 for the ingroup and 1/6 ≤ 0.8 for the outgroup), but not the exact match criterion.
* Ingroup sequence 4 is an exact match with outgroup sequence 6.  This meets the frequency cutoff criteria (1/4 ≥ 0.2 for the ingroup and 1/6 ≤ 0.8 for the outgroup), including the exact match criterion.  This is the match identified under the 'hapmatch' heading.


INSTALLATION

*hapmatch* requires a working installation of [Kevin Thornton](http://www.molpopgen.org/markdown/krthornt)'s [*libsequence*](https://molpopgen.github.io/libsequence/) library.  The current code has been updated to require *libsequence* v2.

The easiest way to install *libsequence* is via [bioconda](https://bioconda.github.io): 

```
conda install -c bioconda libsequence
```

The *hapmatch* source code can then be compiled with:

```
g++ -o hapmatch hapmatch.cc -lsequence -std=c++11 -Wall
```


USAGE

Usage information can be found by running the command:

```
hapmatch
```

The program expects five input values: the sample size and frequency cutoff of the ingroup, the sample size and frequency cutoff of the outgroup, and a threshold for the number of allowed mismatches.  A threshold of 0 returns exact matches only.

```
hapmatch 4 0.444 6 0.948 0
```


EXAMPLE

The following command simulates 3 datasets, each containing 10 chromosome copies, with the first 4 copies deriving from population 1 (the ingroup) and the remaining 6 copies from population 2 (the outgroup).  *hapmatch* is set to find exact matches only, where haplotypes also meet the frequency cutoffs of ≥44.4% for population 1 and ≤94.8% for population 2.

```
ms 10 3 -t 4 -I 2 4 6 1 | hapmatch 4 0.444 6 0.948 0
```

The coalescent genealogies generated by *ms* are random, but the output formatting looks like the following, with each line containing summary values for a single input dataset:

```
S    Sin    Sout    hapmatch
27    26    24    0
17    12    15    1
14    8    9    0
```
where *S* is the total number of segregating sites, *S<sub>in</sub>* is the number of segregating sites in the ingroup, *S<sub>out</sub>* is the number of segregating sites in the outgroup, and *hapmatch* is the number of haplotypes that meet the required criteria.

