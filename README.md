This script finds the most common chromosomal segments
shared identical-by-descent (IBD). It was used to find common 
segments for the paper by Szydlowski: 'A clue to the etiology of
disorders of sex development from identity-by-descent analysis
in dogs with cryptic relatedness' published in Animal Genetics.

To search for most common chromosomal segment (of IBD2 type:
both alleles shared identical-by-descent), the algorithm
considers all subsets of individuals, starting from the largest
subsets of size n-1, n-2, n-3 etc. For example, in the full set
of n=11 individuals, the algorithm searches for common segments 
shared by all 11 individuals and then generates 11 subsets of 10
individuals, 55 subsets of 9 individuals, 165 subsets of 8 individuals,
and so on. To find chromosomal segments shared identical-by-descent
in a set of individuals, the shared fragments (of IBD2 type) are first
detected in all pairwise comparisons using KING (version 2.2.7).

You need to have BEDTOOLS installed on your computer.
https://bedtools.readthedocs.io/en/latest/

Before using the script, please run KING (version 2.2.7). 
You need the output file 'king.segments'.
Reference: Manichaikul et al. 2010, Bioinformatics 26: 2867-2873.

You need to specify IDs of individuals to consider.  

The output file 'cibd.out.segments.bed' gives:
Chrom, Start, End, Segment size, IDS of individuals sharing the segment (joined with '&').

The output file 'cibd.out.sets.txt' gives:
Segment size, Subset size, IDs in the subset.
