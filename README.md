# pSTR
Protein Short Tandem Repeats (pSTR) standalone code. The web tool associated with this code is available in http://cbdm-01.zdv.uni-mainz.de/~munoz/pSTR/.

Protein Short Tandem Repeats (pSTR) are units directly adjacent to a conserved equal unit. When they are detected in a comparison between two protein sequences, they are checked for the variation of their unit number. pSTR is a code to look for the protein Short Tandem Repeats (pSTR) of a protein dataset. 

pSTR are considered from length 2 amino acids.

As input, one must use a protein dataset in FASTA format. From the protein dataset, the code aligns pairwise (with MUSCLE v3.8.1551) all sequences, and calculates the pSTR with unit number variation (epSTR). The epSTR are mapped back to a query sequence (the one provided as input, or the first one of the provided dataset), for comparison purposes.
If no name is provided for the execution, a random number will be selected.

USAGE -> perl pSTR.pl FILE NAME

EXAMPLE -> perl pSTR.pl example1.fasta CIRBP

DEPENDENCIES
1) Perl version >= v5
2) BioPerl version >= 1.7.6
3) Bio::SeqIO package. Can be installed using the cpan repository.
4) R version >= 3.6.1
5) R package ggplot2 >= v.3.2.1

OUTPUTS -> They will be stored in folder ./NAME/, and are:
1) NAME_pSTR.txt -> pSTR found per protein
2) NAME_raw.txt -> epSTR found per protein pair, and mapped position in the query
3) NAME_epSTR.txt -> Unique epSTR found from all protein pairs, and mapped position in the query
4) NAME_results.pdf -> Heatmap with all protein pairs, showing pairs with >=1 epSTR (green) or no epSTR (yellow)
5) NAME_epSTRdistribution.pdf -> Distribution of epSTR in the query sequence
