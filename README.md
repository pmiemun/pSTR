# pSTR
Protein Short Tandem Repeats (pSTR) standalone

Protein Short Tandem Repeats (pSTR) are units directly adjacent to a conserved equal unit. When they are detected in a comparison between two protein sequences, they are checked for the variation of their unit number. pSTR is a code to look for the protein Short Tandem Repeats (pSTR) of a protein dataset. In the following example, units in square brackets are taken as pSTR with unit number variation. They are fully conserved N- or C-terminally to themselves.

AACD[ACD]EFGH -- IIKL[KL]KL\n
AACD --- EFGH[II]IIKL -- KL\n
**** --- **** -- **** -- **\n

pSTR are considered from length 2 amino acids.

As input, one must use a protein dataset in FASTA format. From the protein dataset, the code aligns pairwise (with MUSCLE v3.8.1551) all sequences, and calculates the pSTR with unit number variation (epSTR). The epSTR are mapped back to a query sequence (the one provided as input, or the first one of the provided dataset), for comparison purposes.
If no name is provided for the execution, a random number will be selected.

USAGE
	perl pSTR.pl FILE NAME

EXAMPLE
	perl pSTR.pl example1.fasta CIRBP

DEPENDENCIES
	Perl version >= v5
	BioPerl version >= 1.7.6
		Bio::SeqIO package. Can be installed using the cpan repository.
	R version >= 3.6.1
	R packages (& dependencies):
		ggplot2 >= v.3.2.1

OUTPUTS
They will be stored in folder ./NAME/, and are:
	NAME_pSTR.txt -> pSTR found per protein
	NAME_raw.txt -> epSTR found per protein pair, and mapped position in the query
	NAME_epSTR.txt -> Unique epSTR found from all protein pairs, and mapped position in the query
	NAME_results.pdf -> Heatmap with all protein pairs, showing pairs with >=1 epSTR (green) or no epSTR (yellow)
	NAME_epSTRdistribution.pdf -> Distribution of epSTR in the query sequence
