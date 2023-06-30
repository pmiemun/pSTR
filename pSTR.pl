#!/usr/bin/perl
use strict;
use Bio::SeqIO;  
my $list = "./db/mammals_proteomes.txt";	my $rlist = ""; my %seqs = '';	my %mapped = ''; my $file2 = ''; my $qq = ''; my $queryid = '';
###############################################################################################################################################
################################################# 1) Get & check input file and project name ##################################################
###############################################################################################################################################
my ($input,$group) = @ARGV;
if (length($input) < 2) {	
	print "Execution halted. At least one protein in FASTA format must be supplied\n\nUsage: perl pSTR.pl FILE NAME\n"; exit; 
}
if (length($group) < 1) {	$group = int(rand(10000000));	}
my $folder = "./$group"; `mkdir $folder`; 
my $file1 = "$folder/orthologs.txt"; my $output = "$folder/repeat1.txt";
my $check = `grep ">" $input | wc -l`; chomp $check;
if		($check < 1)  { print "Execution halted. At least one protein in FASTA format must be supplied\n"; exit; }
elsif ($check == 1) {
	######################################################### 1a) Look for orthologs ############################################################
	open (OUT,">>$file1");	print OUT "Ortholog\tInput\n";	
	$file2 = "$folder/orthologs.fasta";	open (OUT2,">>$file2");	open (IN,"<$list");	my $index = -1;		print "Looking for orthologs...\n";
	while (<IN>) {		chomp $_; $index++; next if ($index == 0);		
		my @info = split(/\t/,$_);	my $proteomeid = $info[0]; my $taxid = $info[1];
		my $proteome = "./db/$proteomeid\_$taxid.fasta";		my $multi = Bio::SeqIO->new( -file  => "<$input" , '-format' => 'fasta');
		while (my $seqq = $multi->next_seq()) {	$queryid = $seqq->display_id; $qq = $seqq->seq;	} 									#Save query's information
		my $id = `blastp -query $input -db $proteome -outfmt 6 -max_target_seqs 5 | head -1 | cut -f2`;	chomp $id; 	#Search for the ortholog
		if (length($id) < 2) 	{	print OUT "-\t$proteomeid\_$taxid.fasta\n";			} 																	#Check the result
		else {		
			my $xx = `blastdbcmd -db ./db/$proteomeid\_$taxid.fasta -dbtype prot -entry '$id'`; chomp $xx;						#Get the ortholog sequence
			print OUT2 "$xx\n"; 	print OUT "$id\t$proteomeid\_$taxid.fasta\n";		
		}	
	}	close IN; close OUT; close OUT2;
	my $multi = Bio::SeqIO->new( -file  => "<$file2" , '-format' => 'fasta');																			#Store orthologs seqs in memory
	while (my $seqq = $multi->next_seq()) {	my $na = $seqq->display_id;	my $sss = $seqq->seq; $seqs{$na} = $sss; 	$rlist .= "\"$na\"\,";	}
#
print "Searching for pSTR...\n";
&indiv_pstr($file2,$folder,$group);	#Calculate individual pSTR of orthologs independently
#	
} else {
	######################################################## 1b) Store given orthologs ##########################################################
	open (OUT,">>$file1");	print OUT "Ortholog\n";
	my $multi = Bio::SeqIO->new( -file  => "<$input" , '-format' => 'fasta');
	while (my $seqq = $multi->next_seq()) {		
		my $na = $seqq->display_id;	my $sss = $seqq->seq; $seqs{$na} = $sss; print OUT "$na\n";	$rlist .= "\"$na\"\,";	
	} 	#Store orthologs seqs in memory
	close OUT;
#
print "Searching for pSTR...\n";
&indiv_pstr($input,$folder,$group);	#Calculate individual pSTR of orthologs independently
#
}
############################################# 1c) Align the query with the rest of the sequences ##############################################
my %qqseq = '';	my %ssseq = ''; my $ii = 0;	
my $new = '';	if ($check == 1) {$new = $file2;}	else {	$new = $input;}	my $xx1 = Bio::SeqIO->new( -file  => "<$new" , '-format' => 'fasta');
print "Aligning sequences...\n";
while (my $se = $xx1->next_seq()) {		$ii++;
	if (($ii == 1) && ((length($qq) == 0) || (length($queryid) == 0))) {	
		$queryid = $se->display_id; $qq = $se->seq; $qqseq{$queryid} = $qq;	$ssseq{$queryid} = $qq;
	} 
	else {
		my $seq = $se->seq; my $id = $se->display_id;
		open (OMN,">>aux5.txt"); print OMN ">$queryid\n$qq\n>$id\n$seq\n"; close OMN;	my @aabb = `muscle -in aux5.txt -quiet`;	unlink "aux5.txt";
		my $number = $#aabb + 1; my $gre = ">"; my $y = 0; my $qqqq = ''; my $sstr = ''; 
		for (my $x = 1; $x <= $number; $x++) {	chomp $aabb[$x];							#Retrieve aligned seqs
			if (index($aabb[$x], $gre) != -1) {	$y++;	} else {	if ($y == 0) {	$qqqq .= $aabb[$x];	} else {	$sstr .= $aabb[$x];	}	}
		}
		$qqseq{$id} = $qqqq;	$ssseq{$id} = $sstr;
	}
}
###############################################################################################################################################
###################################################### 2) Calculate interesting regions #######################################################
###############################################################################################################################################
chop($rlist);	open (OUT2,">>$output");	print OUT2 "AC1\tAC2\tValue\tepSTR\n";	open (IN,"<$file1");	my $index = -1;	
my $position = "$folder/site_position.txt"; 	open (POS,">>$position");		print POS "epSTR\tposition\n"; close POS;
print "Searching for epSTR...\n";
while (<IN>) {		chomp $_; $index++;  next if ($index == 0);
	my @info = split(/\t/,$_);		my $query = $info[0];
	open (IN2,"<$file1");	my $index2 = -1;		
	while (<IN2>) {		chomp $_; $index2++; next if ($index2 == 0); 
		my @info2 = split(/\t/,$_);		my $subject = $info2[0];
		if ($index2 <= $index) {print OUT2 "$query\t$subject\t-\tNA\n"; } else {
			if (($query eq "-") || ($subject eq "-")) {	print OUT2 "$query\t$subject\tNo ortholog\t-\n"; }
			else {		my $aux = "aux.txt"; open (AUX,">>$aux");	
				print AUX ">$query\n$seqs{$query}\n>$subject\n$seqs{$subject}"; close AUX; #Save both seqs in a file
				#Alignment is done with MUSCLE by default; MUSCLE is 3x faster than ClustalO; ClustalO is 4x faster than MAFFT
				my @abc = `muscle -in aux.txt -quiet`; 	#my @abc = `clustalo -i aux.txt`;	#my @abc = `mafft --quiet aux.txt`;	#Align
				my $number = $#abc + 1;	my $queryseq = ''; my $subjectseq = '';	my $gre = ">"; my $y = 0;
				for (my $x = 1; $x <= $number; $x++) {	chomp $abc[$x];							#Retrieve aligned seqs
					if (index($abc[$x], $gre) != -1) {	$y++;		} else {	if ($y == 0) { $queryseq .= $abc[$x];	} else {	$subjectseq .= $abc[$x];	}	}
				}
				my $interesting = '';		my $length_query = length($queryseq);		my $length_subject = length($subjectseq);
				#Check gaps in query vs what's aligned in subject		//	#Check gaps in subject vs what's aligned in query				
				$interesting = &sites($length_query,$queryseq,$subjectseq,$interesting,$position,$qqseq{$subject},$ssseq{$subject});
				$interesting = &sites($length_subject,$subjectseq,$queryseq,$interesting,$position,$qqseq{$query},$ssseq{$query});
				if (length($interesting) < 2) { print OUT2 "$query\t$subject\tNo result\t+\n"; } 
				else {	print OUT2 "$query\t$subject\tResult\t$interesting\n";	}	#Print results
	unlink "aux.txt";		}	}	}			
} close IN; close OUT2;
###############################################################################################################################################
############################################################## 3) Plot results ################################################################
###############################################################################################################################################
my $ouou = "$folder/$group\_results.pdf"; 					my $title = "All epSTR"; &plot($rlist,$output,$ouou,$title);
my $sssi = "$folder/$group\_epSTRdistribution.pdf"; 	my $pstr = "$folder/$group\_epSTR.txt";
`sort -n $position | uniq | sort -k 2 -n >> $pstr`;	unlink $position;		#Unique pSTR
open(AUX3,">>R_distrib");	print AUX3 "library(ggplot2)\n";	print AUX3 "df2 = read.table('$pstr',sep='\t',header=TRUE,dec='.')\n";
print AUX3 "d<-ggplot(df2, aes(y=epSTR, x=position)) + scale_x_continuous(limits = c(0, NA)) + \n";
print AUX3 "labs(x=\"Mapped position in the query ($queryid)\",y=\"epSTR\") + \n";
print AUX3 "theme_light() + geom_point(shape=16) + theme(legend.position=\"none\")\n";	print AUX3 "pdf(\'$sssi\');d;dev.off()\n";	
close AUX3;	`Rscript R_distrib`;	unlink "R_distrib"; 
my $ou2 = "$folder/repeat_evolutio.txt";	open (OOA,">>$ou2");	open (IN,"<$output"); 	##Filter out unnecessary lines in output file
while (<IN>) {	chomp $_;	my @info = split(/\t/,$_); next if (($info[2] eq "\-") && ($info[3] eq "NA"));print OOA "$_\n";} close IN; close OOA;

`grep -v "+" $ou2 | awk \'{ print \$1,\$2,\$4 }\' >> $folder/$group\_raw.txt`; unlink $ou2;
unlink $output; if ($check > 1) { unlink $file1; }
print "\nYour results are in folder \"$folder/\"\n";
exit;
###############################################################################################################################################
########################################################## SUBROUTINE: Plot results ###########################################################
###############################################################################################################################################
sub plot () {
	my ($rlist,$inp,$output,$title) = @_;
	open(AUX3,">>R_sites");		print AUX3 "library(ggplot2)\n";	print AUX3 "dat = read.table('$inp',sep='\t',header=TRUE)\n";
	print AUX3 "AC3 <- factor(dat\$AC1, levels = c($rlist))\nAC4 <- factor(dat\$AC2, levels = c($rlist))\n";
	print AUX3 "d1 <- ggplot(dat, aes(x = AC3, y = AC4, fill = Value)) + geom_tile(color = \"gray\")+  coord_fixed() +\n"; 
	print AUX3 "theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) + \n";
	print AUX3 "theme(text = element_text(size = 6.5), axis.title.x=element_blank(), axis.title.y=element_blank()) + \n";
	print AUX3 "labs(title=\"$title\")  + scale_fill_manual(name=\"\", values = c(\"\#ff8080ff\",\"\#ffd42aff\",\"\#88aa00ff\")) +";
	print AUX3 "theme(legend.position = \"none\") + scale_x_discrete(position = \"top\")\n";
	print AUX3 "pdf(\'$output\');d1;dev.off()\n";
	close AUX3;	`Rscript R_sites`; unlink "R_sites";
}
###############################################################################################################################################
################################################# SUBROUTINE: Detect pSTR using protein pairs #################################################
###############################################################################################################################################
sub sites () {				#Look for pSTR
	my ($length,$seq1,$seq2,$interesting,$position,$exec_query_seq,$subjectseq) = @_;		open (POS,">>$position");
	for (my $a = 0; $a <= $length; $a++) {									#Check gap regions
		my $residue = substr($seq1,$a,1);
		if ($residue eq "-") {		my $rr = "-";	
			my $start = $a; while ($rr eq "-") {	$a++;	$rr = substr($seq1,$a,1);	}	my $end = $a-1;	#Subject has a gap region from $start to $end
			my $reg = substr($seq2,$start,$end-$start+1);				#Extract the region in the query aligning with the gap region	
			
			next if (length($reg) == 1);
							
			my $lele = $end-$start+1; 													#Look for same region right before or right after in the subject & in the query
			my $regA = substr($seq1,$start-$lele,$lele);		my $regB = substr($seq1,$end+1,$lele);			#			A---B
			my $regC = substr($seq2,$start-$lele,$lele);		my $regD = substr($seq2,$end+1,$lele);			#			CxxxD
			my $init = $start-$lele; if ($init < 0) 		{ $regA = "-"; $regC = "+";		}									#If upstream region is out of the prot
			my $fin = $end+1+$lele; if ($fin > $length) { $regB = "-"; $regD = "+";		}									#If downstream region is out of the prot
			if ((($reg eq $regA) && ($reg eq $regC)) || (($reg eq $regB) && ($reg eq $regD)))	{ 
			#Now that a pSTR has been found, we calculate its position in the protein that contains it
				my $aaal = substr($seq2,0,$start);	#Retrieve all protein from the beginning to the start of the pSTR
				my $count1 = $aaal =~ tr/-//;	my $count2 = $seq2 =~ tr/-//;		#Count gaps in region & in complete protein
				my $middle = int((($end-$count1+1)+($start-$count1+1))/2);		#pSTR position is the middle between the start and end, minus the gaps
			#Map the position of the pSTR to the alignment of the subject with the execution query
				my $lesu = length($subjectseq); my $mapped = ''; my $real = 0;
				for (my $x = 0; $x <= $lesu; $x++) {							
					my $rere = substr($subjectseq,$x,1);	if ($rere ne "-") { $real++; }	if ($real eq $middle) { $mapped = $x; $x = $lesu+1; } 
				}
			#Map the position of the pSTR to the execution query itself
				my $gap = 0; my $pos_exec_query = '';
				for (my $z = 0; $z <= $lesu; $z++) {
					my $erer = substr($exec_query_seq,$z,1); if ($erer eq "-") {	$gap++; }
					if ($z == $mapped) {	$pos_exec_query = $mapped-$gap+2; $z = $lesu+1;	print POS "[$reg]\t$pos_exec_query\n";	}
				}
			$interesting .= "[$reg($pos_exec_query)]";					#Save results	
			}	
		}	
	}	close POS;	return($interesting); 
}
###############################################################################################################################################
##################################################### SUBROUTINE: Detect pSTR per protein #####################################################
###############################################################################################################################################
sub indiv_pstr () {				#Look for pSTR
	my ($file,$folder,$group) = @_;
	my $output = "$folder/$group\_pST.txt"; open (OUT,">>$output");	print OUT "AC\tpSTR\tUnits\tStart\tEnd\n";
	my $multi = Bio::SeqIO->new( -file  => "<$file" , '-format' => 'fasta');
	while (my $seqq = $multi->next_seq()) {		
		my $na = $seqq->display_id;	my $sss = $seqq->seq;	my $len = length($sss);
		for (my $a = 0; $a <= $len; $a++) {		#Go over the full sequence, taking each aa as the start of pSTRs of all possible size
			my $max_pstr = ($len-$a)/2;		#Theoretically, the longest pSTR has a length of half of the protein
	
			for (my $b = 2; $b <= $max_pstr; $b++) {	#Check all possible lengths for pSTR starting from the current aa, from length = 2 aa
	
				my $rere = substr($sss,$a,$b);		#From current aa ($a), get the possible unit of the pSTR, of length $b
				my $units = 0; my $lll = $a+$b;
				my $rere2 = $rere;
				while ($rere eq $rere2) {	$rere2 = substr($sss,$lll,$b);	$units++; $lll = $lll+$b;	}
				my $start = $a + 1; my $end = $lll-$b;
				if ($units > 1) {	print OUT "$na\t$rere\t$units\t$start\t$end\n";	}
			}
		}
	} close OUT;
	open (OUT2, ">>$folder/$group\_pSTR.txt");	open (IN,"<$output");	my $indd = 0;	my $ya = '';
	while (<IN>) { chomp $_; $indd++;
		if ($indd == 1) {	print OUT2 "$_\n";	} 
		else {
			my @info = split(/\t/,$_);	my $actual = "($info[0]|$info[1]|$info[4])";	#AC, Unit, and End position
			next if (index($ya, $actual) != -1);			#Next if we already have a pSTR with the same unit from this protein ending in same position
			print OUT2 "$_\n";	$ya .= $actual;
		}
	}	close IN;	unlink $output;
}
