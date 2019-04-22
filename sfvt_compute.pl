#!/usr/bin/perl

use strict;
use DBI;
use Cwd;
use Cwd 'chdir';
use Getopt::Std;
use Data::Dumper;

use vars qw($opt_i $opt_a $opt_f $opt_o $opt_s $opt_h);
getopts('i:a:f:o:s:h');

my $usage = "
     DESCRIPTION:
	This program compute sfvt protein position on the target sequence for flu.

     REQUIRED ARGUMENTS:
	-i input file
	-a alignment file name
	-f sequence features file name
	-o output file name
	-s sequence id

     OPTIONS:	
	-h help

";


if ($opt_h) { die $usage; }
#print "output file:$opt_o\n";

my $inputFile = $opt_i;
my $alignFile = $opt_a;
my $featureFile = $opt_f;
my $outputfile = $opt_o;
my $seq_id = $opt_s;
#print "input file: $inputFile\nalignment file: $alignFile\nfeature file: $featureFile\noutput file: $outputfile\nseq_id: $seq_id\n";

open( RESULT, ">>$outputfile" ) or die();

#get ref seq
my $ref_seq = &getRefseq;
#print "ref_seq:$ref_seq\n";

#read feature file
my @features = &readFeatures;

#run alignment
my @aln_seqs = &RunMSA;

my $aligned_ts = $aln_seqs[1];
my $aligned_rs = $aln_seqs[2];
my %pos_map_ref = &map_refseqpos_to_alignedrefseqpos($aligned_rs);
my %pos_map_target = &map_alignedtargetpos_to_targetpos($aligned_ts);


foreach my $feature (@features) {
	my ($fid, $marker) = split( /%/, $feature );
	my @all_markers = split( /,/, $marker );
	my $variant = "";
	my $displayedVariant = "";
	my $reference = "";
	my $proteinPos   = "";
	my $marker_out_flag = 0;

	foreach my $single_marker (@all_markers) {
		$single_marker =~ s/\s+//ig;
		my @pos = split( /-/, $single_marker );
		$marker_out_flag = &marker_out_of_range(\%pos_map_ref, \@pos);
		if($marker_out_flag){
			print RESULT "$seq_id\t$fid\tmarker out of range\tERR\tERR\n";
			last;
		}
		if ( $#pos > 0 ) {
			my $position = $pos[0];
			my $length   = $pos[1] - $pos[0] + 1;
			$reference = $reference . "," . &Print_seq( $position, (substr $ref_seq, $position-1, $length) );
		}
		else {
			my $position = $pos[0];
			my $length   = 1;
			$reference = $reference . "," .$position.(substr $ref_seq, $position-1, 1);
		}
					
		my $result = &get_protpos_and_alignpos(\%pos_map_ref, \%pos_map_target, $aligned_ts, \@pos);
		if(defined $result && $result ne "N/A"){
			my @tmp = split("&",$result);
			my $prot_pos = $tmp[0];
			my $align_pos = $tmp[1];
			my $formattedseq = &format_seq($aligned_rs, $aligned_ts, $align_pos);
			$proteinPos = $proteinPos . "," . $prot_pos;
			$variant = $variant . $formattedseq;
			my $formattedDisplayedVariant = &format_display_variant(\%pos_map_ref, $aligned_ts, \@pos);
			$displayedVariant = $displayedVariant . $formattedDisplayedVariant;
		}
		else{
			#variant can not be determined
			$proteinPos = $proteinPos . "," . "N/A";
			$variant = $variant . "unknown";
		}
					
	} #foreach marker

	$reference =~ s/^,//;
	$proteinPos =~ s/^,//;

	#if variant contains "unknown", then the proteinPos must be "N/A"
	if(($variant =~ m/unknown/g) ||($variant =~ m/X/gi)) {
		$variant = "unknown";
		$proteinPos = "N/A";
		$displayedVariant = &format_display_unknown_variant(\@all_markers);
	}
	$displayedVariant =~ s/,$//;

	my $tmp = $variant;
	if($tmp ne ""){
		$tmp =~ s/-//ig;
		if($tmp eq ""){
			$variant = $reference;
			$variant =~ s/,/_/ig;
		}
	}
	
	
	print RESULT "$seq_id\t$fid\t$variant\t$proteinPos\t$displayedVariant\n";

		
} #foreach my $feature

#close RESULT;

exit(0);

#read input file, @input_seqs[1] is target sequence, @input_seqs[2] is ref sequence
sub getRefseq {
	open( IN, "<$inputFile" ) or die "cannot open input file\n";
	my @input_seqs;
	my $tmp = "";
	while (<IN>) {
		chomp;
		s/\n//g;
		if (/>/) {
			push( @input_seqs, $tmp );
			$tmp = "";
		}
		else {
			$tmp = $tmp . $_;
		}
	}
	push( @input_seqs, $tmp );
	close IN;
	return $input_seqs[2];
}

sub readFeatures {
	my (@features );
	open( IN, "<$featureFile" ) or die "cannot open feature file\n";
	while (<IN>) {
		chomp;
		my ($sfn_seq_id, $pos) = split /\t/;
		push( @features, "$sfn_seq_id%$pos" );
	}
	close IN;
	return @features;
}

sub RunMSA {
	my $aa_file   = $inputFile;
	my $fasta_out = $alignFile;
#	my $clw_out   = "flu_sfvt_temp/output.aln";
#	my $cmd       = 
#	  "muscle -in $aa_file -fastaout $fasta_out -clwout $clw_out -quiet ";

	my $cmd       =
	  "clustalw -infile=$aa_file -type=protein -align -output=fasta -outfile=$fasta_out";
	my $re = &run_command( $cmd, 0 );
	if ( $re != 0 ) {
		print "\nclustalw failed.\n\n";
		print RESULT "$seq_id\t0\tclustalw failed\tERR\tERR\n";
		exit(0);
	}
	open( IN, "<$fasta_out" ) or die "cannot open aln file\n";
	my @aln_seqs;
	my $tmp = "";
	while (<IN>) {
		chomp;
		s/\n//g;
		if (/>/) {
			push( @aln_seqs, $tmp );
			$tmp = "";
		}
		else {
			$tmp = $tmp . $_;
		}
	}
	push( @aln_seqs, $tmp );
	close IN;
	return @aln_seqs;
}

sub run_command {
	my ( $cmd, $fatal ) = @_;

	if ($cmd) {
		my $re = system($cmd);

		if ( $fatal && ( $re != 0 ) ) {
			print MAIN_LOG "Fatal error: cannot excute cmd: $cmd\n";
			print "\nFatal error: cannot excute cmd: $cmd\n";
			exit(-1);
		}
		elsif ( !$fatal && ( $re != 0 ) ) {
			print MAIN_LOG "Waring: cannot excute cmd: $cmd\n";
			print "\nWarning: cannot excute cmd: $cmd\n";
			return (-1);
		}
	}
	return 0;
}

## map positions on refseq to their corresponding positions on aligned refseq
sub map_refseqpos_to_alignedrefseqpos {
	my ($aligned_refseq) = @_;
	my %map = ();
	my $ar_len = length($aligned_refseq);
	my @ar = split("", $aligned_refseq);
	my $ref_pos = 1;
	for(my $i=0; $i<$ar_len; $i++){
		if($ar[$i] ne "-"){
			$map{$ref_pos} = $i+1;
			$ref_pos++;
		}
	}
	
	return %map;
}

## map positions on aligned target sequence to their corresponding positions on target sequence
sub map_alignedtargetpos_to_targetpos {
	my ($aligned_targetseq) = @_;
	my %map = ();
	my $at_len = length($aligned_targetseq);
	my @at = split("", $aligned_targetseq);
	my $target_pos = 1;
	for(my $i=0; $i<$at_len; $i++){
		if($at[$i] ne "-"){
			$map{$i+1} = $target_pos;
			$target_pos++;
		}
		else{
			$map{$i+1} = "na";
		}
	}
	
	return %map;
}

## check if the given marker in the refseq out of range 
sub marker_out_of_range {

	my (%pos_map_ref) = %{$_[0]};
	my @pos = @{$_[1]};
	
	if( $#pos > 0 ){
		my $start = $pos[0];
		my $end = $pos[1];
		if( (!exists $pos_map_ref{$start}) || (!exists $pos_map_ref{$end}) ){
			return 1;
		}
	}
	else{
		my $position = $pos[0];
		if(!exists $pos_map_ref{$position}){
			return 1;
		}
	}

	return 0;
}

#returns prot pos on the target and the aligned pos, & separated
#if the variant cannot be determined, then only return "N/A"
#if position(s) mapped to only dash(s) in the aligned target, then return "N/A&aligned pos"
sub get_protpos_and_alignpos {
	my (%pos_map_ref) = %{$_[0]};
	my (%pos_map_target) = %{$_[1]};
	my ($aligned_ts) = $_[2];
	my @pos = @{$_[3]};

	my ($start, $end, $position);
	my ($prot_pos_start, $prot_pos_end, $prot_pos);
	my ($align_start, $align_end, $align_pos);
	my $retval;
	
	#if the given positions on the refseq is in the format of a range i.e. start-end
	if( $#pos > 0 ){
		$start = $pos[0];
		$end = $pos[1];

		$align_start = $pos_map_ref{$start};
		$align_end = $pos_map_ref{$end};

		#check if the variant can be determined
		my $tmp = substr($aligned_ts, 0, $align_start);
		my $tmp2 = substr($aligned_ts, $align_end-1);
#		print "tmp: $tmp; tmp2: $tmp2\n";
		$tmp =~ s/^-*//ig;
		$tmp2 =~ s/^-*//ig;
		if($tmp eq "" || $tmp2 eq ""){
			return "N/A";
		}

		$prot_pos_start = $pos_map_target{$align_start};
		for(my $i=$align_start; $i<$align_end; $i++){
			if($prot_pos_start ne "na"){
				last;
			}
			else{
				$prot_pos_start = $pos_map_target{$i+1};
			}
		}
		$prot_pos_end = $pos_map_target{$align_end};		
		for(my $i=$align_end; $i>$align_start; $i--){
			if($prot_pos_end ne "na"){
				last;
			}
			else{
				$prot_pos_end = $pos_map_target{$i-1};
			}
		}
		if($prot_pos_start eq "na" || $prot_pos_end eq "na"){
			$prot_pos = "N/A";
		}
		else{
			$prot_pos = $prot_pos_start . "-" . $prot_pos_end;
		}

		$retval = $prot_pos . "&" . "$align_start-$align_end";
			
	} #if pos in the format of start-end
	else{
		$position = $pos[0];

		$align_pos = $pos_map_ref{$position};

		#check if the variant can be determined
		my $tmp = substr($aligned_ts, 0, $align_pos);
		my $tmp2 = substr($aligned_ts, $align_pos-1);
#		print "tmp: $tmp; tmp2: $tmp2\n";
		$tmp =~ s/^-*//ig;
		$tmp2 =~ s/^-*//ig;
		if($tmp eq "" || $tmp2 eq ""){
			return "N/A";
		}

		$prot_pos = $pos_map_target{$align_pos};
		if($prot_pos eq "na"){
			$prot_pos = "N/A";
		}

		$retval = $prot_pos . "&" . "$align_pos";
		
	} #else single pos

	return $retval;
	
}

sub Print_seq {
	my ( $pos, $ori ) = @_;
	my $seq = "";
	my @AAs = split( //, $ori);
	for ( my $i = 0 ; $i <= $#AAs; $i++ ) {
		$seq = $seq . ( $pos + $i ) . $AAs[$i] ."_";
	}
	$seq =~ s/_$//;
	return $seq;
}

sub format_seq {
	my ($aligned_rs, $aligned_ts, $align_pos) = @_;
	my ($aligned_feat_rs, $aligned_feat_ts, $feat_len);
	my @pos = split("-", $align_pos);
	my $retval = "";

	if( $#pos > 0 ){
		$feat_len = $pos[1]-$pos[0]+1;
	}
	else{
		$feat_len = 1;
	}

	$aligned_feat_rs = substr($aligned_rs, $pos[0]-1, $feat_len);
	$aligned_feat_ts = substr($aligned_ts, $pos[0]-1, $feat_len);
#	print "feat_rs: $aligned_feat_rs\nfeat_ts: $aligned_feat_ts\n";
		
	my @rs = split("", $aligned_feat_rs);
	my @ts = split("", $aligned_feat_ts);
	for(my $i=0; $i<scalar(@rs); $i++){
		if($ts[$i] eq $rs[$i]){
			$retval = $retval . "-"; 
		}
		else{
			if($ts[$i] eq "-"){
				$retval = $retval . "?";
			}
			if($rs[$i] eq "-"){
				$retval = $retval . "[" . $ts[$i] . "]";
			}
			if($rs[$i] ne "-" && $ts[$i] ne "-"){
				$retval = $retval . $ts[$i];
			}
		}
	} 
	
	return $retval;
}

sub format_display_variant {
	my (%pos_map_ref) = %{$_[0]};
	my ($aligned_ts) = $_[1];
	my @pos = @{$_[2]};

	my $retval = "";
	my ($feat_len, $start, $align_pos, $residue);

	$start = $pos[0];
	if( $#pos > 0 ){
		$feat_len = $pos[1]-$pos[0]+1;
	}
	else{
		$feat_len = 1;
	}
	for(my $i=0; $i<$feat_len; $i++){
		$align_pos = $pos_map_ref{$start};
		$residue = substr($aligned_ts, $align_pos-1, 1);
		$retval = $retval . $residue . $start . ",";
		$start = $start + 1;
	}

	return $retval;
} #sub ends

sub format_display_unknown_variant {

	my @all_markers = @{$_[0]};
	my ($feat_len, $start, $residue);
	my $retval = "";

	foreach my $single_marker (@all_markers) {
		$single_marker =~ s/\s+//ig;
		my @pos = split( /-/, $single_marker );
		$start = $pos[0];
		if( $#pos > 0 ){
			$feat_len = $pos[1]-$pos[0]+1;
		}
		else{
			$feat_len = 1;
		}
		for(my $i=0; $i<$feat_len; $i++){
			$residue = "?";
			$retval = $retval . $residue . $start . ",";
			$start = $start + 1;
		}
		
	} #foreach


	return $retval;
} #sub ends

