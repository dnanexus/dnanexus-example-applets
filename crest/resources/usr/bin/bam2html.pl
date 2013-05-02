#!/usr/bin/perl -w
use strict;
use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use Data::Dumper;
use Bio::Assembly::IO;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use File::Spec;
use File::Path;
use File::Copy;
use File::Basename;

use constant HTML_COLOR => {
	'MATCH'		=> '#000000', #black
	'MISMATCH'	=> '#800000', #brown
	'LOWQUAL'	=> '#C0C0C0', #grey
	'INDEL'		=> '#FF8040', #Orange
	'SCLIP'		=> '#6960EC', #Slate Blue2
};
my $lowqual_cutoff = 20;

my ($bam_d, $bam_g,  $file);
my $ref_genome;
my $range;

my ( $help, $man, $version, $usage );
my $output;
my $RNASeq;
my $optionOK = GetOptions(
	'd=s'			=> \$bam_d,
	'g=s'			=> \$bam_g,
	'f|i|file|input=s'	=> \$file,
	'r|ref_genome=s'	=> \$ref_genome,
	'range=s'		=> \$range,
	'RNASeq'		=> \$RNASeq,
	'o|output=s'	=> \$output,
	'h|help|?'      => \$help,
	'man'           => \$man,
	'usage'         => \$usage,
	'v|version'     => \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );
croak "you must specify at least one bam file" unless ($bam_d || $bam_g);
warn "Missing reference genome, the view is not correct" unless($ref_genome);
croak "You need either provide an input file or a range" unless ($file || $range);

my @bams;
push @bams, Bio::DB::Sam->new( -bam => $bam_d, -fasta => $ref_genome) if($bam_d);
push @bams, Bio::DB::Sam->new( -bam => $bam_g, -fasta => $ref_genome) if($bam_g);
open STDOUT, ">$output" if $output;

if($range) {
	print STDOUT bam2html($range, @bams);
	close STDOUT if $output;
	exit(0);
}

open my $FILE, "<$file" or croak "can't open $file:$OS_ERROR";
while( my $line = <$FILE>) {
	chomp $line;
	my @fields = split /\t/, $line;
	print "<h3 style=\"color:blue;text-align:center\">$line</h3>\n";
	if(scalar(@fields) == 3) { # deal with file: chr	start	end format
		my ($chr, $s, $e) = @fields;
		print STDOUT bam2html("$chr:$s-$e", @bams);
		next;
	}
	# deal with CREST output file
#	if($fields[0] =~ /^chr/) {
#		$fields[0] = substr($fields[0], 3);
#		$fields[2] = substr($fields[2], 3);
#	}
	my ($chr, $s, $e) = ($fields[0], $fields[1], $fields[1]);
	print STDOUT bam2html("$chr:$s-$e", @bams);
	($chr, $s, $e) = ($fields[4], $fields[5], $fields[5]);
	print STDOUT bam2html("$chr:$s-$e", @bams);
}
close STDOUT if($output);
exit(0);
#print bam2html("12:11676858-11676858", $bam1, $bam2);


sub get_input_bam {
    my ($raw_bam_dir, $sample, $group) = @_;

    $raw_bam_dir = File::Spec->catdir($raw_bam_dir, $group);
    opendir(my $dh, $raw_bam_dir) or croak "can't open directory: $raw_bam_dir: $OS_ERROR";
    my @files = grep { /^$sample-.*bam$/ } readdir($dh);
	close $dh;
	return $files[0];
}


# this function returns a <pre> </pre> block for the specific contig in the ace file
sub bam2html {
	my ($range, @bams) = @_;
	my ($chr, $start, $end) = $range =~ m/^(.*?):(\d+)-(\d+)/;
	if(!$chr or !$start or !$end) {
		croak "The range format is not correct, use chr:start-end";
	}
	my ($r_start, $r_end) = ($start, $end);
	my $name_len = length($range);
	for my $bam (@bams) {
		my $segment = $bam->segment(-seq_id => $chr, -start => $start, -end => $end);
		my @alignments = $segment->features;
	#	return "There are too many reads in this region, this is must be a highly repetitive region, abort!"
	#		if(scalar(@alignments) > 500);
		for my $a(@alignments) {
			my $l = length($a->query->name);
			$name_len = $l if($name_len < $l);
		}
	}
	if($start == $end) { # we want to check a specific point
		if($RNASeq) { # for RNASeq we just extend a little bit instead of dynamic do it
			$r_start = $start - 100;
			$r_end = $end + 100;
		}
		else {	
			for my $bam (@bams) {
				my $segment = $bam->segment(-seq_id => $chr, -start => $start - 1, -end => $end + 1);
				my @alignments = $segment->features;
				for my $a (@alignments) {
					next if($a->start > $end || $a->end < $start);
					$r_start = $a->start if $a->start < $r_start;
					$r_end = $a->end if $a->end > $r_end;
				}
			}
		}
		$range = "$chr:$r_start-$r_end";
	}
	my $rtn_str = "<pre>\n";
	my($ref, $pos2padded) = get_padded_ref($range, @bams);
	$rtn_str .= print_bam_ruler($r_start, $r_end, $pos2padded, $name_len) . "\n";
	my $line = print_ref($ref, $range, $name_len);
	my $line_len = length($line);
	$rtn_str .= $line . "\n" . ' ' x $line_len . "\n";
	$ref =~ s/\*//g; #remove * from ref sequence
	for my $bam (@bams) {
		my $segment = $bam->segment(-seq_id => $chr, -start => $start, -end => $end);
		my @alignments = $segment->features;
		for my $a (@alignments) {
			#next if($a->cigar_str !~ m/S/);
			if($a->strand == 1) {
				my $tmp_str = print_bam_seq($a, $r_start, $r_end, $ref, $pos2padded, $name_len);
				$rtn_str .= $tmp_str . "\n" if($tmp_str); 
			}
		}
		for my $a (@alignments) {
			#next if($a->cigar_str !~ m/S/);
			if($a->strand == -1) {
				my $tmp_str = print_bam_seq($a, $r_start, $r_end, $ref, $pos2padded, $name_len);
				$rtn_str .= $tmp_str . "\n" if($tmp_str);
			}
		}
		$rtn_str .= '_' x $line_len . "\n";
	}
	$rtn_str .= "\n</pre>";
	return $rtn_str;
}

# only print the unpadded position, which is meaningful
sub print_bam_ruler {
	my($start, $end, $pos2padded, $name_len) = @_;
	my ($string, $mark);
	$string = $mark = " "x($name_len + 2);
	
	#print | at 20 and . at 10 
	for( my $i = $start; $i <= $end; $i++) {
		$mark .= $i % 10 == 0 ? ($i % 20 == 0 ? "|" : ".") : " ";
		last if($i == $end);
		$mark .= ' ' x ($pos2padded->[$i - $start + 2] - $pos2padded->[$i - $start + 1] - 1);
	}
	# print numbers above | for padded consensus
	my $padded_i = 0;
	for( my $i = 1; $i <= $end - $start + 1; $i++ ) {
		my $j = $i + $start - 1;
		my $overhead = $pos2padded->[$i] - $padded_i;
		if($overhead > 0 ) {
			$string .= ' ' x $overhead;
			$padded_i += $overhead;
		}
		if( $j % 20 == 0) {
			my $l = length($j);
			my $half_l = int(($l+1)/2);
			$string = substr($string, 0, length($string) - $half_l );
			$string .= $j;
			$padded_i += ($l - $half_l);
		}
	}
	return join("\n", ($string, $mark));
}

sub print_ref {
	my( $ref, $chr, $name_len) = @_;
	return $chr . " " x ($name_len + 2 - length($chr)) . $ref;
}

sub get_padded_ref {
	my ($range, @bams) = @_;
	my $ref_str; #padded ref genome
	my @pos2padded;
	my($chr, $start, $end) = $range =~ m/(.*?):(\d+)-(\d+)/;
	
	@pos2padded = 0 .. ($end - $start + 1);
	push @pos2padded, 100 + $end - $start;
	
	foreach my $bam (@bams) {
		my ($d, $cum_pad) = (0, 0);
		my $padded_fun = sub {
			my ($seqid, $pos, $p) = @_;
			return if($pos < $start || $pos > $end);
			$d++; 
			$pos2padded[$d] += $cum_pad;
			my $max_ins = 0;
			for my $pileup (@$p) {
				$max_ins = $pileup->indel if($max_ins < $pileup->indel);
			}
			if($max_ins > $pos2padded[$d+1] + $cum_pad - $pos2padded[$d] - 1) {
				$cum_pad += $max_ins - ($pos2padded[$d+1] + $cum_pad-  $pos2padded[$d] - 1);
			}
		};
		$bam->pileup($range, $padded_fun);
		if($d <= $end - $start + 1) {
			for(my $i = $d + 1; $i <= $end-$start+1; $i++) {
				$pos2padded[$i] += $cum_pad;
			}
		}
	}
	my $refseq = $bams[0]->segment($chr, $start, $end)->dna;
	my @str = split //, "*" x $pos2padded[$end - $start + 1];
	for my $i ( 1 .. ($end - $start + 1)) {
		$str[$pos2padded[$i]] = substr($refseq, $i-1, 1);
	}
	shift @str; #remove the 0th base
	return (join( '', @str), \@pos2padded );
}

sub print_bam_seq { # it's complicated
	my ($align, $start, $end, $ref, $pos2padded, $name_len) = @_;
	my $rtn = $align->query->name;
	$rtn .= ' ' x ($name_len - length($rtn));
	$rtn .= $align->strand == -1 ? '- ' : '+ ';
	my $p_s = $align->start;	
	my $p_e = $align->end;
	my $seq = $align->query->dna;
	my @qual = $align->qscore;
	my $repeat = 0;
	if($align->has_tag("XT")) {
		$repeat = 1 if($align->aux_get("XT") ne "U");
	}

	my $r_cigar = $align->cigar_array;
	my ($s, $e) = ($align->query->start, $align->query->end);
	
	# reads partial in the region
	my $leading = "";
	if($p_s < $start) {
		($r_cigar, $p_s, $s) = find_new_start($r_cigar, $p_s, $s, $start);
	}

	my $tailing = "";
	if($p_e > $end) {
		($r_cigar, $p_e, $e) = find_new_end($r_cigar, $p_e, $e, $end);
	}
	return if($s >= $e);	
	my @cigar = @{$r_cigar};

	# deal with leading softclip
    my	$op = shift @cigar;
	my $l =  $pos2padded->[$p_s - $start + 1] - $pos2padded->[1];
	if($op->[0] eq 'S') {
		my $ss = $op->[1] - $l;
		if($ss < 0) {
			$ss = -$ss;
			$leading .= ' ' x $ss;
			$ss = 0;
		}
		$leading .= '<font color="' . HTML_COLOR->{'SCLIP'} . '">';
		for( my $i = $ss; $i < $op->[1]; $i++) {
			my $c = substr($seq, $i, 1);
			$c = lc $c if($qual[$i] < $lowqual_cutoff);
			$leading .= $c;
		}
		$leading .= '</font>' ;
	}
	else {
		$leading .= ' ' x $l;	
		unshift @cigar, $op;
	}

	#dealing with tailing softclip
	$op = pop @cigar;
	$l = $pos2padded->[$end - $start + 1] - $pos2padded->[$p_e - $start + 1];
	if($op->[0] eq 'S') {
		my $ee = $op->[1] - $l;
		$ee = 0 if($ee < 0 );
		$tailing = '<font color="' . HTML_COLOR->{'SCLIP'} . '">';
		for( my $i = 0; $i < $op->[1] - $ee; $i++){
			my $c = substr($seq, $i + $e, 1);
			$c = lc $c if($qual[$i + $e] < $lowqual_cutoff);
			$tailing .= $c;
		}
		$tailing .= '</font>';
	}
	else{
		push @cigar, $op;
	}

	# generate the alignment part, only M, I, D and N
	my $mid = '';
	my $mode;
	foreach $op (@cigar) {
		my $l = $op->[1];
		if($op->[0] eq 'M') {
			my $line = '';
			while($l > 0 ){
				$l--;
				my $c = substr($seq, $s - 1, 1); #seq is 0 based
				my $newmode = "MISMATCH";
				my $cc = substr($ref, $p_s - $start, 1); 
				#last unless ($c && $cc);
				$newmode = "MATCH" if($cc eq $c);
				if($qual[$s-1] < $lowqual_cutoff) {
					$c = lc $c;
					$newmode = 'LOWQUAL' if($newmode eq 'MATCH');
				}
				if(!$mode) {
					$line .= '<font color="' . HTML_COLOR->{$newmode} . '">';
				}
				elsif($mode ne $newmode) {
					$line .= '</font>' . '<font color="' . HTML_COLOR->{$newmode} . '">';
				}
				$mode = $newmode;
				$line .= $c;
				$s++; $p_s++;
				# dealing with padded * in reference genome
				if($p_s < $p_e && $s < $e && $l > 0) { 
					my $tmp = $pos2padded->[$p_s - $start + 1]-$pos2padded->[$p_s-1 - $start + 1];
					if( $tmp > 1) {	
						$line .= '</font>';
						$line .= '<font color="' . HTML_COLOR->{'INDEL'} . '">';
						$line .= '*' x ($tmp - 1);
						$mode = 'INDEL';
					}
				}
			}
			$mid .= $line;
		}
		if($op->[0] eq 'D' || $op->[0] eq 'I' || $op->[0] eq 'N') {
			my $newmode = 'INDEL';
			my $tmp; #extra padded * after indel
			$mid .= '</font>' if($mode && $mode ne $newmode);
			$mid .= '<font color="' . HTML_COLOR->{'INDEL'} . '">';

			if($op->[0] eq 'D' || $op->[0] eq 'N') {
				$mid .= ($op->[0] eq 'D' ? '*' : '=') x $l;
				$p_s += $l;
				$tmp = $pos2padded->[$p_s - $start + 1]-$pos2padded->[$p_s - $l - $start ] - $l if($p_s < $p_e);
			}
			else{
				$tmp = $pos2padded->[$p_s - $start + 1] - $pos2padded->[$p_s - 1 - $start + 1] - $l if($p_s < $p_e);
				while($l > 0 ) {
					my $c = substr($seq, $s - 1, 1);
					$c = lc $c if($qual[$s - 1] < $lowqual_cutoff);
					$mid .= $c;
					$l--; $s++;
				}
			}
			$mode = $newmode;
			if($p_s < $p_e && $tmp > 1) {	
				$mid .= '*' x ($tmp - 1);
			}
		}
	}
	$mid .= '</font>';
	$rtn .= $leading . $mid . $tailing;		
	$rtn = '<b><i>' . $rtn . '</i></b>' if($repeat);
	return $rtn;
}

sub find_new_start {
	my ($r_cigar, $p_s, $s, $start) = @_;
	my @cigar = @{$r_cigar};
 
	while(1) {
		my $op = shift @cigar;
		next if( $op->[0] eq 'S' || $op->[0] eq 'H');
		if( $op->[0] eq 'I') { 
			$s += $op->[1];	next;
		}
		if( $p_s + $op->[1] < $start ) { 
			$p_s += $op->[1];	
			$s += $op->[1] if $op->[0] eq 'M'; 
		}
		else {
			$s += ($start - $p_s) if $op->[0] eq 'M'; 
			unshift @cigar, [$op->[0], $op->[1] - ($start - $p_s)];
			$p_s = $start;
			return (\@cigar, $p_s, $s);
		}
	}
}

sub find_new_end {
	my ($r_cigar, $p_e, $e, $end) = @_;
	my @cigar = @{$r_cigar};

	while(1) {
		my $op = pop @cigar;
		next if( $op->[0] eq 'S' || $op->[0] eq 'H');
		if( $op->[0] eq 'I') { 
			$e -= $op->[1];	next;
		}
		if( $p_e - $op->[1] > $end ) { 
			$p_e -= $op->[1];	
			$e -= $op->[1] if $op->[0] eq 'M'; 
		}
		else {
			$e -= ($p_e - $end) if $op->[0] eq 'M';
			push @cigar, [$op->[0], $op->[1] - ($p_e - $end)];
			$p_e = $end;
			return (\@cigar, $p_e, $e);
		}
	}
}

=head1 NAME

bam2html.pl - a bam file viewer that just simple display part of the alignment as HTML file.


=head1 VERSION

This documentation refers to bam2html.pl version 0.0.1.


=head1 USAGE
	
	Display part of a bam file:
	    bam2html.pl -r hg18.fa -d diag.bam --range 1:123566-123766 -o diag.html
	Display part of two bam files, one diagnositc, one germlie for comparison.
	    bam2html.pl -r hg18.fa -d diag.bam -g germline.bam --range 1:123566-123766 -o diag.html
	Display part of two bam files, one diagnositc, one germlie for comparison from
	a list of positions in a file, each line should be tab sepearted as: chr, start, and end.
	    bam2html.pl -r hg18.fa -d diag.bam -g germline.bam -o diag.html -f position.txt
	Display part of two bam files, one diagnositc, one germlie for comparison from a file
	generated by CREST.pl.
	    bam2html.pl -r hg18.fa -d diag.bam -g germline.bam -o predSV.html -f predSV.txt

=head1 REQUIRED ARGUMENTS

	To run the program, several parameter must specified.
	-d:					The input (diagnositic) bam file
	-g:					The input (germ line) bam file
	-r, --ref_genome: 	The reference genome file in fa format

=head1 OPTIONS

	The	options that can be used for the program.
	-o: 		The output file, default to STDOUT if missing.
	--range:	The range where SV will be detected, using chr1:100-200 format.
	-f, -i:		The input file from either CREST.pl or tab seperated chr, start, end.

	-h, --help	 Help information
	--man		 Man page.
	--usage		 Usage information.
	--version	 Software version.


=head1 DESCRIPTION

 This is a bam file viewer that just simple display part of the alignment as 
 HTML file. The program is developed to view Structure Varitions around break 
 point. So manual review will a breeze.  Any way this program can be used
 otherway, but don't put a too big range to display as the view will not
 be pretty as each read will occupy a line.


=head1 DIAGNOSTICS

If the program does not respond for minutes, please be a little bit patient as
sometimes generating the html file could take longer time.
If you provide a range as "chr1:50-100" and nothing was output, please make sure
the reference genomes for mapping and display are exact the same and the chrom
name is chr1 not 1.


=head1 DEPENDENCIES

The program depend on several packages:
1. Bioperl perl module.
2. Bio::DB::Sam, version 1.5 or later, it requires samtools lib installed.


=head1 BUGS AND LIMITATIONS

There are no known bugs in this module, but the method is limitted to bam file 
that has soft-clipping cigar string generated.Please report problems to 
Jianmin Wang  (Jianmin.Wang@stjude.org)
Patches are welcome.

=head1 AUTHOR

Jianmin Wang (Jianmin.Wang@stjude.org)



=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010 by St. Jude Children's Research Hospital.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
