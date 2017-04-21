#!/usr/bin/perl -w
use strict;
#use lib '/home/jwang2/AssembleTest/Detect/nonSJsrc/dev';
use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Sam;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use File::Basename;
use Cwd;
use List::MoreUtils qw/ uniq /;
#custom packages
use SCValidator qw($lowqual_cutoff $min_percent_id $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use SVUtil qw($use_scratch $work_dir get_input_bam get_work_dir parse_range is_PCR_dup);

use constant FQ_BASE_NUMBER => 33;

my $rmdup = 0;
my $print_read = 1;
$use_scratch = 0;
$min_percent_id = 90;
my ($work_dir, $out_dir);
my $ref_genome;
# input/output
my ($out_prefix, $range, $input_bam );
my $out_suffix = "sclip.txt";
my $paired = 1;
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'o|out_dir=s'	=> \$out_dir,
	'ref_genome=s'  => \$ref_genome,
	'p=s'			=> \$out_prefix,
	'scratch!'		=> \$use_scratch,
	'paired!'		=> \$paired,
	'rmdup!'		=> \$rmdup,
	'lq_cutoff=i'	=> \$lowqual_cutoff,
	'min_pct_id=i'	=> \$min_percent_id,
	'min_pct_hq=i'	=> \$min_percent_hq,
	'print_read!'	=> \$print_read,
	'r|range=s'		=> \$range,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

my $start_dir = getcwd;
if($input_bam) {
	croak "The bam file you specified does not exist!\n" unless(-e $input_bam);
	$input_bam = File::Spec->rel2abs($input_bam);
}
else{
	croak "You must specify the input bam file or sample name";
}
croak "You must provide the reference genome in fasta format!" if(!$ref_genome);
croak "The reference genome file you speicified does not exist!\n"
	unless(-e $ref_genome);

my $input_base = fileparse($input_bam);

#setup output dir and workind directory
$out_dir = getcwd if(!$out_dir);
mkdir $out_dir if(!-e $out_dir || ! -d $out_dir);
$work_dir = get_work_dir() if($use_scratch);
$work_dir = $out_dir if(!$work_dir);
$use_scratch = undef if($work_dir eq $out_dir);
chdir($work_dir); 

# figure out output prefix
$out_prefix = $input_base if(!$out_prefix);

my $sam = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);

my $output_file;
my $validator = SCValidator->new();
if(!$paired) {
	$validator->remove_validator("strand_validator");
}
if($range) {
	my ($chr, $start, $end) = parse_range($range);
	my $tmp = $chr;
	$tmp = $tmp . ".$start" if($start);
	$tmp = $tmp . ".$end" if($end);
	$output_file = join('.', $out_prefix, $tmp, $out_suffix);
	$output_file = File::Spec->catfile($out_dir, $output_file);

	my($pcover, $ncover) = extract_range_sclip(
		-SAM => $sam, 
		-RANGE => $range, 
		-WORK_DIR => $work_dir, 
		-OUTPUT => $output_file,  
		-VALIDATOR => $validator);
	$output_file = join('.', $out_prefix, $tmp, "cover");
	$output_file = File::Spec->catfile($out_dir, $output_file);
	warn "$output_file file exists and it will be replaced!" if(-e $output_file);
	$rmdup = 0;
	#$start = 0 unless $start;
    #my ($tid) = $sam->header->parse_region($chr);
    #my $chr_len = $sam->header->target_len->[$tid];
	#$end = $chr_len unless $end;
	
	#my ($coverage) = $sam->features(-type => 'coverage', -seq_id => $chr,
	#	-start => $start, -end => $end);
	#my @cc = $coverage->coverage;

	$rmdup = 0;
	my($fh, $fname) = tempfile( DIR => $work_dir);
	foreach my $p (keys(%{$pcover})) {
		my $c = ($rmdup ? scalar(@{$pcover->{$p}}) : $pcover->{$p});
		print $fh join("\t", $chr, $p, "+", $c, count_coverage($sam, $chr, $p)), "\n";
	}
	foreach my $p (keys(%{$ncover})) {
		my $c = ($rmdup ? scalar(@{$ncover->{$p}}) : $ncover->{$p});
		print $fh join("\t", $chr, $p, "-", $c , count_coverage($sam, $chr, $p)), "\n";
	}
	system("sort -k 2 -n $fname -o $output_file");
	system("rm $fname");
}
else{
	$output_file = join('.', $out_prefix, $out_suffix);
	$output_file = File::Spec->catfile($out_dir, $output_file);
    my $header = $sam->header;
    my $target_names = $header->target_name;
    $output_file = join('.', $out_prefix, "cover");
	my $read_file = join('.', $out_prefix, "sclip.txt");
    $output_file = File::Spec->catfile($out_dir, $output_file);
    $read_file = File::Spec->catfile($out_dir, $read_file);

	if(-e $output_file) {
		warn "$output_file file exists and it will be replaced!";
		system("rm $output_file");
	}
	my @t_names = uniq @{$target_names};
    foreach my $chr (@t_names){
		my($fh, $fname) = tempfile( DIR => $work_dir);
        my($pcover, $ncover) = extract_range_sclip(
			-SAM => $sam,
			-RANGE =>$chr, 
			-WORK_DIR => $work_dir, 
			-OUTPUT => $read_file, 
			-VALIDATOR => $validator);
		foreach my $p (keys(%{$pcover})) {
			my $c = ($rmdup ? scalar(@{$pcover->{$p}}) : $pcover->{$p});
			print $fh join("\t", $chr, $p, "+", $c, count_coverage($sam, $chr, $p) ), "\n";
		}
		foreach my $p (keys(%{$ncover})) {
			my $c = ($rmdup ? scalar(@{$ncover->{$p}}) : $ncover->{$p});
			print $fh join("\t", $chr, $p, "-", $c, count_coverage($sam, $chr, $p)), "\n";
		}
		system("sort -k 2 -n $fname -o $fname.sorted");
		system("cat $fname.sorted >> $output_file");
		system("rm $fname");
		system("rm $fname.sorted");
    }
}
chdir $start_dir;
exit(0);

sub count_coverage {
	my ($sam, $chr, $pos, $clip) = @_;
	if($rmdup) {
		my @pairs;
		my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
		my $n = 0;
		my $itr = $seg->features(-iterator => 1);
		while( my $a = $itr->next_seq) {
			my $sclip_len = 0;
			if($clip) {
				my @cigar_array = @{$a->cigar_array};
				$sclip_len = $cigar_array[0]->[1] if($cigar_array[0]->[0] eq 'S' && $clip == RIGHT_CLIP);
				$sclip_len = $cigar_array[$#cigar_array]->[1] if($cigar_array[$#cigar_array]->[0] eq 'S' && $clip == RIGHT_CLIP);
			}
			next if(@pairs > 0 && is_PCR_dup($a, \@pairs, $sclip_len));
			$n++;
			#return $n if( $n > $max_repetitive_cover);
			if($a->mpos) {
				push @pairs, [$a->start, $a->end, $a->mate_start, $a->mate_end, $sclip_len];
			}
			else {
				push @pairs, [$a->start, $a->end, 0, 0, $sclip_len];
			}

		}
		return $n;
	}
	else{
		my ($c) = $sam->features(-type => 'coverage', -seq_id=> $chr, 
			-start => $pos, -end => $pos);
		my @c_d = $c->coverage;
		return $c_d[0];
	}
}

sub extract_range_sclip {
	my %arg = @_;
	my $sam = $arg{-SAM} || croak "missing -SAM";
	my $range = $arg{-RANGE} || croak "missing -RANGE";
	my $output_file = $arg{-OUTPUT} || croak "missing -OUTPUT";
	my $validator = $arg{-VALIDATOR} || croak "missing -VALIDATOR";

	my($fh, $fname) = tempfile( DIR => $work_dir);
	my (%plus_cover, %neg_cover);
	$sam->fetch($range, 
		sub {
			my $a = shift;
			my $cigar_str = $a->cigar_str;
#			print STDERR $a->qname, "\t", $cigar_str, "\n";
			my @cigar_array = @{$a->cigar_array};
			return if($a->cigar_str !~ m/S/);
			if($paired && !$a->proper_pair) { #paired but mate is not mapped 
				$validator->remove_validator("strand_validator");
			}
			#return if(!$a->proper_pair && $paired);
			#return if($paired && !$a->mpos); 
			my ($sclip_len, $ort, $pos, $seq, $qual_str, $qual);
			$qual_str = join( "",  (map { chr $_ + FQ_BASE_NUMBER } $a->qscore));
			if($cigar_array[0]->[0] eq 'S' && $validator->validate($a, LEFT_CLIP) ) {
				$sclip_len = $cigar_array[0]->[1]; $ort = "-"; $pos = $a->start;
				$seq = substr($a->query->dna, 0, $sclip_len );
				$qual = substr($qual_str, 0, $sclip_len);
				
				my $print = 1;
				if($rmdup) {
					if(exists $neg_cover{$pos}) {
						$print = 0 if(is_PCR_dup($a, $neg_cover{$pos}, $sclip_len));
					}
					else {
						$neg_cover{$pos} = [];
					}
					if($print == 1) {
						if($a->mpos) {
							push @{$neg_cover{$pos}}, [$a->start, $a->end, $a->mate_start, $a->mate_end, $sclip_len];
						}
						else {
							push @{$neg_cover{$pos}}, [$a->start, $a->end, 0, 0, $sclip_len];
						}
					}
				}
				else {
					if(exists $neg_cover{$pos}) {
						$neg_cover{$pos}++;
					}
					else{
						$neg_cover{$pos} = 1;
					}
				}
				print $fh join("\t", $a->seq_id, $pos, $ort, $a->qname, $seq, $qual), "\n"
					if($print_read && $print == 1);
			}

			if($cigar_array[$#cigar_array]->[0] eq 'S' &&  $validator->validate($a, RIGHT_CLIP)) {
				$sclip_len = $cigar_array[$#cigar_array]->[1]; $ort = '+'; $pos = $a->end;
				my $l = length($a->query->dna);
				$seq = substr($a->query->dna, $l - $sclip_len);
				$qual = substr($qual_str, $l - $sclip_len);

				my $print = 1;
				if($rmdup) {
					if(exists $plus_cover{$pos}) {
						$print = 0 if(is_PCR_dup($a, $plus_cover{$pos}, $sclip_len));
					}
					else {
						$plus_cover{$pos} = [];
					}
					if($print ==1 ) {
						if($a->mpos) {
							push @{$plus_cover{$pos}}, [$a->start, $a->end, $a->mate_start, $a->mate_end, $sclip_len];
						}
						else {
							push @{$plus_cover{$pos}}, [$a->start, 0, 0, $a->mate_end, $sclip_len];
						}
					}
				}
				else {
					if(exists $plus_cover{$pos}) {
						$plus_cover{$pos}++;
					}
					else{
						$plus_cover{$pos} = 1;
					}
				}
				$print = is_PCR_dup($a, $plus_cover{$pos}, $sclip_len) if($rmdup);
				print $fh join("\t", $a->seq_id, $pos, $ort, $a->qname, $seq, $qual), "\n"
					if($print_read && $print);
			}
			if($paired && !$a->proper_pair) { #paired but mate is not mapped, add back
				$validator->add_validator("strand_validator"); 
			}
		}
	);
	if($print_read) {
		system("sort -k 2 -n $fname -o $fname.sorted");
		system("cat $fname.sorted >> $output_file");
		system("rm $fname");
		system("rm $fname.sorted");
	}
	return(\%plus_cover, \%neg_cover);
}


=head1 NAME

extractSClip.pl - extract positions with soft clipped read in bam file.


=head1 VERSION

This documentation refers to extractSClip.pl version 0.0.1.


=head1 USAGE

	# extract all positions with soft clipped reads in whole genome:
	./extractSClip.pl -i sample.bam -g hg18.fa 
	# extract chr1 positions with soft clipped reads
	./extractSClip.pl -i sample.bam -g hg18.fa -r chr1


=head1 REQUIRED ARGUMENTS

	-i: Input bam file.
	--ref_genome: The genome file in fa file, must be the same used to map reads.
	

=head1 OPTIONS

	-r: The range of positions need to be extracted. Format: chr1 or chr1:500-5000.
	-o: The output directory, default is current directory.
	--scratch: use scracth space, default off.
	--rmdup: remove PCR dumplicate reads, default on, use --normdup to turn it off.
	--lq_cutoff: low quality cutoff value, default 20. 
	--min_pct_id: minimum percent identify for the aligned high qual part,default 90.
	--min_pct_hq: minimum percent high quality for soft clipped part, default 80.
	--print_read: individual soft-clipped read will be printed, default off.
	-h, --help: The help page.
	--man: Print the man page.
	--usage: Print usage information.
	-v, --version: print version information.


=head1 DESCRIPTION

This is a program to extract all soft-clipped positions such that for each position
a list of requirements need to be satisfied. More specifically, the orientaion of 
pair-end read should be satisfied, the minimum percent identify of aligned part 
need to be satisfied, and the minimum percent of hiqh quality soft-clipped part 
should be satisfied.


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
