#!/usr/bin/perl -w
#use lib '/home/jwang2/AssembleTest/Detect/nonSJsrc/dev';
use strict;
use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use Bio::SearchIO;
use Bio::SeqIO;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use File::Path;
use File::Copy;
use File::Basename;
use List::MoreUtils qw/ uniq /;
use Cwd;
use SCValidator qw($lowqual_cutoff $min_percent_id $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use SVUtil qw($use_scratch $work_dir clip_fa_file prepare_reads parse_range get_input_bam 
	find_smallest_cover get_work_dir is_PCR_dup read_fa_file get_sclip_reads);
use StructVar;
require SVExtTools;

use Transcript;
use Gene;
use GeneModel;

my $debug = 1;
$min_percent_id = 90;
# input/output
my ( $out_dir, $out_prefix, $range, $bam_d,	$bam_g, $sclip_file);
my $out_suffix = "predSV.txt";
my $read_len = 100;
my $ref_genome ;

my $hetero_factor = 0.4;
my $triger_p_value = 0.05;

#RNASeq support
my $RNASeq = 0;
my $gene_model_file;
my $gm_format = "REFFLAT";
my $gm;
my $max_sclip_reads = 50; # we will consider the position if we have enough sclipped reads
my $min_one_side_reads = 5;
my $min_pct_sclip_reads = 5; # we require at least 5% reads have softclipping

# external programs related variable
# 1. cap3 related variables
my $cap3 = "cap3";
my $cap3_options=" -h 70 -y 10 > /dev/null";

# 2 blat related variables, using blat server and blat standalone
my $blat_client_exe = "gfClient";
my $blat_client_options = ' -out=psl -nohead -minIdentity=95 -maxIntron=5';  
my $target_genome = "hg18.2bit";
my $dir_2bit = "/";
my $blat_server = "sjblat";
my $blat_port = 50000;
my $blat_exe = "blat";
my $blat_options = " -tileSize=7 -stepSize=1 -out=psl -minScore=15 -noHead -maxIntron=1 ";

$use_scratch = 0;
my $paired = 1;

my $rescue = 1;
# other options 
my ($min_sclip_reads, $max_repetitive_cover, $min_sclip_len );
my $min_hit_reads;
my $rmdup;
my $min_clip_len;
my ($max_score_diff, $max_num_hits);
my ($min_dist_diff, $min_hit_len) = (15, 15);

#SV filter related parameters
my $max_bp_dist = 15;
my $min_percent_cons_of_read = 0.75;
my $germline_seq_width = 100;
my $germline_search_width = 50;
my $rm_tandem_repeat = 1; #tandem repeat mediated events
my $tr_max_indel_size = 100;   #tandem repeat mediated indel size
my $tr_min_size = 2;
my $tr_max_size = 8; 	  #tandem repeat max_size of single repeat
my $tr_min_num = 4;		  #tandem repeat minimum number of occurence

# verification pupuse
my $verify_file;
my $sensitive;
my $cluster_size;

# common help options
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	# SV validation related parameters
	'max_bp_dist=i'		=> \$max_bp_dist,
	'min_percent_cons_of_read=f'	=> \$min_percent_cons_of_read,
	'germline_seq_width=i'	=> \$germline_seq_width,
	'germline_search_width=i'	=> \$germline_search_width,
	# tandem repeat related parameters
	'rm_tandem_repeat!'	=> \$rm_tandem_repeat,
	'tr_max_indel_size=i'	=> \$tr_max_indel_size,
	'tr_min_size=i'	=> \$tr_min_size,
	'tr_max_size=i'	=> \$tr_max_size,
	'tr_min_num=i'	=> \$tr_min_num,
	'hetero_factor=f'	=> \$hetero_factor,
	'triger_p_value=f'	=> \$triger_p_value,
	
	# input/output 
	#'s|sample=s'	=> \$sample,
	'd|input_d=s'	=> \$bam_d,
	'g|inpug_g=s'	=> \$bam_g,
	'f|sclipfile=s'	=> \$sclip_file,
	'p|prefix=s'	=> \$out_prefix,
	'o|out_dir=s'	=> \$out_dir,
	'l|read_len=i'	=> \$read_len,
	'ref_genome=s'	=> \$ref_genome,
	'v|verify_file=s'	=> \$verify_file,
	'sensitive'		=> \$sensitive,	
	'paired!'		=> \$paired,
	'rescue!'		=> \$rescue,
	'cluster_size=i'	=> \$cluster_size,
	# external programs location and options
	'scratch!'		=> \$use_scratch,
	'cap3=s'		=> \$cap3,
	'cap3opt=s'		=> \$cap3_options,
	'blatclient=s'	=> \$blat_client_exe,
	'blatclientopt=s'	=> \$blat_client_options,
	'blat=s'		=> \$blat_exe,
	't|target_genome=s'	=> \$target_genome,
	'blatopt=s'		=> \$blat_options,
	'blatserver=s'	=> \$blat_server,
	'blatport=i'	=> \$blat_port,
	'2bitdir=s'		=> \$dir_2bit,
	#RNAseq support
	'RNASeq'		=> \$RNASeq,
	'genemodel=s'		=> \$gene_model_file,
	'gmformat=s'		=> \$gm_format,
	#other related parameters
	'r|range=s'		=> \$range,
	'max_score_diff=i'		=> \$max_score_diff,
	'm|min_sclip_reads=i'		=> \$min_sclip_reads,
	'min_one_side_reads=i'		=> \$min_one_side_reads,
	'max_rep_cover=i'	=> \$max_repetitive_cover,
	'min_sclip_len=i'	=> \$min_sclip_len,
	'min_hit_len=i'		=> \$min_hit_len,
	'min_dist_diff=i'	=> \$min_dist_diff,
	'min_hit_reads=i'	=> \$min_hit_reads,
	'rmdup!'		=> \$rmdup,

	# soft_clipping related parameters (from SVDector package)
	'min_percent_id=i'	=> \$min_percent_id,
	'min_percent_hq=i'	=> \$SCValidator::min_percent_hq,
	'lowqual_cutoff=i'	=> \$lowqual_cutoff,

	# common help parameters
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'version'		=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

my $start_dir = getcwd;

# figure out input file
if(!$bam_d) {
	pod2usage(1);
	croak "You need specify input bam file(s)";
}

if(!$sclip_file) {
	pod2usage(1);
	croak "You need to specify the softclipping file";
}
$sclip_file = File::Spec->rel2abs($sclip_file);

if(!$ref_genome) {
	pod2usage(1);
	croak "You need to specify the reference genome used by bam file";
}

croak "The file you sepcified does not exist" unless (
	-e $sclip_file && -e $bam_d && -e $ref_genome && -e $target_genome);
my $input_base;
$bam_d = File::Spec->rel2abs($bam_d);
$input_base = fileparse($bam_d);
$bam_g = File::Spec->rel2abs($bam_g) if($bam_g);

#RNASeq support
if($RNASeq) {
	croak "You need specify the input gene model file" unless ($gene_model_file);
	StructVar->add_RNASeq_filter();
	$min_sclip_reads = 10 unless $min_sclip_reads;
	$max_repetitive_cover = 5000 unless $max_repetitive_cover;
	$min_sclip_len = 20 unless $min_clip_len;
	$max_score_diff = 5 unless $max_score_diff;
	$max_num_hits = 3 unless $max_num_hits;
	$min_hit_reads = 5 unless $min_hit_reads;
	$cluster_size = 5 unless $cluster_size;
}
else {
	$min_sclip_reads = 3 unless $min_sclip_reads;
	$max_repetitive_cover = 500 unless $max_repetitive_cover;
	$min_sclip_len = 20 unless $min_clip_len;
	$max_score_diff = 5 unless $max_score_diff;
	$max_num_hits = 10 unless $max_num_hits;
	$min_hit_reads = 1 unless $min_hit_reads;
	$cluster_size = 1 unless $cluster_size;
}

if($gene_model_file) {
	$gm = GeneModel->new if($gene_model_file);
	$gm->from_file($gene_model_file, $gm_format);
	StructVar->gene_model($gm);
}

# set up the external programs and validators
# Those variable will be global
my $assembler = Assembler->new( 
	-PRG => $cap3,
	-OPTIONS => $cap3_options
);

my $mapper = Mapper->new(
	-PRG => join(' ', ($blat_client_exe, $blat_server, $blat_port)),
	-OPTIONS => $blat_client_options,
	-BIT2_DIR => $dir_2bit,
	-MAX_SCORE_DIFF => $max_score_diff,
	-MAX_NUM_HITS => $max_num_hits,
);

my $aligner = Aligner->new(
	-PRG => $blat_exe,
	-OPTIONS => $blat_options,
);

StructVar->assembler($assembler);
StructVar->read_len($read_len);
StructVar->aligner($aligner);
StructVar->mapper($mapper);
StructVar->RNASeq(1) if($RNASeq);
StructVar->genome($target_genome);
StructVar->remove_filter("tandem_repeat") unless($rm_tandem_repeat);
StructVar->tr_max_indel_size($tr_max_indel_size);
StructVar->tr_min_size($tr_min_size);
StructVar->tr_max_size($tr_max_size);
StructVar->tr_min_num($tr_min_num);
StructVar->max_bp_dist($max_bp_dist) if($max_bp_dist);
StructVar->germline_seq_width($germline_seq_width) if($germline_seq_width);
StructVar->germline_search_width($germline_search_width) if($germline_search_width);


my $validator = SCValidator->new();
$validator->remove_validator('strand_validator') if(!$paired);

#setup output and working directory
$out_dir = getcwd if(!$out_dir);
mkdir $out_dir if(!-e $out_dir || ! -d $out_dir);
$work_dir = get_work_dir(-SCRATCH => $use_scratch);
$work_dir = $out_dir if(!$work_dir);
$use_scratch = undef if($work_dir eq $out_dir); # don't erase the out_dir
chdir($work_dir); 

# figure out output prefix
if(!$out_prefix) {
	$out_prefix = $input_base;
}

my $sam_d = Bio::DB::Sam->new( -bam => $bam_d, -fasta => $ref_genome);
StructVar->sam_d($sam_d);
my $sam_g;
if($bam_g) {
	$sam_g = Bio::DB::Sam->new( -bam => $bam_g);
	StructVar->sam_g($sam_g);
}

#my $output_file = File::Spec->catfile($outdir, $out_prefix . $out_suffix);
my $output_file;
$output_file = join('.', $out_prefix, $out_suffix);
$output_file = File::Spec->catfile($out_dir, $output_file);

# the softclip file is sorted, so no need to re-sort it
open my $SCLIP, "<$sclip_file" or croak "can't open $sclip_file:$OS_ERROR";
my %sclip;
my $pre_p;
while( my $line = <$SCLIP> ) {
	chomp $line;
	my ($chr, $pos, $ort, $cover, $C) = split /\t/, $line;
	if( ! exists($sclip{$chr})) {
		$sclip{$chr} = [];
		$pre_p = $pos;
	}
	$C = 30 unless $C;
	if($pos < $pre_p) {
		print STDERR "The input file is not sorted!";
		exit(1);
	}
	push @{$sclip{$chr}}, [$pos, $ort, $cover, $C];
}
close $SCLIP;

my @final_SV;
my @s_cover = @{find_smallest_cover($min_sclip_reads, $hetero_factor, $triger_p_value)};
if($range) {
	@final_SV = detect_range_SVs($sam_d, $range, \%sclip, \@s_cover);
}
else {
	foreach my $chr (keys %sclip) {
		my @tmpSV = detect_range_SVs($sam_d, $chr, \%sclip, \@s_cover);
		@final_SV = @tmpSV unless (@final_SV);
		foreach my $sv (@tmpSV) {
			push @final_SV, $sv if($sv && ! is_dup_SV(\@final_SV, $sv));
		}
	}
}
undef %sclip;
open my $OUT, ">$output_file" or croak "can't open $output_file:$OS_ERROR";
foreach my $SV (@final_SV) {
	if($SV->filter) {
		print $OUT $SV->to_string ;
		if($RNASeq) {
			print $OUT "\t", join("\t", @{$SV->get_genes});
		}
		print $OUT "\n";
	}
}
chdir $start_dir;
exit(0);

sub detect_range_SVs {
	my ($sam_d, $range, $r_sclips, $r_scover) = @_;
	my ($chr, $start, $end) = parse_range($range);
	my ($tid) = $sam_d->header->parse_region($chr);
	my $chr_len = $sam_d->header->target_len->[$tid];
	my $r_range_sclips = search_sclip($r_sclips, $range);
	return unless $r_range_sclips;
	my @scover = @{$r_scover};
	my @rtn;

	my ($f_tree, $r_tree);
	if($RNASeq) {
		my $full_chr = $chr;
		my ($gm_chr) = $gm->get_all_chr;	
		$full_chr = "chr" . $chr if($chr !~ m/chr/ && $gm_chr =~ /chr/);
		($f_tree, $r_tree) = ($gm->sub_model($full_chr, "+"), $gm->sub_model($full_chr, "-"))	if($RNASeq);
	}
	my $n = scalar @{$r_range_sclips};
	push @{$r_range_sclips}, [$r_range_sclips->[$n-1][0] + $cluster_size + 1, '+', 1, 50];
	push @{$r_range_sclips}, [$r_range_sclips->[$n-1][0] + $cluster_size + 1, '-', 1, 50];
	$n = $n + 2;

	# try to save the 1 base off problem of soft-clipping 
	my (@ss1, $p1, $c1);
	my (@ss2, $p2, $c2);
	my ($n1, $n2, $cover1, $cover2) = (0, 0, 0, 0);
	for( my $i = 0; $i < $n; $i++) {
		my $s = $r_range_sclips->[$i];
#		print join("\t", @{$s}), "\n" if($debug);
		my $clip = $s->[1] eq '+' ? RIGHT_CLIP : LEFT_CLIP;
		next if($s->[0] >= $chr_len );
		if($clip == RIGHT_CLIP) {
			$p1 = $c1 = $s->[0] unless $p1;
			if($s->[0] < $c1) {
				print STDERR "The soft-clipping file has problem! It's either not sorted or
					the file has mutliple parts for a chromsome!";
				last;
			}
			if($s->[0] - $c1 < $cluster_size) {
				$n1 += $s->[2];
				$cover1 = $s->[3] if($cover1 < $s->[3]);
				push @ss1, $s; 
				if($s->[0] < $c1) {
					$p1 += $s->[0];
					$c1 = $p1 / scalar @ss1;
				}
				next;
			}
		}
		else {
			$p2 = $c2 = $s->[0] unless $p2;
			if($s->[0] < $c2) {
				print STDERR "The soft-clipping file has problem! It's either not sorted or
					the file has mutliple parts for a chromsome!";
				last;
			}
			if($s->[0] - $p2 < $cluster_size) {
				$n2 += $s->[2];
				$cover2 = $s->[3] if($cover2 < $s->[3]);
				push @ss2, $s; 
				if($s->[0] < $c2) {
					$p2 += $s->[0];
					$c2 = $p2 / scalar @ss2;
				}
				next;
			}
		}
				
		my @cc = $clip == LEFT_CLIP ? @ss2 : @ss1;	
		my $pos = $clip == LEFT_CLIP ? $ss2[$#ss2]->[0] : $ss1[0]->[0]; 
		my ($n_s, $cover_s) = $clip == LEFT_CLIP ? ($n2, $cover2) : ($n1, $cover1);
		my $tmp_range = $clip == LEFT_CLIP ? [$ss2[0]->[0], $ss2[$#ss2]->[0] + $cluster_size] : 
			[$ss1[0]->[0] - $cluster_size, $ss1[$#ss1]->[0]];

		if($clip == LEFT_CLIP) {
			@ss2 = (); $p2 = $c2 = $s->[0];
			$n2 = $s->[2];
			$cover2 = $s->[3];
			push @ss2, $s;
		}
		else {
			@ss1 = (); $p1 = $c1 = $s->[0];
			$n1 = $s->[2];
			$cover1 = $s->[3];
			push @ss1, $s;
		}

	    my @s_reads;
		if($RNASeq) {
#			print $n_s, "\n";
			next if($n_s < $min_sclip_reads && $cover_s > $s_cover[$n_s] ) ;
			next if($n_s < $max_sclip_reads && $n_s * 100 <= $cover_s * $min_pct_sclip_reads);
			my @f_genes = $f_tree->intersect($tmp_range);
			my @r_genes = $r_tree->intersect($tmp_range);
			next if(scalar @f_genes == 0 && scalar @r_genes == 0);
		}
		else {
			if($n_s < $min_sclip_reads) { #too few covered
				if($cover_s > $s_cover[$n_s]) {
					next if(!$sensitive);
					foreach my $c (@cc) {
						next if($c->[0] >= $chr_len);
						push @s_reads, get_sclip_reads(-SAM => $sam_d,
										   -CHR =>$chr, 
										   -START => $c->[0], 
										   -END => $c->[0], 
										   -CLIP => $clip, 
										   -MINCLIPLEN => 0,
										   -VALIDATOR => $validator,
										   -PAIRED => $paired,
										   -RMDUP => $rmdup,
										   -EXTRA => abs($c->[0] - $pos) );
					}
				}
			}
		}
		if(! @s_reads) {
			foreach my $c (@cc) {
				next if($c->[0] >= $chr_len);
				push @s_reads, get_sclip_reads(-SAM => $sam_d,
										   -CHR =>$chr, 
										   -START => $c->[0], 
										   -END => $c->[0], 
										   -CLIP => $clip, 
										   -MINCLIPLEN => 0,
										   -VALIDATOR => $validator,
										   -PAIRED => $paired,
										   -RMDUP => $rmdup,
										   -EXTRA => abs($c->[0] - $pos) );
			}
		}
		my $l = 0;
		foreach my $r (@s_reads) {
			my $len = length $r->{seq};
			$l = $len if($l < $len);
		}
		next if ($l < $min_sclip_len);
		print STDERR join("\t", ($chr, $pos, $clip == LEFT_CLIP ? "-" : "+", scalar @s_reads)), "\n";
		my @SV =  detect_SV(
			-SAM => $sam_d, 
			-CHR => $chr, 
			-POS => $pos, 
			-ORT => $clip == LEFT_CLIP ? "-" : "+", 
			-SCLIP => \%sclip,
			-COVER => $cover_s,
			-SREADS => \@s_reads);
		foreach my $sv (@SV) {
			$sv->{type} = $sv->type;
			push @rtn, $sv if($sv && ! is_dup_SV(\@rtn, $sv));
		}
	}

	if($RNASeq) {
		push @rtn,	find_del_SVs($sam_d, $chr, $start, $end);
	}
	return @rtn;
}

sub is_dup_SV {
	my($r_SVs, $sv) = @_;
	foreach my $s (@{$r_SVs}) {
		return 1
		if( $s->first_bp->{pos} == $sv->second_bp->{pos} &&
			$s->first_bp->{chr} eq $sv->second_bp->{chr} &&
			$s->second_bp->{pos} == $sv->first_bp->{pos} &&
			$s->second_bp->{chr} eq $sv->first_bp->{chr} &&
			$s->{type} == $sv->{type} );
		return 1
		if( $s->first_bp->{pos} == $sv->first_bp->{pos} &&
			$s->first_bp->{chr} eq $sv->first_bp->{chr} &&
			$s->second_bp->{pos} == $sv->second_bp->{pos} &&
			$s->second_bp->{chr} eq $sv->second_bp->{chr} &&
			$s->{type} == $sv->{type} );

	}
	return;
}

sub search_sclip {
	my ($r_sclip, $range) = @_;
	my ($chr, $start, $end) = parse_range($range);
	return $r_sclip->{$chr} if(!$start);
	my $start_i = bin_search($r_sclip->{$chr}, $start);	
	my $end_i = bin_search($r_sclip->{$chr}, $end);
	my @tmp = @{$r_sclip->{$chr}};
	if($start_i <= $end_i) {
		@tmp = @tmp[$start_i .. $end_i];
		if($start_i == $end_i) {
			return undef if($tmp[0]->[0] < $start || $tmp[0]->[0] > $end);
		}
		return \@tmp;
	}
	else {
		return undef;
	}
}

sub bin_search {
	my ($a, $p) = @_;
	my ($s, $e) = (0, scalar(@{$a})-1);
	my $m = int( ($s + $e)/2);
	while(1) {
		return $s if($a->[$s][0] >= $p);	
		return $e if($a->[$e][0] <= $p);
		return $m if($a->[$m][0] == $p || ($a->[$m-1][0] < $p && $a->[$m+1][0] > $p));
		if($a->[$m][0] > $p) {
			$e = $m;
		}
		else {
			$s = $m;
		}
		$m = int( ($s+$e)/2 );
	}
}

sub count_coverage {
	my ($sam, $chr, $pos, $clip) = @_;
	if($rmdup && !$RNASeq) {
		my @pairs;
		my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
		my $n = 0;
		return 0 unless $seg;
		my $itr = $seg->features(-iterator => 1);
		while( my $a = $itr->next_seq) {
			next unless($a->start && $a->end); #why unmapped reads here?
			my $sclip_len = 0;
			if($clip) {
				my @cigar_array = @{$a->cigar_array};
				#$sclip_len = $1 if($a->cigar_str =~ m/S(\d+)$/ && $clip == RIGHT_CLIP); 
				$sclip_len = $cigar_array[0]->[1] if($cigar_array[0]->[0] eq 'S' && $clip == RIGHT_CLIP);
				#$sclip_len = $1 if($a->cigar_str =~ m/^S(\d+)/ && $clip == LEFT_CLIP); 
				$sclip_len = $cigar_array[$#cigar_array]->[1] if($cigar_array[$#cigar_array]->[0] eq 'S' && $clip == RIGHT_CLIP);
			}
			next if(@pairs > 0 && is_PCR_dup($a, \@pairs, $sclip_len));
			$n++;
			return $n if( $n > $max_repetitive_cover);
			push @pairs, [$a->start, $a->end, $a->mate_start ? $a->mate_start : 0, 
				$a->mate_end ? $a->mate_end : 0, $sclip_len];
		}
		return $n;
	}
	else{
		my ($c) = $sam->features(-type => 'coverage', -seq_id=> $chr, 
			-start => $pos, -end => $pos);
		return 0 unless $c;
		my @c_d = $c->coverage;
		return $c_d[0];
	}
}

sub detect_SV {
	my %param = @_;
	my @s_reads = @{$param{-SREADS}};
	my($sam, $chr, $ort, $r_sclip, $coverage) = (	$param{-SAM}, 
		$param{-CHR},  $param{-ORT},	$param{-SCLIP},	$param{-COVER});
	my $clip = $ort eq '+' ? RIGHT_CLIP : LEFT_CLIP;
	return if($coverage > $max_repetitive_cover);
	my $fa_name = prepare_reads(\@s_reads, $clip);
	my($contig_file, $sclip_count, $contig_reads) = $assembler->run($fa_name); 
	return if(!$contig_file or !(-e $contig_file));
	
	$contig_file = clip_fa_file($contig_file, $clip);
	my $contig_seqs = read_fa_file($contig_file);	
	if(scalar @s_reads == 1) {
		$contig_file = clip_fa_file($fa_name, $clip);
		$contig_seqs = read_fa_file($contig_file);
		my @seq_names = keys %{$contig_seqs};
		$sclip_count->{$seq_names[0]} = 1;
	}

	my @SV;
	my $where = ($clip == LEFT_CLIP ? "right" : "left");
	my $mapping = $mapper->run( -QUERY => $contig_file );
	my (%reads, %quals, %s_lens, %pos);
	foreach my $r (@s_reads) {
		$reads{$r->{name}} = $r->{full_seq};
		$quals{$r->{name}} = $r->{full_qual};
		$s_lens{$r->{name}} = length($r->{seq});
		$pos{$r->{name}} = $r->{pos};
	}

	foreach my $s (keys(%{$mapping})) {
		next if($sclip_count->{$s} < $min_sclip_reads);

		# try to find a better mapping
		my @tmp_reads;
		my $selected_read;
		my @tmp_pos;
		foreach my $n (@{$contig_reads->{$s}}) {
			push @tmp_pos, $pos{$n};
		}
		@tmp_pos = uniq @tmp_pos;
		my $real_pos = $clip == LEFT_CLIP ? $tmp_pos[$#tmp_pos] : $tmp_pos[0];
			
		my $n_new_SV = 0;
		my $tmp_bp;
		foreach my $t (@{$mapping->{$s}}) {
			my $qort = $t->{qstrand};
				my $t_pos = $t->{tstart};
			$t_pos = $t->{tend} if( ($qort eq "+" && $clip == LEFT_CLIP) ||
				($qort eq '-' && $clip == RIGHT_CLIP));
			if($t->{tchr} eq $chr && (abs($t_pos - $real_pos) < $min_dist_diff)) { # the soft clipping is incorrect
				@SV = @SV[0 .. ($#SV - $n_new_SV + 1)] if($n_new_SV > 0);
				last;		
			}
			my $first_bp = {};
			$first_bp->{sc_reads} =  $sclip_count->{$s};
			$first_bp->{$where . "_chr"} = $chr;
			$first_bp->{$where . "_pos"} = $real_pos;
			$first_bp->{all_pos} = \@tmp_pos;
			$first_bp->{$where . "_ort"} = "+";
			$first_bp->{pos}  = $real_pos;
			$first_bp->{cover} = $coverage;
			$first_bp->{chr} = $chr;
			my $new_where = ($where eq "right" ? "left" : "right");
			
			$first_bp->{$new_where . "_chr"} = $t->{tchr};
			$first_bp->{$new_where . "_pos"} = $t_pos;
			$first_bp->{$new_where . "_ort"} = $qort;
			$first_bp->{sc_seq} = $contig_seqs->{$s};
			$first_bp->{reads} = $contig_reads->{$s};

			$tmp_bp = $first_bp unless $tmp_bp;
			my $second_bp = check_sclip(-SAM => $sam, -TARGET => $t, -CHR => $chr, 
				-POS => $real_pos, -SCLIP => $r_sclip, -CLIP => $clip);	
			if(@{$second_bp}) {
				foreach my $tmp_bp (@{$second_bp}) {
					push @SV, StructVar->new(-FIRST_BP => $first_bp, -SECOND_BP => $tmp_bp);
					$n_new_SV++;
				}
			}
		}
		# save some one side good soft-clipping SV
		if( $rescue == 1 &&
			$n_new_SV == 0 && 
			$sclip_count->{$s} >= $min_one_side_reads && 
			(scalar(@{$mapping->{$s}}) == 1 || ($mapping->{$s}[0]{perfect} == 1 && $mapping->{$s}->[1]{perfect} == 0)) &&
			$tmp_bp->{chr} &&
			length($contig_seqs->{$s}) * 0.95 < $mapping->{$s}[0]{matches}	) {

			my %second_bp = %{$tmp_bp};
			$second_bp{sc_reads} = 0;
			$second_bp{sc_seq} = '';
			($second_bp{chr}, $second_bp{pos}) = $second_bp{pos} == $second_bp{left_pos} ? 
				($second_bp{right_chr}, $second_bp{right_pos}) : ($second_bp{left_chr}, $second_bp{left_pos});
			$second_bp{cover} = count_coverage($sam, $second_bp{chr}, $second_bp{pos});
			$second_bp{reads} = [];
			push @SV, StructVar->new(-FIRST_BP => $tmp_bp, -SECOND_BP => \%second_bp);			
		}

	}
	system("rm $fa_name");	system("rm $fa_name*");
	return @SV;
}

sub get_gene_range {
	my ($chr, $pos, $clip) = @_;
	my $r = $clip == LEFT_CLIP ? [$pos, $pos+5] : [$pos - 5, $pos];
	my @genes; 
	my $tmp = $chr;
	$tmp  = "chr" . $tmp if($tmp  !~ m/chr/); 
	my ($f_tree, $r_tree) = ($gm->sub_model($tmp, "+"), $gm->sub_model($tmp, "-"));
	push @genes, $f_tree->intersect($r);
	push @genes, $r_tree->intersect($r);
	my ($s, $e);
	if(scalar @genes == 0) { #no gene here
		return [$pos - $read_len, $pos + $read_len];
	}
	foreach my $g (@genes)	 {
		my $start = $g->val->get_start($pos, $read_len);
		my $end = $g->val->get_end($pos, $read_len);
		$s = $start if(!$s || $s > $start);
		$e = $end if(!$e || $e < $end);
	}
	return [$s, $e];
}

sub check_sclip {
	my %arg = @_;
	my ($sam, $chr, $pos, $target, $r_sclip, $clip) = 
		( $arg{-SAM}, $arg{-CHR}, $arg{-POS}, $arg{-TARGET}, $arg{-SCLIP}, $arg{-CLIP} );

	my @bp;

	# identify the searching region for partner soft cliping
	my $extension;
	if($RNASeq) {
		$extension = 50;
	}
	else {
		$extension  = $read_len - ($target->{qend} - $target->{qstart});
	}
	my ($tpos, $ort) = ($target->{tend}, $target->{qstrand}); 
	$tpos = $target->{tstart} if( ($clip == LEFT_CLIP && $ort eq '-')
		|| ($clip == RIGHT_CLIP && $ort eq '+'));
	$extension = abs($tpos - $pos) - 1 
		if($chr eq $target->{tchr} && abs($tpos - $pos) <= $extension); # don't consider itself
	my $range =$target->{tchr} . ":";
	if($tpos < $extension) {
		$range .= '1';
	}
	else{
		$range .= $tpos - $extension;
	}
	
	$range .= "-"; 	$range .= $tpos + $extension;
	
	# the orginal genome sequence where we find the soft_clipping
	my $orig;
	my $tmp = $chr;
	#$tmp = 'chr' . $tmp if($tmp !~ m/^chr/);
	$orig = $target_genome . ":" . $tmp . ":";
	my $base_pos;
	if($RNASeq) { #let's do a wild guess!
		my $r = get_gene_range($tmp, $pos, $clip);
		return \@bp unless $r;
		$base_pos = $r->[0];
		$orig = join("", ($orig, $r->[0], "-", $r->[1]));
	}
	else {
		$base_pos = $pos - $read_len;
		$orig = join("", ($orig, $pos - $read_len, "-", $pos +  $read_len));
	}

	return \@bp if(!exists( $r_sclip->{$target->{tchr}})); #mapped to chrY or other minor chrs
	# the soft clipping happens in highly repetitive region
	# return \@bp if($coverage > $max_repetitive_cover);

	my $r_pp = search_sclip($r_sclip, $range);
	return  \@bp if(!$r_pp);
	
	my $real_pos;
	my $found = 0;
	my @r_reads;
	my @l_reads;
	foreach my $s (@{$r_pp}) {
		#next if($s->[2] < $min_hit_reads); 
		next if($chr ne $target->{tchr} && (
			($clip == LEFT_CLIP  &&  $ort ne $s->[1] ) ||
			($clip == RIGHT_CLIP && $ort eq $s->[1] )) );

		next if($chr eq $target->{tchr} && ($s->[0] < $tpos - $extension 
			|| $s->[0] > $tpos + $extension));
		my $tort = $s->[1];
		next if($s->[3] > $max_repetitive_cover);
		my $tclip = $tort eq '+' ? RIGHT_CLIP : LEFT_CLIP;
		my @reads = get_sclip_reads(
			-SAM => $sam, 
			-CHR => $target->{tchr}, 
			-START => $s->[0], 
			-END => $s->[0], 
			-CLIP => $tclip,
			-MINCLIPLEN => $min_hit_len,
		   	-VALIDATOR => $validator,
		   	-PAIRED => $paired,
		   	-RMDUP => $rmdup);
		next if(!@reads || scalar(@reads) == 0);
#		push @r_reads, @reads if($tclip == RIGHT_CLIP);
#		push @l_reads, @reads if($tclip == LEFT_CLIP);
#	}

#	foreach my $tclip ( (RIGHT_CLIP, LEFT_CLIP)) {
#		my $s_reads = $tclip == RIGHT_CLIP ? \@r_reads : \@l_reads;
#		next if(scalar @{$s_reads} == 0);
		my %count;
#		my $fa_name = prepare_reads($s_reads, $tclip);
		my $fa_name = prepare_reads(\@reads, $tclip);
		
		my($contig_file, $sclip_count, $contig_reads, $singlet_file) = $assembler->run($fa_name); 
		my (%reads, %quals, %s_lens, %pos);
#		foreach my $r (@{$s_reads}) {
		foreach my $r (@reads) {
			$reads{$r->{name}} = $r->{full_seq};
			$quals{$r->{name}} = $r->{full_qual};
			$s_lens{$r->{name}} = length($r->{seq});
			$pos{$r->{name}} = $r->{pos};
		}

		$contig_file = clip_fa_file($contig_file, $tclip);
		if( -s $singlet_file )  {
			$singlet_file = clip_fa_file($singlet_file, $tclip);
			system("cat $singlet_file >> $contig_file");
			system("rm $singlet_file");
		}
		my $seqs = read_fa_file($contig_file);	
		my $hits = $aligner->run( -TARGET => $orig, -QUERY => $contig_file);
		foreach my $t (keys %{$hits}) {
			my $h = $hits->{$t};
			my $tmp_bp;
			my @all_pos;
			my $all_reads_name;
			if(exists $contig_reads->{$t}) {
				foreach my $n (@{$contig_reads->{$t}}) {
					push @all_pos, $pos{$n};
				}
				$all_reads_name = $contig_reads->{$t};
				@all_pos = uniq @all_pos;
				@all_pos = sort {$a <=> $b} @all_pos;

			}
			else {
				push @all_pos, $pos{$t};
				$all_reads_name = [$t];
			}

			if(is_good_hit($h, $clip, $tclip)) {
				my $hit_ort = ($h->strand('query') == 1 ? "+" : "-");
				my $real_pos = $tclip == RIGHT_CLIP ? $all_pos[0] : $all_pos[$#all_pos];
				if($tclip == RIGHT_CLIP ) {
					$tmp_bp = {
						left_ort	=> '+', 
						left_pos 	=> $real_pos, 
						left_chr 	=> $target->{tchr}, 
						right_ort 	=> $hit_ort,  
						right_pos 	=> ($hit_ort eq "+" ? $h->start("hit") : $h->end("hit")) + $base_pos,
						right_chr 	=> $chr, 
					}
				}
				else {
					$tmp_bp = {
						left_ort 	=> $hit_ort, 
						left_chr 	=> $chr, 
						left_pos 	=> ($hit_ort eq "+" ? $h->end("hit") : $h->start("hit")) + $base_pos,
						right_ort 	=> "+", 
						right_chr 	=> $target->{tchr}, 
						right_pos 	=> $real_pos,
					}
				}
				$tmp_bp->{chr}		= $target->{tchr};
				$tmp_bp->{pos} 		= $real_pos;
				$tmp_bp->{all_pos} 	= \@all_pos;
				$tmp_bp->{sc_seq} 	= $seqs->{$t};
				$tmp_bp->{cover} 	= count_coverage($sam, $target->{tchr}, $real_pos);
				$tmp_bp->{sc_reads} = exists $sclip_count->{$t} ? $sclip_count->{$t} : 1;
				$tmp_bp->{reads} 	= exists $contig_reads->{$t} ? $contig_reads->{$t} : [$t];

				push @bp, $tmp_bp;
			}
		}
		system("rm $fa_name"); system("rm $fa_name*");
	}
	return \@bp;
}

# many more filter can be added here
sub is_good_hit {
	my ($hit, $clip, $tclip) = @_;

	return if($hit->length_aln('query') < $min_hit_len);
	return if($hit->frac_identical * 100 < $min_percent_id );
	return 1;

	# do we want to do the check?
	my $ort = ($hit->start('query') < $hit->end('query') ? "+" : "-");
	my ($dist, $t_dist, $qstart, $qend);
	$qstart = $ort eq '+' ? $hit->start('query') : $hit->end('query') ;
	$qend = $ort eq '+' ? $hit->end('query') : $hit->start('query');
	$t_dist = $qstart;
	$t_dist = $hit->query_length - $qend if($tclip == LEFT_CLIP);

}

# RNASeq support to identify deletions

sub cigar2hit {
	my ($s, @cigar_array) = @_;
	my @pos;
	my $p = $s; 
	foreach my $c (@cigar_array){
		my($op, $len) = @{$c};
		if($op eq 'N') {
			push @pos, [$s, $p];
			$s = $p + $len;
			$p = $s;
			next;
		}
		$p += $len if($op eq 'M' || $op eq 'D');
	}
	push @pos, [$s, $p];
	return \@pos;
}

sub	find_del_SVs  {
	my ($sam, $chr, $start, $end) = @_;
	my $max_sclip_len = 0;
	my $range = $chr;
	if($start && $end) {
		$range = $chr . ":" . $start . "-" . $end;
	}
	my %SV;
	my $tmp = $chr;
	$tmp = "chr" . $chr if($chr !~ m/chr/);
 
	my($f_tree, $r_tree) = ($gm->sub_model($tmp, "+"), $gm->sub_model($tmp, "-"));

	$sam->fetch($range, sub {
		my $a = shift;
		my $cigar_str = $a->cigar_str;
#		print $a->qname, "\t", $cigar_str, "\n";
		return unless $cigar_str =~ m/N/; # the read cross two genes must have N
#		return unless $a->score >= 97;
		my $sv = identify_del_SV($f_tree, $a);
		if($sv){
			my $tmp1 = $sv->[0] . "+" . $sv->[1];
			if(!exists($SV{$tmp1})) {
				$SV{$tmp1} = {};
				$SV{$tmp1}->{num} = 0;
				$SV{$tmp1}->{left} = 0;
				$SV{$tmp1}->{right} = 0;
			}
			$SV{$tmp1}->{num}++;
			$SV{$tmp1}->{left} = $sv->[5] if($SV{$tmp1}->{left} < $sv->[5]);
			$SV{$tmp1}->{right} = $sv->[6] if($SV{$tmp1}->{right} < $sv->[6]);
			$SV{$tmp1}->{gene_ort}  = "+";
			$SV{$tmp1}->{left_gene} = $sv->[3];
			$SV{$tmp1}->{right_gene} = $sv->[4];
			push @{$SV{$tmp1}->{reads}}, $a->qname;

		}
		$sv = identify_del_SV($r_tree, $a);
		if($sv){
			my $tmp1 = $sv->[0] . "-" . $sv->[1];
			if(!exists($SV{$tmp1})) {
				$SV{$tmp1} = {};
				$SV{$tmp1}->{num} = 0;
				$SV{$tmp1}->{left} = 0;
				$SV{$tmp1}->{right} = 0;
				$SV{$tmp1}->{reads} = [];
			}
			$SV{$tmp1}->{num}++;
			$SV{$tmp1}->{left} = $sv->[5] if($SV{$tmp1}->{left} < $sv->[5]);
			$SV{$tmp1}->{right} = $sv->[6] if($SV{$tmp1}->{right} < $sv->[6]);
			$SV{$tmp1}->{gene_ort}  = "-";
			$SV{$tmp1}->{left_gene} = $sv->[3];
			$SV{$tmp1}->{right_gene} = $sv->[4];
			push @{$SV{$tmp1}->{reads}}, $a->qname;
		}
	}
	);

	my @rtn;
	foreach my $sv (keys %SV) {
		my( $pos1, $pos2 ) = $sv =~ m/(\d+)[+|-](\d+)/;
		print $pos1, "\t", $pos2, "\n";
		my( $cover1, $cover2) = (count_coverage($sam, $chr, $pos1 - 1), count_coverage($sam, $chr, $pos2));
		next if($SV{$sv}->{num} < $min_sclip_reads);
		my $cover = $cover1 < $cover2 ? $cover1 : $cover2;
		next if($SV{$sv}->{num} * 100 < $cover *  $min_pct_sclip_reads && $SV{$sv} < 50);
		push @rtn, StructVar->new (
			-FIRST_BP => {	left_ort => '+',
							left_pos => $pos1 - 1,
							left_chr => $chr,
							right_ort => '+',
							right_pos => $pos2,
							right_chr => $chr,
							chr => $chr,
							cover => count_coverage($sam, $chr, $pos1 - 1), 
							pos => $pos1 - 1,
							sc_reads => $SV{$sv}->{num},
							left_len => $SV{$sv}->{left},
							right_len => $SV{$sv}->{right},
							gene_ort => $SV{$sv}->{gene_ort},
							left_gene => $SV{$sv}->{left_gene},
							right_gene => $SV{$sv}->{right_gene},
							reads => $SV{$sv}->{reads},
						},
			-SECOND_BP => {	left_ort => '+',
							left_pos => $pos1 - 1,
							left_chr => $chr,
							right_ort => '+',
							right_pos => $pos2,
							right_chr => $chr,
							chr => $chr,
							pos => $pos2,
							cover => count_coverage($sam, $chr, $pos2),
							sc_reads => $SV{$sv}->{num},
							reads => [],
						}
		);
	}
	return @rtn;
}

sub identify_del_SV {
	my ($tree, $a) = @_;
	my @cigar_array = @{$a->cigar_array};
	return unless($a->start && $a->end);
	return if(scalar $tree->intersect([$a->start, $a->end]) <= 1);

	my @hits = @{ cigar2hit($a->start, @cigar_array) };
	my @genes;
	my $last_gene;
	my @change;
	my ($left, $right) = (0, 0);
	for( my $i = 0; $i < scalar @hits; $i++)  {
		my $h = $hits[$i];
		my @g = $tree->intersect($h);
		my $g_size = scalar @g;
		if($g_size == 0) {
#			print STDERR $a->qname, " is mapped outside of gene\n";
			next;
		}
		if( scalar @g == 1) {
			push @genes, $g[0] unless ($last_gene && $last_gene->val->name eq $g[0]->val->name);
			if($last_gene && $last_gene->val->name ne $g[0]->val->name) {
				$right += ($h->[1] - $h->[0] );
				push @change, $i;	
			}
			$left += ($h->[1] - $h->[0] ) if($right == 0); 
			$last_gene = $g[0];
			next;
		}
#		print STDERR $a->qname, " connect two genes into one?\n";
		return;
	}
	return if(scalar @genes <= 1);
	if(scalar @genes > 2) {
#		print STDERR $a->qname, "crossed more than 2 genes!\n";
		return;
	}
	# now we down to 2 genes
	my ($left_gene, $right_gene) = @genes;	
	return unless $left_gene->val->overlap([$hits[$change[0]-1]->[0], $hits[$change[0]-1]->[1]]);
	return unless $right_gene->val->overlap([$hits[$change[0]]->[0], $hits[$change[0]]->[1]]);
	my $p = $left_gene->val->end;
	print join("\t", ( $a->score, 
		$a->cigar_str, $left, $right, $hits[$change[0]-1]->[1], $hits[$change[0]]->[0], 
		$left_gene->val->name, $right_gene->val->name)), "\n";	
	return [$hits[$change[0]-1]->[1], $hits[$change[0]]->[0], $a, $left_gene->val, $right_gene->val, $left, $right];	
}

=head1 NAME

CREST.pl - a Structure Variation detection tools using softclipping for
whole genome sequencing.


=head1 VERSION

This documentation refers to CREST.pl version 0.0.1.


=head1 USAGE
	
	This program depends on several things that need to be installed and/or
	specified.  The program uses BioPerl and Bio::DB::Sam module to parse 
	the files and bam files.  Also it uses Blat software suits to do genome
	mapping and alignment.  To make the program efficient, it also requires
	a blat server setup.  And the program uses CAP3 assembler.

	Identify somatic SVs from input:
	    CREST.pl -f somatic.cover -d tumor.bam -g germline.bam 
	Identify SVs from a single bam compared wtih reference genome:
		CREST.pl -f sclip.cover -d sample.bam 
	Identify SVs from a single bam on chr1 only
		CREST.pl -f sclip.cover -d sample.bam -r chr1

=head1 REQUIRED ARGUMENTS

	To run the program, several parameter must specified.
	-d, --input_d:		 The input (diagnositic) bam file
	-f, --sclipfile:	 The soft clipping information generated by extractSClip.pl
	--ref_genome:		 The reference genome file in fa format
	-t, --target_genome: The 2bit genome file used by blat server and blat
	--blatserver:		 The blat server name or ip address
	--blatport:			 The blat server port

=head1 OPTIONS

	The options that can be used for the program.
	-g, --input_g		 The germline (paired) bam file
	-p, --prefix		 The output prefix, default is the input bam file name
	-o, --out_dir 		 Output directory, default is the working directory
	-l, --read_len		 The read length of the sequencing data, defaut 100
	--sensitive          The program will generate more SVs with higher false positive rate.

	--(no)scratch 		 Use the scratch space, default is off.
	--(no)paired		 Use paired reads or not, defafult is on, so change to --nopaired for unpaired reads.
	--cap3				 CAP3 executable position, default cap3
	--cap3opt			 CAP3 options, single quoted, Default ' > /dev/null'
	--blatclient 		 gfClient excutable position, default gfClient
	--blatclientopt		 gfClient options, single quoted, default '-out=psl -nohead'
	--blat				 blat executable potion, default balt
	--blatopt			 blat options, single quoted, default '-tileSize=7 -stepSize=1 -out=psl -minScore=15'
	
	-r, --range			 The range where SV will be detected, using chr1:100-200 format.
	--max_score_diff	 The maximum score difference when stopping select hit, default 10.
	--min_sclip_reads	 Minimum number of soft clipping read to triger the procedure, default 3.
	--max_rep_cover		 The min number of coverage to be called as repetitive and don't triger
						 the procedure, default 500.
	--min_sclip_len		 The min length of soft clipping part at a position to triger the detection,
						 default 20.
	--min_hit_len		 Min length of a hit for genome mapping
	--min_dist_diff		 Min distance between the mapped position and the soft clipping position, default 20.
	--(no)rmdup			 Remove PCR dumplicate.  Default remove.
	
	--min_percent_id	 Min percentage of identity of soft clipping read mapping, default 90
	--min_percent_hq	 Min percentage of high quality base in soft clipping reads, default 80
	--lowqual_cutoff	 Low quality cutoff value, default 20.

	-h, --help, -?		 Help information
	--man				 Man page.
	--usage				 Usage information.
	--version			 Software version.

	--(no)rm_tandem_repeat	Remove tandem repeat caused SV events, default is ON. When it's on ptrfinder program
							need to be on the path.
	--tr_max_indel_size		Maximum tandem repeat mediated INDEL events, default 100
	--tr_min_size			Minimum tandem reapet size, default 2
	--tr_max_size			Maximum tandem repeat size, default 8
	--tr_min_num			Minimum tandem repeat number, defaut 4
	--min_percent_cons_of_read Minimum percent of consensus length of read length, default 0.75
	--max_bp_dist			Maximum distance between two idenfitifed break points, default 15
	--germline_seq_width	Half window width of genomic sequence around break point for germline SV filtering,
							default 100
	--germline_search_width	Half window width for seaching soft-clipped reads around breakpoint for germline SV
							filtering, default 50.

	--hetero_factor			The factor about the SV's heterogenirity and heterozygosity, default 0.4.
	--triger_p_value		The p-value that will triger the SV detection when number of soft-clipped reads is small,
							defaut 0.05.

	--(no)rescue			Use rescue mode, when it's on, the a SV with only 1 side with enough soft-clipped reads
							is considered as a valid one instead of rejecting it.  Default on.
	--min_one_side_reads	When rescure mode is on, the minimum number of soft-clipped reads on one side, default 5.

	--RNASeq				RNAseq mode, default off
	--genemodel				A gene model file, currently only refFlat format (BED) is supported. Requried for RNASeq.
	--cluster_size			The soft-clipped reads within cluster_size will be considered together, default is 3, RNAseq mode
							only.


=head1 DESCRIPTION

This is a program designed to identify Structure Variations (SVs) using soft
clipping reads.  With the improvement of next-gen sequencing techniques,
both the coverage and length of pair-end reads have been increased steady.
Many SV detection software availalbe uses pair-end information due to the 
limitted read length.  Now 100bp is pretty common and many reads will cross
the break points. Some mapping programs ( bwa, etc ) have the ability to
identify so called soft-clipping reads.  A soft-clipping read is a read that
different part can be mapped to different genomic posiiton, but the read can
be uniquely positioned using the mate-pair position and the insert length.
So this program use those reads to do an assembly-mapping-assembly-alignment
procedure to identify potential structure variiations.


=head1 DIAGNOSTICS

A list of every error and warning message that the application can generate
(even the ones that will "never happen"), with a full explanation of each
problem, one or more likely causes, and any suggested remedies. If the
application generates exit status codes (e.g., under Unix), then list the exit
status associated with each error.

=head1 CONFIGURATION AND ENVIRONMENT

The program is designed to use under a cluster or high performance computing
environment since it's dealing with over 100G input data.  The program can be 
used as highly parallellized.  And the program is developped under linux/unix.


=head1 DEPENDENCIES

The program depend on several packages:
1. Bioperl perl module.
2. Bio::DB::Sam, version 1.5 or later, it requires samtools lib installed.
3. blat suits, include blat, gfClient, gfServer.
4. CAP3 assembly program.


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
