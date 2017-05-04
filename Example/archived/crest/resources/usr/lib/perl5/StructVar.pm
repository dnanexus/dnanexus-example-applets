package StructVar;

# $Id: StructVar.pm, v 1.00 2010/09/24 jianmin.wang Exp $

use strict;
use Bio::DB::Sam;
use Carp;
use Data::Dumper;
use English;
use Bio::SearchIO;
use GeneModel;
use Cwd;
use SVExtTools;
use List::MoreUtils qw/ uniq /;
use SVUtil qw( prepare_reads read_fa_file get_sclip_reads rev_comp);
use SCValidator qw(LEFT_CLIP RIGHT_CLIP);

use constant UNK => 0;
use constant CTX => 1;
use constant INV => 2;
use constant INS => 3;
use constant DEL => 4;
use constant ITX => 5;

# type of SVs and thier string
my %type2str = (
	0 => 'UNKNOWN',  # unknow type, not single event 
	1 => 'CTX',		 # cross chromosome translocation
	2 => 'INV',		 # inversion, requires both breakpoints available
	3 => 'INS',	     # (tandem) insertion
	4 => 'DEL',		 # deletion
	5 => 'ITX',	     # within chromosome translocation
);
# class static variables
my $sam_d;
my $sam_g;
my $gm;
my $assembler;
my $aligner;
my $mapper;
my @filters = (\&_germline_sclip_filter, \&_type_distance_filter, \&_germline_indel_filter);
my $RNASeq;
my $read_len;
my $genome;
my $tr_max_indel_size;   #tandem repeat mediated indel size
my $tr_min_size;
my $tr_max_size; 	  #tandem repeat max_size of single repeat
my $tr_min_num;		  #tandem repeat minimum number of occurence

my $max_bp_dist = 15;
my $germline_seq_width = 100;
my $germline_search_width = 50;
my $min_percent_cons_of_read = 0.75;

my %default_filters = (
	'germline_sclip'	=> \&_germline_sclip_filter,
	'type_distance'		=> \&_type_distance_filter,
	'germline_indel'	=> \&_germline_indel_filter,
#	'mapping_artifact'	=> \&_mapping_artifact_filter,
	'mapping_quality'	=> \&_mapping_quality_filter,
	'tandem_repeat'		=> \&_tandem_repeat_filter,
);

sub new {
	my $class = shift;
	my %param = @_;
	my $self = {};
	$self->{first_bp} = {};
	$self->{second_bp} = {};
	$self->{first_bp} = $param{-FIRST_BP} if($param{-FIRST_BP}); 
	$self->{second_bp} = $param{-SECOND_BP} if($param{-SECOND_BP});
	$self->{consensus} = $param{-CONSENSUS} ? $param{-CONSENSUS} : "";
	$self->{type} = "";
	bless $self, ref($class) || $class;
	return $self;
}

sub add_filter {
	my $self = shift;
	my $name = shift;
	my $f = shift;
	croak "you must specify a filter" if(!$f);
	croak "the filter must be a subroutine" if(ref($f) ne 'CODE');
	$default_filters{$name} = $f;
}

sub add_RNASeq_filter {
	$default_filters{'RNASeq_strand'} = \&_RNASeq_strand_filter;
	$default_filters{'RNASeq_INS'}	= \&_RNASeq_INS_filter;
}

sub remove_filter {
	my $self = shift;
	my $name = shift;
	if($default_filters{$name}){
		delete $default_filters{$name};
	}
	else {
		print STDERR "Filter $name is not in filter list";
	}
} 

sub tr_max_indel_size {   #tandem repeat mediated indel size
	my $self = shift;
	my $value = shift;
	$tr_max_indel_size = $value if($value);
	return $tr_max_indel_size;
}

sub min_percent_cons_of_read {
	my $self= shift;
	my $value = shift;
	$min_percent_cons_of_read = $value if($value);
	return $min_percent_cons_of_read;
}

sub tr_min_size {
	my $self = shift;
	my $value = shift;
	$tr_min_size = $value if($value);
	return $tr_min_size;
}
sub germline_seq_width {
	my $self = shift;
	my $value = shift;
	$germline_seq_width = $value if($value);
	return $germline_seq_width;
}

sub germline_search_width {
	my $self = shift;
	my $value = shift;
	$germline_search_width = $value if($value);
	return $germline_search_width;
}
	
sub max_bp_dist {
	my $self = shift;
	my $value = shift;
	$max_bp_dist = $value if($value);
	return $max_bp_dist;
}

sub tr_max_size { 	  #tandem repeat max_size of single repeat
	my $self = shift;
	my $value = shift;
	$tr_max_size = $value if($value);
	return $tr_max_size;
}

sub tr_min_num {	  #tandem repeat minimum number of occurence
	my $self = shift;
	my $value = shift;
	$tr_min_num = $value if($value);
	return $tr_min_num;
}
sub genome {
	my $self = shift;
	my $value = shift;
	$genome = $value if($value);
	return $genome;
}
sub sam_d {
	my $self = shift;
	my $value = shift;
	$sam_d = $value if($value);
	return $sam_d;
}

sub sam_g {
	my $self = shift;
	my $value = shift;
	$sam_g = $value if($value);
	return $sam_g;
}

sub assembler {
	my $self = shift;
	my $value = shift;
	$assembler = $value if($value);
	return $assembler;
}

sub aligner {
	my $self = shift;
	my $value = shift;
	$aligner = $value if($value);
	return $aligner;
}

sub RNASeq {
	my $self = shift;
	my $value = shift;
	$RNASeq = 1 if($value);
	return $RNASeq;
}

sub mapper {
	my $self = shift;
	my $value = shift;
	$mapper = $value if($value);
	return $mapper;
}

sub gene_model {
	my $self = shift;
	my $value = shift;
	$gm = $value if($value);
	return $gm;
}

sub read_len {
	my $self = shift;
	my $value = shift;
	$read_len = $value if($value);
	return $read_len;

}
sub first_bp {
	my $self = shift;
	my %param = @_;
	if(%param) {
		$self->{first_bp}{left_chr}  = $param{-LEFT_CHR};
		$self->{first_bp}{left_pos}  = $param{-LEFT_POS};
		$self->{first_bp}{left_ort}  = $param{-LEFT_ORT};
		$self->{first_bp}{right_chr} = $param{-RIGHT_CHR};
		$self->{first_bp}{right_pos} = $param{-RIGHT_POS};
		$self->{first_bp}{right_ort} = $param{-RIGHT_ORT};
		$self->{first_bp}{sc_reads}  = $param{-SC_READS};
		$self->{first_bp}{pos}       = $param{-POS};
		$self->{first_bp}{cover}     = $param{-COVER};
		$self->{first_bp}{sc_seq}	 = $param{-SC_SEQ};
		$self->{first_bp}{chr}		 = $param{-CHR};
	}
	return $self->{first_bp};
}

sub second_bp {
	my $self = shift;
	my %param = @_;
	if(%param) {
		$self->{second_bp}{left_chr}  = $param{-LEFT_CHR};
		$self->{second_bp}{left_pos}  = $param{-LEFT_POS};
		$self->{second_bp}{left_ort}  = $param{-LEFT_ORT};
		$self->{second_bp}{right_chr} = $param{-RIGHT_CHR};
		$self->{second_bp}{right_pos} = $param{-RIGHT_POS};
		$self->{second_bp}{right_ort} = $param{-RIGHT_ORT};
		$self->{second_bp}{sc_reads}  = $param{-SC_READS};
		$self->{second_bp}{pos}       = $param{-POS};
		$self->{second_bp}{cover}     = $param{-COVER};
		$self->{second_bp}{sc_seq}	 = $param{-SC_SEQ};
		$self->{second_bp}{chr}		 = $param{-CHR};
	}
	return $self->{second_bp};
}

sub type {
	my $self = shift;
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	if($bp1->{left_chr} ne $bp1->{right_chr}) {
		return CTX if($bp1->{left_ort} eq $bp2->{left_ort} &&
					  $bp1->{right_ort} eq $bp2->{right_ort});

		return UNK;
	}
	if( $bp1->{left_ort} eq $bp1->{right_ort} &&
		$bp2->{left_ort} eq $bp2->{right_ort} ) { # Insertion or deletion
		return DEL if( $bp1->{left_pos} <= $bp1->{right_pos} &&
					   $bp2->{left_pos} <= $bp2->{right_pos});
		return INS if( $bp1->{left_pos} >= $bp1->{right_pos} &&
					   $bp2->{left_pos} >= $bp2->{right_pos});
		return UNK;
	}
	if( $bp1->{left_ort} ne $bp1->{right_ort} &&
		$bp2->{left_ort} ne $bp2->{right_ort} ) {
		return INV if( ($bp1->{left_ort} eq '+' && $bp2->{left_ort} eq '-')
					|| ($bp2->{left_ort} eq '+' && $bp1->{left_ort} eq '-'));
		return ITX if( ($bp1->{left_ort} eq '+' && $bp2->{left_ort} eq '+')  
					|| ($bp1->{left_ort} eq '-' && $bp2->{left_ort} eq '-') );
		return UNK;
	}
	return UNK;
}

sub to_string {
	my $self = shift;
	my $type = $self->type;
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	my ($left_len, $right_len) = (0, 0);
	if(exists $bp1->{left_len}) {
		($left_len, $right_len ) = ($bp1->{left_len}, $bp1->{right_len});
	}
	else {
		$left_len = length($bp1->{sc_seq}) if($bp1->{sc_seq}); 
		$right_len = length($bp2->{sc_seq}) if($bp2->{sc_seq}); 
	}
	if($type == INV && $bp1->{pos} > $bp2->{pos}) {
		($bp1, $bp2) = ($bp2, $bp1);
		($left_len, $right_len) = ($right_len, $left_len);
	}
	my ($cover1, $cover2) = ($bp1->{cover}, $bp2->{cover});

	my ($chrA, $chrB, $ortA, $ortB) = ($bp1->{left_chr}, $bp1->{right_chr},
		$bp1->{left_ort}, $bp1->{right_ort});
	my ($posA, $posB, $countA, $countB) = ($bp1->{pos}, $bp2->{pos}, $bp1->{sc_reads}, $bp2->{sc_reads});
	if($bp1->{chr} eq $bp1->{right_chr} && $bp1->{pos} == $bp1->{right_pos}) {
		$posA = $bp2->{pos};
		$posB = $bp1->{pos};
		$countA = $bp2->{sc_reads};
		$countB = $bp1->{sc_reads};
		($cover1, $cover2) = ($cover2, $cover1);
		($left_len, $right_len) = ($right_len, $left_len);
	}

	my $rtn = join("\t", ($chrA, $posA, $ortA, $countA, $chrB, $posB, $ortB, $countB,
		$type2str{$type}, $cover1, $cover2, $left_len, $right_len, $self->get_statistics(100),
		$self->{c_start}, $self->{start_chr}, $self->{t_start}, 
		$self->{c_end}, $self->{end_chr}, $self->{t_end}, $self->{consensus}));
	if($self->{type} == INV) {
		$rtn = join("\t", ($rtn, $self->{c_start2}, $self->{start_chr2}, $self->{t_start2},
			$self->{c_end2}, $self->{end_chr2}, $self->{t_end}, $self->{consensus2}));
	}
	return $rtn;
}

sub set_consensus {
	my $self = shift;
	my ($validator, $paired ) = @_;
	if(!$assembler) {
		print STDERR "No assembler set, no consensus generated\n";
		return;
	}
	my %all_reads;
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	my $N = 0;

	my @s_reads1 = $self->_get_reads($sam_d, $bp1);
	my @s_reads2 = $self->_get_reads($sam_d, $bp2);
	
	if($self->{type} == INV) {
		$self->{consensus} = _generate_consensus(\@s_reads1);
		$self->{consensus2} = _generate_consensus(\@s_reads2);
	}
	else {
		push @s_reads1, @s_reads2;
		$self->{consensus} = _generate_consensus(\@s_reads1);
	}
}

sub _generate_consensus {
	my $s_reads = shift;
	my $N = scalar(@{$s_reads});
	return $s_reads->[0]{seq} if($N == 1);
	my $fa_name = prepare_reads($s_reads);
	my ($contig_file, $sclip_count) = $assembler->run($fa_name);
	if(!$contig_file || -s $contig_file == 0) { 
		system ("rm $fa_name"); system("rm $fa_name*");
		return;
	}
	my $n = 0;
	my $contig_name;
	foreach my $c (keys %{$sclip_count}) {
		if($sclip_count->{$c} > $n) {
			$n = $sclip_count->{$c};
			$contig_name = $c;
		}
	}
	return if($N * 8 > $n * 10);
	my $contig_seqs = read_fa_file($contig_file);
	return $contig_seqs->{$contig_name};
}

sub _get_reads{
	my $self = shift;
	my ($sam, $bp) = @_;

	my ($chr, $pos) = ($bp->{chr}, $bp->{pos});
	my ($start, $end) = ($pos, $pos);
	if($bp->{all_pos}) {
		my @tmp = @{$bp->{all_pos}};
		@tmp = sort {$a <=> $b} @tmp;
		($start, $end) = ($tmp[0], $tmp[$#tmp]);
	}
	my %reads;
	return if(scalar @{$bp->{reads}} == 0);
	foreach my $r (@{$bp->{reads}}){
		$reads{$r} = 1;
	}

	#my ($sam, $chr, $start, $end, $reads) = 
	#	($args{-SAM}, $args{-CHR}, $args{-START}, $args{-END}, $args{-READS}) ;

	my @rtn;
	my $range = $chr;
	$range = $chr . ":" . $start . "-" . $end;
	my($s, $e) = ($start, $end);
	$sam->fetch($range, sub {
		my $a = shift;
		return unless (exists $reads{$a->qname});
		return unless ($a->start && $a->end);
		return unless (($a->start >= $start && $a->start <= $end) 
			|| ($a->end >= $start && $a->end <= $end));
		$s = $a->start if($s > $a->start);
		$e = $a->end if($e < $a->end);
		my $qscore = $a->qscore;
		delete $reads{$a->qname};
		push @rtn, {
			name => $a->qname,
			seq => $a->query->dna,
			qual => $qscore,
		};
	} );
#	$bp->{range} = [$s, $e];
	if($self->{type} == INV) {
		if($start == $bp->{left_pos} || $end == $bp->{left_pos}) {
			$bp->{left_range} = [$s, $e];
		}
		else {
			$bp->{right_range} = [$s, $e];
		}
		return @rtn;
	}
	if($chr eq $self->{first_bp}{left_chr} && ($start == $self->{first_bp}{left_pos} || $end == $self->{first_bp}{left_pos})){
		$self->{first_bp}{left_range} = [$s, $e];
	}
	if($chr eq $self->{first_bp}{right_chr} && ($start == $self->{first_bp}{right_pos} || $end == $self->{first_bp}{right_pos})){
		$self->{first_bp}{right_range} = [$s, $e];
	}

	return @rtn;	
}

my $start_dir;
BEGIN {
	$start_dir = getcwd;
}

sub _cat_genes {
	my ($ort, @genes) = @_;
	my $rtn;
	my @names;
	foreach my $g (@genes) {
		push @names, $g->val->name;
	}
	@names = uniq @names;
	$rtn = join(";", @names);
	$rtn = $rtn . "($ort)" if($rtn);
	return $rtn;
}

sub get_genes {
	my $self = shift;
	return ['NA', 'NA', 'NA'] unless $gm;

	my ($gene1, $gene2);
	my $dist;
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	my ($ort1, $ort2) = ($bp1->{left_ort}, $bp1->{right_ort});
	my ($chr1, $chr2) = ($bp1->{left_chr}, $bp1->{right_chr});
	my ($pos1, $pos2) = ($bp1->{pos}, $bp2->{pos});
	($pos1, $pos2) = ($pos2, $pos1) if($bp1->{right_pos} == $pos1);
	my $range1 = $ort1 eq "+" ? [$pos1 - 5, $pos1] : [$pos1, $pos1 + 5];
	my $range2 = $ort2 eq "+" ? [$pos2, $pos2 + 5] : [$pos2 - 5, $pos2];

	my ($tmp1, $tmp2) = ($chr1 =~ m/chr/ ? $chr1 : "chr" . $chr1, 
		$chr2 =~ m/chr/ ? $chr2 : "chr" . $chr2);
	my ($r_tree1, $f_tree1) = ($gm->sub_model($tmp1, "-"), 	$gm->sub_model($tmp1, "+"));
	my ($r_tree2, $f_tree2) = ($gm->sub_model($tmp2, "-"), 	$gm->sub_model($tmp2, "+"));

	my ($g_ort1, $g_ort2);
	my @genes;
	@genes = $r_tree1->intersect($range1) if($r_tree1);
	$gene1 = _cat_genes("-", @genes) if(scalar @genes > 0 );
	undef @genes;
	@genes = $f_tree1->intersect($range1) if($f_tree1);
	if(scalar @genes > 0) {
		my $tmp = _cat_genes("+", @genes);
		if($gene1) {
			$gene1 .= ";" .  $tmp;
		}
		else {
			$gene1 = $tmp;
		}
	}
	undef @genes;
	@genes = $r_tree2->intersect($range2) if($r_tree2);
	$gene2 = _cat_genes("-", @genes) if(scalar @genes > 0 );
	undef @genes;
	@genes = $f_tree2->intersect($range2) if($f_tree2);;
	if(scalar @genes > 0) {
		my $tmp = _cat_genes("+", @genes);
		if($gene2) {
			$gene2 .= ";" .  $tmp;
		}
		else {
			$gene2 = $tmp;
		}
	}

	$gene1 = 'NA' unless $gene1;
	$gene2 = 'NA' unless $gene2;
	my $type = $self->{type};
	if( $type == INS  or $type == DEL ) {
		($pos1, $pos2) = ($pos2, $pos1)  if($pos1 > $pos2);
		@genes =  $r_tree1->intersect([$pos1 - 5, $pos2 + 5]);
		if(scalar @genes >= 2) {
			$dist = scalar(@genes) - 2;
			$dist .= "(-)";
		}
		@genes =  $f_tree1->intersect([$pos1 - 5, $pos2 + 5]);
		if(scalar @genes >= 2) {
			if($dist) {
				$dist .= ":" . (scalar (@genes) - 2);
			}
			else {
				$dist = scalar (@genes) - 2;
			}
			$dist .= "(+)";
		}
	}
	$dist = 'NA' unless $dist;
	return [$gene1, $gene2, $dist];
}

sub to_full_string {
	my $self = shift;
	my $type = $self->{type};
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	return join("\t", ( $bp1->{left_chr}, $bp1->{left_pos}, $bp1->{left_ort},
		$bp1->{right_chr}, $bp1->{right_pos}, $bp1->{right_ort}, $bp1->{sc_reads},
		$bp2->{left_chr}, $bp2->{left_pos}, $bp2->{left_ort},
		$bp2->{right_chr}, $bp2->{right_pos}, $bp2->{right_ort}, $bp2->{sc_reads},
		$type2str{$type}));
}

sub filter {
	my $self = shift;
	print STDERR "SV filter starting....\n";
	$self->{type} = $self->type;
	$self->set_consensus;

	if(!$self->{consensus}  || length $self->{consensus} < $read_len * $min_percent_cons_of_read) {
		$self->{error} = "No Consensus or consensus too short";
		print STDERR "FAILED\n";
		return;
	}
	if($self->{type} == INV && (!$self->{consensus2} || length $self->{consensus2} < $read_len * $min_percent_cons_of_read)) {
		$self->{error} = "No Consensus or consensus too short";
		print STDERR "FAILED\n";
		return;
	}

	foreach my $f (values(%default_filters) ){
		if(! $f->($self)){
			print STDERR "FAILED\n";
			return;
		}
	}
	print STDERR "PASSED\n";
	return 1;
}

sub _mapping_quality_filter {
	my $self = shift;
	if($self->{type} == INV) {
		my $tmp1 = $self->_mapping_quality($self->{first_bp}, $self->{consensus}, 1);
		my $tmp2 = $self->_mapping_quality($self->{second_bp}, $self->{consensus2}, 2);
		return $tmp1 && $tmp2;
	}
	return $self->_mapping_quality($self->{first_bp}, $self->{consensus}, 1);
}

sub _mapping_quality {
	my $self = shift;
	my $bp = shift;
	my $con = shift;
	my $which = shift;
	$which = "" unless $which == 2;
	my $l = length($con);
	print STDERR "Mapping quality filter ... ";

	open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
	print $QUERY ">query\n", $con, "\n";
	close $QUERY;

	my ($range1, $range2) = $self->_find_bp_ref_range($bp, $con);
	my ($s1, $e1, $t_s1, $t_e1) = $self->_map_consensus($range1);
	return unless $e1;
	my ($chr, $s, $e) = split /[:|-]/, $range2;
	my $offset = 0;
	if($s < $t_e1 && $e > $t_s1) { #overlap
		my $tmp = substr $con, 0, $s;
		if($s1 < $l - $e1) {
			$offset = $e1;
			$tmp = substr $con, $e1 + 1;
		}
		open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
		print $QUERY ">query\n", $tmp, "\n";
		close $QUERY;
	}
	my ($s2, $e2, $t_s2, $t_e2) = $self->_map_consensus($range2);
	return unless $e2;
	$s2 += $offset;
	$e2 += $offset;
	if($s1 < $s2) {
		if($e1 < $s2) {
			$self->{"insert" . $which} = substr($con, $e1, $s2 - $e1);
		}
		$self->{"c_start" . $which} = $s1;
		$self->{"t_start" . $which} = $t_s1;
		$self->{"start_chr" . $which} = $bp->{left_chr};
		$self->{"c_end" . $which} = $e2;
		$self->{"t_end" . $which} = $t_e2;
		$self->{"end_chr". $which} = $bp->{right_chr};
		$s1 = 1; 
		($e1, $s2)  = ($s2, $e1);
		$e2 = $l;
	}
	else {
		if($e2 < $s1) {
			$self->{"insert" . $which} = substr($con, $e2, $s1 - $e2);
		}
		$self->{"c_start" . $which} = $s2;
		$self->{"t_start" . $which} = $t_s2;
		$self->{"start_chr" . $which} = $bp->{right_chr};
		$self->{"c_end" . $which} = $e1;
		$self->{"t_end" . $which} = $t_e1;
		$self->{"end_chr" . $which } = $bp->{left_chr};
		$s2 = 1;
		($e2, $s1) = ($s1, $e2);
		$e1 = $l;	
	}
	
	my $n = 1;
	foreach my $r( ($range1, $range2) ) {
		my $r_a = $r eq $range1 ? $range2 : $range1;
		my($chr, $s, $e) = split /[:|-]/, $r;
		my($chr_a, $s_a, $e_a) = split /[:|-]/, $r_a;
		my($s_c, $e_c) = $n == 1 ? ($s1, $e1) : ($s2, $e2);
		$n++;
		open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
		print $QUERY ">query\n", substr($self->{consensus}, $s_c, $e_c - $s_c + 1), "\n";
		close $QUERY;
		my $hit = $aligner->run(-TARGET => $genome . ":" . $r, -QUERY => "query.fa");
		my $h = $hit->{query};
		return unless $h;	
		my $hsp = $h->next_hsp;
		my @matches = $hsp->matches;
		$hit = $mapper->run(-QUERY => "query.fa");	
		foreach $h (@{$hit->{query}}) {
			next if($h->{tchr} eq $chr && $h->{tstart} >= $s - $l && $h->{tend} <= $e + $l); # the same one
			return if($h->{matches} >= $matches[0]); #find a better or as good as match
			return if($h->{tchr} eq $chr_a && $h->{matches} - $matches[0] >= -5 && $self->{type} == CTX);
		}
	}
	print STDERR "PASSED\n";
	return 1;
}

sub _map_consensus {
	my $self = shift;
	my $range = shift;
	my ($s_c, $e_c, $s_t, $e_t);

	my $hits = $aligner->run(-TARGET => $genome . ":" . $range, -QUERY => "query.fa");
	my ($chr, $s, $e) = split /[:|-]/, $range;
	my $h = $hits->{query};
	return unless $h; # can't find a good enough match
	my $hsp = $h->next_hsp;
	my @matches = $hsp->matches;
	($s_c, $e_c, $s_t, $e_t) = ($hsp->start('query'), $hsp->end('query'),
		$hsp->start('hit'), $hsp->end('hit'));
	return if( ($e_c - $s_c + 1) * 0.97 > $matches[0]); # the alignment is not good
	($s_t, $e_t) = $hsp->strand == 1 ?  ($s_t + $s, $e_t + $s) : ($e_t + $s, $s_t + $s);
	return ($s_c, $e_c, $s_t, $e_t);
}

sub _find_bp_ref_range {
	my $self = shift;
	my $bp = shift;
	my $con = shift;

	my $l = int(length($con)/2);
	my $ext = $l * 2;
	$ext = 100000 if($RNASeq);
	if($self->{type} == DEL && abs($bp->{left_pos} - $bp->{right_pos}) < $l ){
		$l = abs($bp->{left_pos} - $bp->{right_pos}) - 1;
	}

	my($range1, $range2);
	my ($chr, $p) = ($bp->{left_chr}, $bp->{left_pos});
	my 	$ort = $bp->{left_ort}; 
	my ($s, $e);
	if($bp->{left_range} ) {
		$s = $bp->{left_range}->[0] - $l;
		$e = $bp->{left_range}->[1] + $l;
	}	
	else {
		$s = $e = $p;
		$s = $ort eq "+" ? $p - $ext : $p - $l;
		$e = $ort eq "+" ? $p + $l : $p + $ext;
	}
	$s = 1 if($s < 1);
	$e = $bp->{left_pos} + $l if($self->{type} == DEL);
	$range1 = $chr . ":" . $s . "-" . $e;

	($chr, $p) = ($bp->{right_chr}, $bp->{right_pos});
	$ort = $bp->{left_ort}; 
	if($bp->{right_range} ) {
		$s = $bp->{right_range}->[0] - $l;
		$e = $bp->{right_range}->[1] + $l;
	}	
	else {
		$s = $e = $p;
		$s = $ort eq "-" ? $p - $ext : $p - $l;
		$e = $ort eq "-" ? $p + $l : $p + $ext;
	}
	$s = $bp->{right_pos} - $l if($self->{type} == DEL);
	$s = 1 if($s < 1);
	$range2 = $chr . ":" . $s . "-" . $e;

	return ($range1, $range2);
}

sub _find_ref_range {
	my $self = shift;
	my $l = int(length($self->{consensus})/2);
	my ($range1, $range2);
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});

	my ($chr, $p) = ($bp1->{chr}, $bp1->{pos});
	my ($s, $e);
	my $ext = $l*2 ;
	$ext = 100000 if($RNASeq);
	if($self->{type} == DEL && abs($bp1->{pos} - $bp2->{pos}) < $l ){
		$l = abs($bp1->{pos} - $bp2->{pos}) - 1;
	}
	
	my $ort = $bp1->{pos} == $bp1->{left_pos} ? 
		($bp1->{left_ort} eq "+" ? "-" : "+") : $bp1->{right_ort};
	if($bp1->{range}) {
		$s = $bp1->{range}->[0] - $l;
		$s = 1 if($s < 1);
		$e = $bp1->{range}->[1] +  $l;
	}
	else {
		$s = $e = $p;
		$s = $ort eq "+" ? $p - $l  : $p - $ext;
		$e = $ort eq "+" ? $p + $ext : $p + $l;
	}
	$s = 1 if($s < 1);
	$e = $bp1->{pos} + $l if($self->{type} == DEL);
	
	$range1 = $chr . ":" . $s . "-" . $e;

	($chr, $p) = ($bp2->{chr}, $bp2->{pos});
	$ort = $bp2->{pos} == $bp2->{left_pos} ? 
		($bp2->{left_ort} eq '+' ? '-' : '+') : $bp2->{right_ort};
	if($self->{type} == INV) {
		$ort = $bp1->{pos} == $bp1->{left_pos} ?
			$bp1->{right_ort} : ($bp1->{left_ort} eq '+' ? '-' : '+');
	}

	if($bp2->{range} && $self->{type} != INV) {
		$s = $bp2->{range}->[0] - $l;
		$e = $bp2->{range}->[1] +  $l;
	}
	else {
		$s = $e = $p;
		$s = $ort eq "+" ? $p - $l : $p - $ext;
		$e = $ort eq "+" ? $p + $ext : $p + $l;
	}
	$s = $bp2->{pos} - $l if($self->{type} == DEL); 
	$s = 1 if($s < 1);

	$range2 = $chr . ":" . $s . "-" . $e;
	return($range1, $range2);
}

my $GAP = -5;
my $MATCH = 2;
my $MISMATCH = -5;
sub _refine_bp {
	my $self = shift;
	my $h = shift;
	my $ref_seq = shift;
	my $con = $self->{consensus};
	my $hsp = $h->next_hsp;
	my ($i, $j) = ($hsp->start("hit"), $hsp->start("query"));
	if($hsp->strand == -1) {
		$con = rev_comp($con);
		$j = length($con) - $hsp->end("query");
	}

	my ($score, $g, $max_score, $s) = (0, 0, 0);
	my ($current_s_r, $current_s_c) = ($i, $j);
	my ($s_r, $s_c, $e_r, $e_c);
	my $p;
	my ($n_matches, $n_max_matches) = (0, 0);
	my @bl_r = @{$hsp->gap_blocks('hit')};
	my @bl_c = @{$hsp->gap_blocks('query')};
	for( my $n = 0; $n < scalar @bl_r; $n++) {
		my $m = $bl_r[$n]->[1];
		if($p) { #gaps
			if (!$RNASeq || $bl_r[$n]->[0] - $p < 25) {
				$score += $GAP * ($bl_r[$n]->[0] - $p);
				$score = 0 if($score < 0 )
			}
			$i += $m;
			$j += $bl_c[$n]->[1];
		}

		while($m) { #matches / mismatches
			($current_s_r, $current_s_c, $n_matches) = ($i, $j, 0) 
				if($score == 0); #reset

			if(substr($ref_seq, $i, 1) eq substr($con, $j, 1)) {
				$score += $MATCH;
				$n_matches++;
			}
			else {
				$score += $MISMATCH;
			}
			if($score > $max_score) {
				$n_max_matches = $n_matches;
				$max_score = $score;
				($s_r, $s_c) = ($current_s_r, $current_s_c);
				($e_r, $e_c) = ($i, $j);
			}
			$n_matches = $score = 0 if($score < 0 );
			$m--; $i++;	$j++;
		}
		$p = $bl_r[$n]->[0] + $bl_r[$n]->[1];
	}
	return($n_max_matches, $s_r, $e_r, $s_c, $e_c);
}

sub _mapping_artifact_filter {
	my $self = shift;
	if($self->{type} == INV) {
		return ($self->_mapping_artifact($self->{consensus})
			|| $self->_mapping_artifact($self->{consensus2}) );
	}
	return $self->_mapping_artifact($self->{consensus});
}

sub _mapping_artifact {
	my $self = shift;
	my $con = shift;
	my $l = length $con;
	print STDERR "Mapping artifact filter ... ";
	$self->{error} = "Mapping Artifact";
	$l = $l - 5;
	my $options = " minIdentity=97  -out=psl -nohead";
	$options = $options . " -maxIntron=5 " unless $RNASeq;
	my $exp = $l;
	if($self->{type} == DEL && $RNASeq)	{
		$exp = $self->{second_bp}{pos} - $self->{first_bp}{pos};
		$options = $options . " -maxIntron=" . ($exp + 100) if($exp > 75000);
	}

	open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
	print $QUERY ">query\n", $con, "\n";
	close $QUERY;
	my $tmp_mapping = $mapper->run(-QUERY => "query.fa", -OPTIONS => $options);
	system("rm query.fa"); system("rm query.fa*");

	$l = $l + 5;
	if($RNASeq) {
		#can't find a decent mapping for DEL event
		return if($self->{type} == DEL && !$tmp_mapping->{query});
		foreach my $h(@{$tmp_mapping->{query}}) {
			#next if($h->{tend} - $h->{tstart} < $exp - 10);
			return if ($self->{type} != DEL && $h->{qend} - $h->{qstart} > $l - 10);
			my $chr = $h->{tchr};
			$chr = "chr" . $chr unless $chr =~ m/^chr/;
			my ($r_tree, $f_tree)  =  ($gm->sub_model($chr, "-"),  $gm->sub_model($chr, "+"));
			return unless ($f_tree && $r_tree);
			my @f_genes = $f_tree->intersect([$h->{tstart}, $h->{tend}]);
			my @r_genes = $r_tree->intersect([$h->{tstart}, $h->{tend}]);
			return if(scalar @f_genes <= 1 && scalar @r_genes <= 1);
			my @f_tmp_genes;
			my @r_tmp_genes;
			
			foreach my $g (@f_genes) {
				foreach my $block (@{$h->{blocks}}) {
					if($g->val->overlap($block)){
						push @f_tmp_genes, $g;
						last;
					}
				}
			}
			foreach my $g (@r_genes) {
				foreach my $block  (@{$h->{blocks}}) {
					if($g->val->overlap($block)){
						push @r_tmp_genes, $g;
						last;
					}
				}
			}
			return if(scalar @f_tmp_genes < 1 && scalar @r_tmp_genes < 1);
		}
	}
	else {
		foreach my $h (@{$tmp_mapping->{query}}){
			return if($h->{tend} - $h->{tstart} < $l + 10 && $h->{qend} - $h->{qstart} > $l - 5);
		}
	}
	$self->{error} = undef;
	return 1;
}

sub _germline_sclip_filter {
	my $self = shift;
	print STDERR "Germline sclip filter\n";
	return 1 if(!$sam_g);
	$self->{error} = "Germline Event";
	my $rtn = _compare_seq($self->{first_bp}) &&
		_compare_seq($self->{second_bp});
	$self->{error} = undef if($rtn);
	return $rtn;
		
}

sub _compare_seq {
	my $bp = shift;
	my $seq;
	if($bp->{pos} == $bp->{left_pos}) {
		my $seg = $sam_d->segment($bp->{right_chr}, $bp->{right_pos} -
			$germline_seq_width, $bp->{right_pos} + $germline_seq_width);
		$seq = $seg->dna;
	}
	else {
		my $seg = $sam_d->segment($bp->{left_chr}, $bp->{left_pos} - $germline_seq_width,
			$bp->{left_pos} + $germline_seq_width);
		$seq = $seg->dna;
	}
	my $sclip_seq  = _get_sclip_seqs($sam_g, $bp->{chr}, $bp->{pos} - $germline_search_width, 
		$bp->{pos} + $germline_search_width);


	open my $TARGET, ">target.fa" or croak "can't open target.fa : $OS_ERROR";
	print $TARGET ">target\n", $seq, "\n";
	close $TARGET;

	if(keys(%{$sclip_seq})) {
		open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
		foreach my $s (keys(%{$sclip_seq})) {
			print $QUERY ">", $s, "\n", $sclip_seq->{$s}, "\n";
		}
		close($QUERY);
		my $hits = $aligner->run(-TARGET => "target.fa", -QUERY => 'query.fa');
		return if(keys(%{$hits}));
	}
	return 1;
}

sub _type_distance_filter {
	my $self = shift;
	my $type = $self->{type};
	$self->{error} = "Type Distance";
	print STDERR "Type distance filter\n";
	if($type == UNK) {  #those events are not sure, always ignore them 
		print STDERR $self->to_full_string, "\n";
		return;
	}
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	return 1 if($bp2->{sc_reads} == 0);
	my ($d1, $d2);
	if($type == INS || $type == DEL) {
		($d1, $d2) = ($bp1->{left_pos} - $bp2->{left_pos}, 
			$bp1->{right_pos} - $bp2->{right_pos});
	}
	if( $type == INV ){
		($d1, $d2) = ($bp1->{left_pos} - $bp2->{left_pos}, 
			$bp2->{right_pos} - $bp1->{right_pos});
	}
	if( $type == ITX ) {
		($d1, $d2) = ($bp1->{left_pos} - $bp2->{right_pos},
			$bp2->{left_pos} - $bp1->{right_pos});
	}
	if($type == CTX ) {
		if($bp1->{left_ort} eq $bp1->{right_ort} || 
			$bp1->{sc_reads} == 0 || $bp2->{sc_reads} == 0) {
			($d1, $d2) = ($bp1->{left_pos} - $bp2->{left_pos}, $bp1->{right_pos} - $bp2->{right_pos}) ;
		}
		else {
			($d1, $d2) = ($bp1->{left_pos} - $bp2->{right_pos},	$bp2->{left_pos} - $bp1->{right_pos});
		}
	}
	print STDERR $self->to_full_string, "\t", $d1, "\t", $d2, "\n";
	if($d1 <= 1 && $d2 <= 1 || abs($d1 - $d2) <= 5) {
		$self->{error} = undef;
		return 1;
	}
	if(abs($d1) > $max_bp_dist || abs($d2) > $max_bp_dist){
		print STDERR "Distance fail\n";
		return;
	}
	$self->{error} = undef;
	return 1; # the distance is small enough
}

sub _germline_indel_filter {
	my $self = shift;
	my $type = $self->{type};
	print STDERR "Germline INDEL FILTER test\n";
	return if(abs($self->{first_bp}{pos} - $self->{second_bp}{pos}) < 40);
	return 1 if(!$sam_g);
	return 1 if( $type != INS && $type != DEL);
	return 1 if( abs($self->{first_bp}{pos} - $self->{second_bp}{pos}) > 50 );
	$self->{error} = "Germline INDEL";
	my ($start, $end ) = ($self->{first_bp}{left_pos}, $self->{first_bp}{right_pos});
	my $indel_len = abs($end - $start);
	($start, $end) = ($end, $start) if($start > $end);

	my $itr = $sam_g->features(-iterator => 1, -seq_id => $self->{first_bp}{chr},
		-start => $start, -end => $end);
	while( my $a = $itr->next_seq ) {
		next if($type == DEL && $a->cigar_str !~ m/D/);
		next if($type == INS && $a->cigar_str !~ m/I/);
		my $cigar_array = $a->cigar_array;
		my $pos = $a->pos;
		foreach my $ca (@{$cigar_array}) {
			if($ca->[0] eq 'M') {
				$pos += $ca->[1];
				next;
			}
			if($ca->[0] eq 'I' && $type == INS ){
				return  if(abs($ca->[1] - $indel_len) <= 20 && abs($pos - $start) <= 20); 
				next;
			}
			if($ca->[0] eq 'D' && $type == DEL) {
				return  if(abs($ca->[1] - $indel_len) <= 20 && abs($pos - $start) <= 20);
				$pos += $ca->[1];
				next;
			}
		}
	}
	$self->{error} = undef;
	return 1;
}

sub _tandem_repeat_filter {
	my $self = shift;
	my $con = $self->{consensus};
	print STDERR "low compexity filter\n";
	open my $QUERY, ">query.fa" or croak "can't open query.fa : $OS_ERROR";
	print $QUERY ">query\n", $self->{consensus}, "\n";
	close $QUERY;
	
	system("ptrfinder -seq query.fa -repsize $tr_min_size,$tr_max_size -minrep $tr_min_num > query.fa.rep");
	if(-e "query.fa.rep") {
		open my $REP, "<query.fa.rep" or croak "can't open query.fa.rep:$OS_ERROR";
		while( my $line = <$REP>) {
			chomp $line;
			my($pattern, $len, $times, $s, $e) = 
				$line =~ m/PATTERN\s(\w+)\sLENGTH\s(\d+)\sTIMES\s(\d+)\sSTART\s(\d+)\sSTOP\s(\d+)\sID/;
			if($self->{type} == DEL || $self->{type} == INS){
				if(abs($self->{first_bp}{left_pos} - $self->{first_bp}{right_pos}) < $tr_max_indel_size) {
					print STDERR "Tandem repeat mediated INDEL!\n";
					close $REP;
					return;
				}
			}
			else{
				if(($s < 5 || length($self->{consensus}) - $e < 5) && $len * $times > 30 ) {
					print STDERR "Tandem repeat mediated events!\n";
					close $REP;
					return;
				}
			}
		}
	}
	return 1;
}

sub _get_sclip_seqs {
	my ($sam, $chr, $start, $end) = @_;
    my %rtn;
    my $range = $chr . ":" . $start . "-" . $end;

    $sam->fetch($range, sub {
        my $a = shift;
        my $cigar_str = $a->cigar_str;
        return if($cigar_str !~ /S/);
        my ($sclip_len, $pos, $seq, $qual);
		my @cigar_array = @{$a->cigar_array};
        if($cigar_array[0]->[0] eq 'S' ) {
            $sclip_len = $cigar_array[0]->[1]; $pos = $a->start;
            return if($pos < $start or $pos > $end); # the softclipping position is not in range
 			$seq = substr($a->query->dna, 0, $sclip_len );
			$rtn{$a->qname} = $seq if($sclip_len >= 10);
		}
		#if($cigar_str =~ m/S(\d+)$/) {
		if($cigar_array[$#cigar_array]->[0] eq 'S') {
			$sclip_len = $cigar_array[$#cigar_array]->[1]; 
			$pos = $a->end;
            return if($pos < $start or $pos > $end); # the softclipping position is not in range
			$seq = substr($a->qseq, $a->l_qseq - $sclip_len );
			$rtn{$a->qname} = $seq if($sclip_len >= 10);
		}
	}
	);
	return \%rtn;
}

sub _RNASeq_strand_filter {
	my $self = shift;
	my $type = $self->type;
	print STDERR "RNASeq strand filter\n";
	return 1 unless $gm;
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	my ($ort1, $ort2) = ($bp1->{left_ort}, $bp1->{right_ort});
	my ($chr1, $chr2) = ($bp1->{left_chr}, $bp1->{right_chr});
	my ($pos1, $pos2) = ($bp1->{pos}, $bp2->{pos});
	($pos1, $pos2) = ($pos2, $pos1) if($bp1->{right_pos} == $pos1);
	my ($tmp1, $tmp2) = ($chr1 =~ m/chr/ ? $chr1 : "chr" . $chr1, 
		$chr2 =~ m/chr/ ? $chr2 : "chr" . $chr2);
	$tmp1 = "chrM" if($tmp1 eq "chrMT");
	$tmp2 = "chrM" if($tmp2 eq "chrMT");
	my ($r_tree1, $f_tree1) = ($gm->sub_model($tmp1, "-"), 	$gm->sub_model($tmp1, "+"));
	my ($r_tree2, $f_tree2) = ($gm->sub_model($tmp2, "-"), 	$gm->sub_model($tmp2, "+"));
	my ($g_ort1, $g_ort2);
	return 1 unless($r_tree1 && $r_tree2 && $f_tree1 && $f_tree2);
	$g_ort1 = "-" if(scalar($r_tree1->intersect([$pos1 - 5, $pos1])) > 0);
	$g_ort1 = "+" if(scalar($f_tree1->intersect([$pos1 - 5, $pos1])) > 0);
	$g_ort2 = "-" if(scalar($r_tree2->intersect([$pos2, $pos2 + 5])) > 0);
	$g_ort2 = "+" if(scalar($f_tree2->intersect([$pos2, $pos2 + 5])) > 0);
	return 1 unless($g_ort1 && $g_ort2);
	return 1 if($g_ort1 eq $ort1 && $g_ort2 eq $ort2);
	return 1 if($g_ort1 ne $ort1 && $g_ort2 ne $ort2);
	return;
}

# INS filter check to make sure that the INS part overlap with any gene
# if the INS part only in intron will return as a false positive
sub _RNASeq_INS_filter {
	my $self = shift;
	my $type = $self->{type};
	return 1 if($type != INS);
	print STDERR "RNAseq INS filter\n";
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	my ($pos1, $pos2) = ($bp1->{pos}, $bp2->{pos});
	($pos1, $pos2) = ($pos2, $pos1) if($pos2 < $pos1);
	my $chr = $bp1->{chr};
	$chr = "chr" . $chr unless $chr =~ m/chr/;
	my ($r_tree, $f_tree) = ($gm->sub_model($chr, "-"), 	$gm->sub_model($chr, "+"));
	foreach my $tree( ($r_tree, $f_tree) ) {
		my @tmp;
		@tmp = $f_tree->intersect([$pos1, $pos2]);
		foreach my $g (@tmp) {
			return 1 if($g->val->overlap([$pos1, $pos2]));	
		}
	}
	return;
}

sub _RNASeq_DEL_filter {
	my $self = shift;
	my $type = $self->type;
	return 1 if($type != DEL);
	print STDERR "RNAseq DEL filter\n";

	my ($gene1, $gene2, $dist) = @{$self->get_genes};
	return if($gene1 eq $gene2);
	my ($d_plus, $d_minus) = (0, 0);
	$d_plus = $1 if($dist =~ m/(\d+)\(\+\)/);
	$d_minus = $1 if($dist =~ m/(\d+)\(\-\)/);
	return if($d_plus == 0 && $d_minus == 0);
	my ($bp1, $bp2) = ($self->{first_bp}, $self->{second_bp});
	return if($bp1->{cover} == 0 || $bp2->{cover} == 0);
	my($left_len, $right_len);
	if(exists $bp1->{left_len}) {
		($left_len, $right_len ) = ($bp1->{left_len}, $bp1->{right_len});
	}
	else {
		$left_len = length($bp1->{sc_seq}) if($bp1->{sc_seq}); 
		$right_len = length($bp2->{sc_seq}) if($bp2->{sc_seq}); 
	}
	return if($left_len < 30 || $right_len < 30);
	return if($bp1->{sc_reads} < 10 || $bp2->{sc_reads} < 10);
	return unless $bp1->{left_gene}->overlap([$bp1->{pos} - $left_len, $bp1->{pos}]);
	return unless $bp1->{right_gene}->overlap([$bp2->{pos} + $right_len, $bp2->{pos}]);
	return 1;
}

sub get_statistics {
	my $self = shift;
	my $half_width = shift;
	my $sam = $self->sam_d;
	my @rtn;

	foreach my $bp ( ($self->{first_bp}, $self->{second_bp})) {
		my ($chr, $pos ) = ($bp->{chr}, $bp->{pos});
		my $range = $chr . ":" . ($pos - $half_width) . "-" . ($pos + $half_width);
		my ($n_seq, $n_rep, $total_pid) = (0, 0, 0);

		$sam->fetch($range, sub {
			my $a = shift;
			return unless ($a->start && $a->end);
			return unless ($a->start >= $pos - $half_width && $a->end <= $pos + $half_width);
			$n_seq++;
			$total_pid += _cal_pid($a);
			if($a->has_tag("XT")) {
			   $n_rep++ if($a->aux_get("XT") ne "U");
		    }
		}
		);
		if($n_seq == 0) {
			push @rtn, (0, 0);
		}
		else {
			push @rtn, ($total_pid/$n_seq, $n_rep/$n_seq);
		}
	}
	return @rtn;
}

sub _cal_pid {
	my $a = shift;
	my ($ref, $matches, $query) = $a->padded_alignment;
	my ($n_match, $n) = (0, 0);
	for( my $i = 0; $i < length($matches); $i++) {
		my $c = substr $matches, $i, 1;
		$n_match++ if($c eq "|");
		$n++;
	}
	return 0 if($n == 0);
	return $n_match/$n;
}


1;
__END__;

=pod

=head1 NAME

