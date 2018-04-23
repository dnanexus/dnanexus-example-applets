package SVUtil;

use strict;
use Carp;
use English;
use Bio::SeqIO;
use Bio::DB::Sam;
use File::Spec;
use File::Path qw(remove_tree);
use File::Temp qw/ tempfile tempdir /;
use Cwd;
use SCValidator qw(LEFT_CLIP RIGHT_CLIP);
require Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(is_new_pair clip_fa_file prepare_reads parse_range
	get_input_bam get_work_dir is_PCR_dup read_fa_file find_smallest_cover
	get_sclip_reads rev_comp);
our @EXPORT_OK = qw($use_scratch $spacer_seq $spcer_seq_len $work_dir);
our $use_scratch = 1;
our  $work_dir;

sub rev_comp {
    my $str = shift;
    $str = reverse $str;
    $str =~ tr/ACTG/TGAC/;
    return $str;
}
# something about binormial distribution
sub choose {
    my ($n,$k) = @_;
    my ($result,$j) = (1,1);
	    
    return 0 if $k>$n||$k<0;
    $k = ($n - $k) if ($n - $k) <$k;
	    
    while ($j <= $k ) {
        $result *= $n--;
        $result /= $j++;
    }
    return $result;
}

sub dbinorm { #desity funtion 
	my($n, $k, $p) = @_;
	return unless($n > 0 && $k >= 0 && $n >= $k && $p >= 0 && $p <=1);
	my $prob = log(choose($n, $k)) + ($n-$k)*log(1-$p) + $k*log($p);
	return exp($prob);
}

sub find_smallest_cover {
	my ($c, $p, $crit) = @_;
	my @s_cover;
	for( my $i = 1; $i < $c; $i++) {
		my $n = 1;
		while(1) {
			my $tmp = 0;
			for( my $j = 0; $j <= $i && $j < $n; $j++) {
				$tmp += dbinorm($n, $j, $p);
			}
			if( $tmp < $crit) {
				$s_cover[$i] = $n - 1;
				last;
			}
			$n++;
		}
	}
	return \@s_cover;
}



sub read_fa_file {
	my $file = shift;
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	my %seqs;
	while( my $seq=$in->next_seq()) {
		$seqs{$seq->display_id} = $seq->seq;
	}
	return \%seqs;
}

# rmdup: remove PCR dumplication related code
# you can replace it using your own criteria 
sub is_new_pair {
    my ($a, $sclip_len, $pairs) = @_;
    foreach my $p (@{$pairs}) {
        return undef if($a->start == $p->[0] && $a->end == $p->[1] &&
            $a->mate_start == $p->[2] && $a->mate_end == $p->[3] && $sclip_len == $p->[4]);
    }
    return 1;
}

#
#clipping sequnce related code
#
# the spacer sequence, you can replace it with your own one
my $spacer_seq = 'ACGTACGT' . 'CGTA' x 6 . 'ATGCATGC';
my $spacer_seq_len = length($spacer_seq);

# clip the reads that with spacer_seq at one end
sub clip_fa_file {
	my ($file, $ort) = @_;
	my $out = $file . ".clip.fa";
	my $in  = Bio::SeqIO->new( -file =>	$file, -format => 'Fasta');
	open my  $OUT, ">$out" or croak "can't open $file.fa:$OS_ERROR";
	while( my $seq=$in->next_seq()) {
		print $OUT ">",  $seq->display_id, "\n";
		if($ort eq RIGHT_CLIP) {
			print $OUT  $seq->subseq($spacer_seq_len + 1, $seq->length), "\n";
		}
		else {
			print $OUT  $seq->subseq(1, $seq->length - $spacer_seq_len), "\n";
		}
	}
	close($OUT);
	return $out;
}

# prepare the reads into a fasta file
# due to most assembler can't handle short read well, we add
# the spacer_seq to the start/end of each sequence arrording to
# the softclip postion
sub prepare_reads {
	my($r_reads, $clip) = @_;
	my ($fh, $fa_name) = tempfile( DIR => $work_dir, SUFFIX => '.fa');
	my $qual_name = $fa_name . ".qual";
	my $spacer_qual = ' 50 ' x $spacer_seq_len;

	open my $qual_fh, ">$qual_name" or croak "can't open $qual_name:$OS_ERROR";
	foreach my $r (@{$r_reads}) {
		print $fh ">", $r->{name}, "\n";
		my $tmp = $r->{seq};
		$tmp = defined $clip ? 
			($clip == RIGHT_CLIP ? $spacer_seq . $tmp : $tmp . $spacer_seq) :
			$tmp;
		print $fh $tmp, "\n";
		next unless $r->{qual};
		print $qual_fh ">", $r->{name}, "\n";
		$tmp =  join(" ", @{$r->{qual}});
		$tmp = defined $clip ? 
			($clip == RIGHT_CLIP ? $spacer_qual . $tmp : $tmp . $spacer_qual) :
			$tmp;
		print $qual_fh $tmp, "\n";
	}
	close $qual_fh;
	return $fa_name;
}

# check PCR duplicate
sub is_PCR_dup {
	my ($a, $pairs, $sclip_len) = @_;
	my ($s, $e, $ms, $me ) = ($a->start, $a->end, $a->mate_start ? $a->mate_start : 0,
		$a->mate_end ? $a->mate_end : 0);
    foreach my $p (@{$pairs}) {
		return 1 if($s == $p->[0] && $e == $p->[1] &&
			$ms == $p->[2] && $me == $p->[3] && $sclip_len == $p->[4]);
	}
	return;
}

# parse the range of input, format is: chr:start-end
# start and end is optional
sub parse_range {
	my $range = shift;
	my ($chr, $start, $end);
	my @field = split /:/, $range;
	$chr = $field[0];
#	$chr = substr($chr, 3) if($chr =~ /^chr/);
	if(@field > 1) {
		@field = split /-/, $field[1];
		$start = $field[0];
		$end = $field[1] if($field[1]);
		if($start !~ m/^\d+$/ or $end !~ m/^\d+$/) {
			croak "wrong range format, need to be: chr:start-end";
		}
	}
	return ($chr, $start, $end);
}

# St. Jude only
# grab the bam file from the bam dir
my $raw_bam_dir = "/nfs_exports/genomes/1/PCGP/BucketRaw";
sub get_input_bam {
    my $sample = shift;
	my ($group) = $sample =~ /(.*?)\d+/; 	
    my $tmp_bam_dir = File::Spec->catdir($raw_bam_dir, $group);
    opendir(my $dh, $tmp_bam_dir) or croak "can't open directory: $raw_bam_dir: $OS_ERROR";
    my @files = grep { /^$sample-.*bam$/ } readdir($dh);
    close $dh;
    return File::Spec->catfile($raw_bam_dir, $group,  $files[0]);
}

#scratch related code
sub get_work_dir {
    my %options = @_;
	my $scratch_dir;
	if($options{-SCRATCH}){
    	$scratch_dir = (exists $options{-DIR} ? $options{-DIR} :
        	(-e '/scratch_space' ? '/scratch_space' : undef) );
	}

    my $dir;
	if($scratch_dir) {
		$scratch_dir = File::Spec->rel2abs($scratch_dir);
	    if($ENV{'JOB_ID'}) {
    	    $dir = File::Spec->catdir($scratch_dir, $ENV{'JOB_ID'})
	    }
		else {
			$dir = tempdir(DIR => $scratch_dir);
		}
	}
    else {
			$dir = tempdir();
    }
	mkdir($dir) if( ! -e $dir);
	$work_dir = $dir;
    return $dir;
}


sub get_sclip_reads {
	my %args = @_;
	my ($sam, $chr, $start, $end, $clip, $min_len, $validator, $paired, $rmdup, $full_seq, $extra_base) = 
		($args{-SAM}, $args{-CHR}, $args{-START}, $args{-END}, $args{-CLIP}, 
		$args{-MINCLIPLEN}, $args{-VALIDATOR}, $args{-PAIRED}, $args{-RMDUP},
		$args{-FULLSEQ}, $args{-EXTRA});
	$extra_base = 0 unless $extra_base;
	my @rtn;
	$min_len = 0 unless $min_len;
	my $max_sclip_len = 0;
	my $range = $chr;
	$range = $chr . ":" . $start . "-" . $end if($start && $end);
	my @pairs;
	my $strand_validator = $paired ? 1 : 0;

	$sam->fetch($range, sub {
		my $a = shift;
		my $cigar_str = $a->cigar_str;
		my @cigar_array = @{$a->cigar_array};
		return if($cigar_str !~ m/S/);
		return unless($a->start && $a->end);	
		#return if(!$a->proper_pair && $paired);
		if($paired && !$strand_validator) {
			$validator->add_validator("strand_validator");
		}

		if($paired && !$a->proper_pair) {
			$validator->remove_validator("strand_validator");
			$strand_validator = 0;
		}
        my ($sclip_len, $pos, $seq, $qual);
		my @tmp = $a->qscore;
		if($cigar_array[0]->[0] eq 'S' && $clip == LEFT_CLIP ) {
			$pos = $a->start;
			return if($pos < $start or $pos > $end); # the softclipping position is not in range
            $sclip_len = $cigar_array[0]->[1]; 
			#print $cigar_str, "\t", scalar @pairs, "\n";
			return unless($validator &&  $validator->validate($a, LEFT_CLIP));
			$max_sclip_len = $sclip_len if($max_sclip_len < $sclip_len);
			return if(@pairs > 0 && $rmdup && is_PCR_dup($a, \@pairs, $sclip_len)); #it's a PCR dumplicate
            $seq = $full_seq ? $a->query->dna :  substr($a->query->dna, 0, $sclip_len + $extra_base);
			@tmp = $full_seq ? @tmp :  @tmp[0 .. ($sclip_len -1 + $extra_base)];

			push @pairs, [$a->start, $a->end, $a->mate_start ? $a->mate_start : 0, 
				$a->mate_end ? $a->mate_end : 0, $sclip_len] if($rmdup);
			my $qscore = $a->qscore;
				
			push @rtn, {
				name => $a->qname,
				seq => $seq,
				qual => \@tmp,
				full_seq => $a->query->dna,
				full_qual => $qscore,
				pos => $pos,
			};
		}
		if($cigar_array[$#cigar_array]->[0] eq 'S' &&  $clip == RIGHT_CLIP ) { 
		  	$pos = $a->end;
			return if($pos < $start or $pos > $end); # the softclipping position is not in range
#			print $cigar_str, "\t", scalar @pairs, "\n";
			return unless $validator->validate($a, RIGHT_CLIP);
            $sclip_len = $cigar_array[$#cigar_array]->[1]; 
			$max_sclip_len = $sclip_len if($max_sclip_len < $sclip_len);
			return if(@pairs > 0 && $rmdup && is_PCR_dup($a, \@pairs, $sclip_len)); #it's a PCR dumplicate
            $seq = $full_seq ? $a->qseq : substr($a->qseq, $a->l_qseq - $sclip_len - $extra_base );
			@tmp = $full_seq ? @tmp : @tmp[($a->l_qseq - $sclip_len - $extra_base) .. ($a->l_qseq - 1)];
			push @pairs, [$a->start, $a->end, $a->mate_start ? $a->mate_start : 0, 
				$a->mate_end ? $a->mate_end : 0, $sclip_len] if($rmdup);
			my $qscore = $a->qscore;
			push @rtn, {
				name => $a->qname,
				seq => $seq,
				qual => \@tmp,
				full_seq => $a->query->dna,
				full_qual => $qscore,
				pos => $pos,
			};
        }
	}
	);
	if($max_sclip_len >= $min_len) {
		return @rtn;
	}
	else{
		undef @rtn; return @rtn;
	}
}

my $start_dir;
BEGIN {
	$start_dir = getcwd;
}
END {
	chdir($start_dir);
    remove_tree($work_dir) if($use_scratch && $work_dir && $start_dir ne $work_dir);
}

1;
