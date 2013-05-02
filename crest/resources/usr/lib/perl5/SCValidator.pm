package SCValidator;
# this is a singleton class

use strict;
use Exporter 'import';
our @EXPORT_OK = qw($lowqual_cutoff $min_percent_id $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use Bio::DB::Sam;
use constant LEFT_CLIP => 0;
use constant RIGHT_CLIP => 1;

our $min_percent_hq = 80;
our $lowqual_cutoff = 20;
our $min_percent_id = 90;

my $SCValidator;
my %default_validators = ( 
	"quality_validator" => \&_quality_validator,
	"align_validator"	=> \&_align_validator,
	"strand_validator"	=> \&_strand_validator,
);

sub new {
	if( $SCValidator) {
		return $SCValidator;
	} else {
		my $type = shift;
		my $this = {
			"validator" => {}	
		};
		foreach my $v ( keys(%default_validators)) {
			$this->{validator}{$v} = $default_validators{$v};
		}
		$SCValidator =  bless $this, $type;
	}
}

sub add_validator {
	my $self = shift;
	my $name = shift;
	$name = lc $name;
	if( exists($self->{validator}{$name})) {
#		print STDERR "Validator: $name is already in the list";
		return;
	}
	if( exists($default_validators{$name})) {
		$self->{validator}{$name} = $default_validators{$name};
		return;
	}
	my $validator = shift;
	die "Validator must be coderef" if(ref($validator) ne "CODE");
	$self->{validator}{$name} = $validator;
}

sub remove_validator {
	my ($self, $name) = @_;
	$name = lc $name;
	if(exists($self->{validator}{$name})) {
		delete $self->{validator}{$name};
	}
}

sub validate {
	my $self = shift;
	my $a = shift;
	my $clip = shift;
	die "The method will only validate bam alignment object"
		if( ref($a) ne "Bio::DB::Bam::Aligment" && 
			ref($a) ne "Bio::DB::Bam::AlignWrapper" );
	return undef if($a->unmapped);
	foreach my $v (keys(%{$self->{"validator"}})) {
		return undef if(!$self->{"validator"}{$v}->($a, $clip));
	}
	return 1;
}

sub _align_validator {
	#we only consider high quality part of the alignment
	my $a = shift;
	my ($ref, $matches, $query) = $a->padded_alignment;
	my @qscores = $a->qscore;
	my $n_matches = 0;
	my $align_len = 0;
	my $j = 0;
	my ($start, $end) = (0, length($ref)-1);
	my @cigar_array = @{ $a->cigar_array };
	$start = $cigar_array[0]->[1] if($cigar_array[0]->[0] eq 'S');
	$end = $end - $cigar_array[$#cigar_array]->[1] if($cigar_array[$#cigar_array]->[0] eq 'S');
	
	for( my $i = $start; $i <= $end; $i++) {
		my $m = substr($matches, $i, 1);
		if($m ne '|') {
			if( substr($query, $i, 1) ne '-') {
		#		$align_len++ if($qscores[$j] >= $lowqual_cutoff);
				$j++;
			}
		#	else {
		#		$align_len++ if($qscores[$j] >= $lowqual_cutoff && $qscores[$j+1] >= $lowqual_cutoff);
		#	}
			next;
		}
		if($qscores[$j] >= $lowqual_cutoff) {
			$n_matches++ if(substr($ref, $i, 1) eq substr($query, $i, 1));
			$align_len++;
		}
		$j++;
	}
	return ($align_len - $n_matches > 1) ? undef : 1 if($align_len < 20);
	return ($n_matches * 100 >= $align_len * $min_percent_id) ? 1 : undef;
}

sub _quality_validator {
	my $a = shift;
	my $clip = shift;
#	my $cigar_str = $a->cigar_str;
	my @cigar_array = @{$a->cigar_array};
	my @qual = $a->qscore;
	if($clip == LEFT_CLIP) {
		#my($sclip_len) = $cigar_str =~ m/^S(\d+)/;
		my $sclip_len = $cigar_array[0]->[1];
		@qual = @qual[0 .. ($sclip_len - 1)];
	}
	elsif ($clip == RIGHT_CLIP) {
		#my($sclip_len) = $cigar_str =~ m/S(\d+)$/;
		my $sclip_len = $cigar_array[$#cigar_array]->[1];
		@qual = @qual[($a->l_qseq - $sclip_len) .. ($a->l_qseq - 1)];
	}
	else {
		print STDERR "Please specify the soft_clipping position:
			0: LEFT_CLIP or 1: RIGHT_CLIP";
	}

	my $n_hq = 0;
	for( my $i = 0; $i <= $#qual; $i++) {
		$n_hq++ if($qual[$i] >= $lowqual_cutoff);
	}
	return ($n_hq * 100 > $#qual * $min_percent_hq) ? 1 : undef;
}

sub _strand_validator {
	my $a = shift;
	my $clip = shift;
	if($clip == LEFT_CLIP ) {
		return ($a->paired && ($a->mtid != $a->tid || $a->mate_start < $a->start)) ? undef : 1;
	}
	elsif($clip == RIGHT_CLIP) {
		return ($a->paired && ($a->mtid != $a->tid || $a->mate_start > $a->end)) ? undef : 1;
	}
	else {
		 print STDERR "Please specify the soft_clipping position:
             0: LEFT_CLIP or 1: RIGHT_CLIP";
	}
}

1;
