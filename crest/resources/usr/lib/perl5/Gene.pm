package Gene;
use strict;
use Transcript;
use Carp;
use Data::Dumper;

# a light weight gene structure is used here

my @Gene_slots;
BEGIN {
	@Gene_slots = qw(NAME CHR START END STRAND EXONS TRANSCRIPTS TYPE);
}
use enum @Gene_slots;

my %attribute = (
    name         => NAME,
	chr			 => CHR,
    start        => START,
    end          => END,
	strand		 => STRAND,
	transcripts  => TRANSCRIPTS,
	type		 => TYPE,
);

#using an array instead of a hash for the node
sub _accessor {
    my $index = shift;
    return sub {
        my $self = shift;
        return undef unless $self;
        if (@_) {
          $self->[$index] = shift;
        }
        return $self->[$index];
    };
}

while(my($at, $idx) = each %attribute) {
    no strict 'refs';
    *$at = _accessor($idx);
}

sub new {
    my $class = shift;
    my $obj = [];
	$obj->[TRANSCRIPTS] = [];
    if (@_) {
		my %arg = @_;
        $obj->[NAME]        = $arg{-NAME} if($arg{-NAME});
        $obj->[CHR]         = $arg{-CHR} if($arg{-CHR});
		$obj->[START]       = $arg{-START} if($arg{-START});
		$obj->[END]         = $arg{-END} if($arg{-END});
		$obj->[STRAND]      = $arg{-STRAND} if($arg{-STRAND});
		$obj->[TRANSCRIPTS] = $arg{-TRANSCRIPTS} if($arg{-TRANSCRIPTS});
    }
    return bless $obj, $class;
}

sub add_transcript {
	my ($self, $fea) = @_;
	croak "You must add a Transcript type into a gene" 
		 unless ($fea->isa('Transcript'));
	if($self->[STRAND] && $self->[STRAND] ne $fea->strand) {
		croak "The transcript has different orientation with the gene";
	}
	if($self->[CHR] && $self->[CHR] ne $fea->chr) {
		croak "The transcript is on different chr with the gene";
	}
#	if($self->[TYPE] && $fea->type ne $fea->type) {
#		croak "The type of the transcript are different from the gene";
#	}
	$self->[STRAND] = $fea->strand;
	$self->[CHR] = $fea->chr;
#	$self->[TYPE] = $fea->type;

	push @{$self->[TRANSCRIPTS]}, $fea;

	$self->[NAME] = $self->[NAME] ? $self->[NAME] . "," . $fea->name : $fea->name;

	
	#update the start and end of the gene
	$self->[START] = $fea->start if(!$self->[START] || $self->[START] > $fea->start);
	$self->[END]   = $fea->end if(!$self->[END] || $self->[END] < $fea->end);
}

sub get_start {
	my ($self, $pos, $ext) = @_;
	my $rtn = $pos;
	foreach my $t (@{$self->[TRANSCRIPTS]}) {
		my $tmp = $t->get_start($pos, $ext);
		$rtn = $tmp if($tmp < $rtn);
	}
	return $rtn;
}

sub get_end {
	my ($self, $pos, $ext) = @_;
	my $rtn = $pos;
	foreach my $t (@{$self->[TRANSCRIPTS]}) {
		my $tmp = $t->get_end($pos, $ext);
		$rtn = $tmp if($tmp > $rtn);
	}
	return $rtn;
}

sub overlap {
	my ($self, $fea) = @_;

	if(ref($fea) eq 'ARRAY') { 
		foreach my $t ( @{$self->[TRANSCRIPTS]} ) {
			return 1 if($t->overlap($fea));
		}
	}
	elsif($fea->isa('Transcript')) {
		return if($fea->strand &&  $self->[STRAND] ne $fea->strand );
		return if($fea->chr && $self->[CHR] ne $fea->chr) ;
		#return if($fea->type && $self->[TYPE] ne $fea->type);
		foreach my $e ( @{$fea->exons}) {
			foreach my $t ( @{$self->[TRANSCRIPTS]} ) {
				return 1 if($t->overlap($e));
			}
		}
	}
	else {
		croak "Not implemented overlap";
	}
	return 0;
}

1;
