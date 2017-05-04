package SVExtToolsI;

# this package has three major tools that the program is going to use:
# 1. Assembler, cap3 is used for this purpose
# 2. Mapper, blat server version is used and highly recommended
# 3. Aligner, any one will work, but the program need to parse the output 
sub new {
	my $class = shift;
	my %param = @_;
	my $self = {};
	$self->{PRG} = undef;
	$self->{OPTIONS} = undef;
	foreach my $k (keys %param) {
		my $k1 = $k;
		$k1 = substr($k, 1) if($k =~ m/^-/);
		$self->{$k1} = $param{$k};
	}

	bless $self, ref($class) || $class;
	return $self;
}

sub program {
	my ($self, $value) = @_;
	$self->{PRG} = $value if($value);
	return $self->{PRG};
}

sub options {
	my ($self, $value) = @_;
	$self->{OPTIONS} = $value if($value);
	return $self->{OPTIONS};
}

sub run {
	print STDERR "you need implement your own run method";
}

package Assembler;
use strict;
use Carp;
use English;
use base qw(SVExtToolsI);

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{PRG} = "cap3" if(!$self->{PRG});
	$self->{OPTIONS} = "" if(!$self->{OPTIONS});
	return $self;
}

# the run function of Assembler returns the contig file and 
# a hashref with the number of reads in each contig

sub run {
	my($self, $file) = @_;
	croak "$self->{PRG} need an input fasta file to do the assembly" if(!$file);
	if(-s $file == 0) {
		print STDERR "$file is of size 0";
		return;
	}
	system(join(" ", ($self->{PRG}, $file, $self->{OPTIONS})));
	my( $r_count, $r_reads ) = _count_reads("$file.cap.ace");
	return ("$file.cap.contigs", $r_count, $r_reads, "$file.cap.singlets");
}

sub _count_reads {
	my $file = shift;
	my %count;
	my %reads;
    open my $ACE, "<$file" or croak "Can't open $file:$OS_ERROR";
	my $contig_name;
    while( my $line = <$ACE> ) {
        if($line =~ m/^CO\s(.*?)\s\d+\s(\d+)/){
			$contig_name = $1;
			$count{$contig_name} = $2;
			$reads{$contig_name} = [];
		}
		if($line =~ m/^RD\s(.*?)\s/) {
			push @{$reads{$contig_name}}, $1;	
		}
	}
	close($ACE);
	return (\%count, \%reads);
}

package Mapper;
use strict;
use Carp;
use Data::Dumper;
use English;
use Bio::SearchIO;
use Bio::SeqIO;
use base qw(SVExtToolsI);

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{PRG} = "gfClient sjblat 50000" if(!$self->{PRG});
	$self->{OPTIONS} = "-out=psl -nohead" if(!$self->{OPTIONS});
	$self->{MAX_SCORE_DIFF} = 20 if(!$self->{MAX_SCORE_DIFF});
	$self->{MAX_NUM_HITS} = 10 if(!$self->{MAX_NUM_HITS});
	return $self;
}

sub max_score_diff {
	my $self = shift;
	my $value = shift;
	$self->{MAX_SCORE_DIFF} = $value if($value);
	return $self->{MAX_SCORE_DIFF};
}

sub max_num_hits {
	my $self = shift;
	my $value = shift;
	$self->{MAX_NUM_HITS} = $value if($value);
	return $self->{MAX_NUM_HITS};
}

sub run {
	my $self = shift;
	my %param = @_;
	croak "Missing QUERY parameter for $self->{PRG}" if(!$param{-QUERY});
	my $output = $param{-OUTPUT} || $param{-QUERY} . ".psl";
	my $options = $param{-OPTIONS} || $self->{OPTIONS};
	system(join(" ", ($self->{PRG}, $self->{BIT2_DIR}, $param{-QUERY}, $output, $options)));
	system("sort -k 10g -k 1nr $output -o $output.sorted");
	open my $SORTED, "<$output.sorted" or croak "can't open $output.sorted:$OS_ERROR";
	open my $PSL, ">$output" or croak "can't open $output:$OS_ERROR";
	my $n = 0;
	while( my $line = <$SORTED>) {
		print $PSL $line;
		$n++;
		last if($n > 10);
	}
	close($PSL);
	close($SORTED);
	return $self->select_target($output);
}


sub select_target {
	my $self = shift;
	my $file = shift;
	my %targets;
	my $parser = Bio::SearchIO->new( -file => $file, -format => 'psl');
	while( my $result = $parser->next_result ) {
		my $qname = $result->query_name;
		$targets{$qname} = [];
	#	$result->sort_hits(sub {$Bio::Search::Result::ResultI::b -> matches('id') <=> 
	#	            $Bio::Search::Result::ResultI::a ->matches('id')});
		my $n_hits = 0;
		my $max_score = 0;
#		my $perfect_only;
		while( my $hit = $result->next_hit ) {
			while( my $hsp = $hit->next_hsp) {
				my ($n_matches) = $hsp->matches;
				last if($max_score - $n_matches > $self->{MAX_SCORE_DIFF});
				#last if($perfect_only && $hit->query_length != $n_matches);
				my @blocks;
				foreach my $bl (@{$hsp->gap_blocks('hit')}) {
					push @blocks, [$bl->[0], $bl->[0] + $bl->[1]];
				}
				push @{$targets{$qname}}, {
					tchr	=> $hit->name,
					tstart	=> $hsp->start('hit'),
					tend	=> $hsp->end('hit'),
					qstart	=> $hsp->start('query'),
					qend	=> $hsp->end('query'),
					qstrand => $hsp->strand('query') == 1 ? '+' : '-',
					matches	=> $n_matches,
					blocks	=> \@blocks,
					perfect => $hit->query_length == $n_matches ? 1 : 0,
				};
				#$perfect_only = 1 if($hit->query_length == $n_matches);
				$max_score = $n_matches if($max_score < $n_matches);
				$n_hits++;
			}
			$hit->rewind;
			last if($n_hits >= $self->{MAX_NUM_HITS});
		}
	}
	return \%targets;
}

package Aligner;
use strict;
use Carp;
use English;
use Bio::SearchIO;
use Bio::SeqIO;
use base qw(SVExtToolsI);

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{PRG} = "blat" if(!$self->{PRG});
	$self->{OPTIONS} = "-tileSize=9 -stepSize=1 -out=psl -nohead -minScore=15" 
		if(!$self->{OPTIONS});
	return $self;
}

# return best hit for each query reads
sub run {
	my $self = shift;
	my %param = @_;
	croak "Missing TARGET parameter for $self->{PRG}" if(!$param{-TARGET});
	croak "Missing QUERY parameter for $self->{PRG}" if(!$param{-QUERY});
	my $output = $param{-OUTPUT} || $param{-QUERY} . ".psl";
	system( join(" ", ($self->{PRG}, $param{-TARGET}, $param{-QUERY}, $output, $self->{OPTIONS})));
	return $output if($param{-FILE});
	return _find_best_hit($output);
}

sub _find_best_hit {
	my $file = shift;
	my $parser = Bio::SearchIO->new( -file => $file, -format => 'psl');
	my %best_hit;
	while( my $result = $parser->next_result ) {
		$result->sort_hits(sub {$Bio::Search::Result::ResultI::b -> matches('id') <=> 
		            $Bio::Search::Result::ResultI::a ->matches('id')});		
		my $hit = $result->next_hit; # the best hit
		my $hsp = $hit->next_hsp;
		$best_hit{$result->query_name} = $hit;
		$hit->rewind;
	}
	return \%best_hit;
}

1;
