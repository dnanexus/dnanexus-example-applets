#!/nfs_exports/apps/64-bit/gnu-apps/perl5.8.9/bin/perl -w
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
use Cwd;

# bam related variables
# input/output
my($file_d, $file_g);

my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	'd|diagnostic=s'	=> \$file_d,
	'g|germline=s'		=> \$file_g,
	'h|help|?'			=> \$help,
	'man'				=> \$man,
	'usage'				=> \$usage,
	'v|version'			=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

my $start_dir = getcwd;

# figure out input file
if(!$file_d or !$file_g ) {
	croak "you need specify diagnostic and germline input file";
}
my ( @cover_D, @cover_diff);
for( my $i = 1; $i < 1000; $i++) {
	$cover_D[$i] = 0;
	$cover_diff[$i] = 0;
}
my($prefix, $path, $suffix) = fileparse($file_d);
my $diff_file = $prefix . ".somatic.cover";
$diff_file = File::Spec->catfile($path, $diff_file);
open my $SOMATIC, ">$diff_file" or croak "can't open $diff_file:$OS_ERROR";
open my $IN_D, "<$file_d" or croak "can't open $file_d:$OS_ERROR";
open my $IN_G, "<$file_g" or croak "can't open $file_g:$OS_ERROR";
my %g_sclip;
while( my $line = <$IN_G> ) {
	chomp $line;
	my ($chr, $pos, $ort, $cover) = split /\t/, $line;
	$g_sclip{$chr} = {} if(!exists($g_sclip{$chr}));
	$g_sclip{$chr}->{$pos} = $ort;
}
close($IN_G);
while( my $line = <$IN_D> ) {
	chomp $line;
	my ($chr, $pos, $ort, $cover) = split /\t/, $line;
	$cover_D[$cover]++;
	next if(exists ($g_sclip{$chr}->{$pos}) && $g_sclip{$chr}->{$pos} eq $ort);
	$cover_diff[$cover]++;
	print $SOMATIC $line, "\n";
}
close($IN_D);
close $SOMATIC;
for( my $i = 1; $i < 1000; $i++) {
	print join("\t", $i, $cover_D[$i], $cover_diff[$i]), "\n";
}

