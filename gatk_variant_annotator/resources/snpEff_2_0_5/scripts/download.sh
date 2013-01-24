#!/bin/sh -e

RELEASE=64

# mkdir download
cd download

#---
# Download
#---

# # Download GTF files (annotations)
# wget -r -A "*gtf.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/gtf/"
# 
# # Download FASTA files (reference genomes)
# wget -r -A "*toplevel.fa.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"
# 
# # Download CDS sequences
# wget -r -A "*cdna.all.fa.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"

#---
# Create directory structure
#---

# # Move all downloaded file to this directory
# mv `find ftp.ensembl.org -type f` .
# 
# # Gene annotations files
# for gtf in *.gtf.gz
# do
# 	short=`../scripts/file2GenomeName.pl $gtf | cut -f 5`
# 	echo ANNOTATIONS: $short
# 
# 	mkdir -p data/$short
# 	cp $gtf data/$short/genes.gtf.gz
# done
# 
# # Reference genomes files
# mkdir -p data/genomes
# for fasta in *.dna.toplevel.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 3`
# 	echo REFERENCE: $genome
# 
# 	cp $fasta data/genomes/$genome.fa.gz
# done
# 
# # CDS genomes files
# for fasta in *.cdna.all.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 3`
# 	echo CDS: $genome
# 
# 	cp $fasta data/$genome/cds.fa.gz
# done

#---
# Config file entries
#---

for fasta in *.cdna.all.fa.gz
do
	full=`../scripts/file2GenomeName.pl $fasta | cut -f 3`
	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 4`
	short=`../scripts/file2GenomeName.pl $fasta | cut -f 5`

	# 'genomes' entry
	echo -e "\t, $short \\"

	# Individual genome entry
	echo -e "# ENSEMBL v$RELEASE : $full"
	echo -e "$short.genome : $genome"
	echo
done

# Back to parent dir
cd - > /dev/null

#---
# Create build queue entries
#---

rm -vf queue_build.txt

# Build from TXT files
for genes in data/*/genes.txt*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "./scripts/snpEffXL.sh build -v $genomeName"
done | sort >> queue_build.txt

# Build from GFF2 files
echo "./scripts/snpEffXL.sh build -v -gff2 amel2"

# Build from GFF3 files
for genes in `ls data/*/genes.gff* | grep -v amel2`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "./scripts/snpEffXL.sh build -v -gff3 $genomeName"
done | sort >> queue_build.txt

# Build from GTF22 files
for genes in data/*/genes.gtf*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "./scripts/snpEffXL.sh build -v -gtf22 $genomeName"
done | sort >> queue_build.txt

#---
# Create test queue entries
#---

# Build from GTF22 files
#for cds in data/*/cds.*
#do
#	dir=`dirname $cds`
#	genomeName=`basename $dir`
#	echo "./scripts/snpEffXL.sh cds -v $genomeName $cds"
#done

