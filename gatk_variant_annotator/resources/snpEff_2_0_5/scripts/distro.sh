#!/bin/sh

#------------------------------------------------------------------------------
# Create a zip file for distribution
# Note: Only binary data is included (no raw gene info / genomes)
#
#                                      Pablo Cingolani 2010
#------------------------------------------------------------------------------

VERSION="2_0_5"
DIR=$HOME/snpEff_$VERSION
rm -rvf $DIR
mkdir $DIR

# Copy core files
cp snpEff.config snpEff.jar $DIR
cp -rvfH galaxy scripts $DIR

cd $DIR
rm -rvf `find . -name "CVS" -type d`
cd -

# Create 'core' zip file
cd $HOME
ZIP="snpEff_v"$VERSION"_core.zip"
rm -f $ZIP 2> /dev/null
zip -r $ZIP snpEff_$VERSION
cd -

# Create ZIP file for each database
for d in `ls data/*/snpEffectPredictor.bin`
do
	DIR=`dirname $d`
	GEN=`basename $DIR`
	
	echo $GEN
	ZIP="snpEff_v"$VERSION"_"$GEN".zip"
	zip -r $ZIP data/$GEN/*.bin
done
