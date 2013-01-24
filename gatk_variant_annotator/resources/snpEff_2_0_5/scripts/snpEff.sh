#!/bin/sh

DIR=$HOME/snpEff/
LIB=$HOME/snpEff/lib

java -Xmx1G \
	-classpath "$LIB/charts4j-1.2.jar:$LIB/flanagan.jar:$LIB/freemarker.jar:$LIB/junit.jar:$LIB/trove-2.1.0.jar:$DIR" \
	ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff \
	$*
