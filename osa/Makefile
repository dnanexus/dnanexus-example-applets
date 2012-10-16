all:
	(cd resources/opt; wget -O OSAv1.8.1.zip "http://www.omicsoft.com/osa/software/OSAv1.8.1.zip")
	(cd resources/opt; unzip -o OSAv1.8.1.zip)
	(cd resources/opt; rm -rf OSAv1.8.1/TestData1/output_linux/*)
	(cd resources/opt/OSAv1.8.1; mono osa.exe --alignrna . Human.B37 RefGene TestData1/secontrol_linux.ini)
