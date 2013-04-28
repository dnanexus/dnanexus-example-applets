PV=2.17.0
PN=BEDTools.v$(PV)
all: bedtools
	cp bedtools-$(PV)/bin/bedtools resources/usr/bin/
	./generate-bedtools-applets

bedtools:
	curl -O "https://bedtools.googlecode.com/files/$(PN).tar.gz"
	tar xzf $(PN).tar.gz
	(cd bedtools-$(PV); make)
