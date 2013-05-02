SHELL=/bin/bash -e

VERSION := ${shell cat dxapp.json | grep upstreamVersion | sed 's/[^0-9.-]*\([0-9.-]*\).*/\1/'}
CIRCOS=circos-$(VERSION)

all: circos

circos:
	wget http://circos.ca/distribution/$(CIRCOS).tgz
	tar -xzf $(CIRCOS).tgz -C resources/opt
	rm $(CIRCOS).tgz

clean:
	rm -r resources/opt/$(CIRCOS)
