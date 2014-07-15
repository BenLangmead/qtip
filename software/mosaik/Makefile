#!gmake

VER=2.2.3
NM=MOSAIK-$(VER)
AR=$(NM)-source.tar
URL=https://mosaik-aligner.googlecode.com/files/$(AR)

MosaikBuild: $(NM)/bin/MosaikBuild
	ln -s -f $< $@

MosaikAligner: $(NM)/bin/MosaikAligner
	ln -s -f $< $@

$(NM)/mosaik: $(NM)
	cd $(NM) && $(MAKE)

$(NM): $(AR)
	tar xvf $<
	mv $(NM)-source $@

$(AR):
	wget --no-check $(URL)

.PHONY: clean
clean:
	rm -f bwa
	rm -rf bwa-*.tar.bz2
	rm -rf bwa-*
