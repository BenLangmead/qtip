#!gmake

VER=2.2.3
NM=MOSAIK-$(VER)
AR=$(NM)-source.tar
URL=https://mosaik-aligner.googlecode.com/files/$(AR)

mosaik: $(NM)/mosaik
	ln -s -f $< $@

$(NM)/mosaik: $(NM)
	make -C $(NM)

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
