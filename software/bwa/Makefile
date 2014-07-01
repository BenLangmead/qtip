#!gmake

VER=0.7.9a
NM=bwa-$(VER)
AR=$(NM).tar.bz2
URL=http://downloads.sourceforge.net/project/bio-bwa/$(AR)

bwa: $(NM)/bwa
	ln -s -f $< $@

$(NM)/bwa: $(NM)
	make -C $(NM)

$(NM): $(AR)
	tar jxvf $<

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -f bwa
	rm -rf bwa-*.tar.bz2
	rm -rf bwa-*
