#!gmake

VER=0.7.12
NM=bwa-$(VER)
AR=$(NM).tar.bz2
URL=http://downloads.sourceforge.net/project/bio-bwa/$(AR)

bwa: $(NM)/bwa
	ln -s -f $< $@

$(NM)/bwa: $(NM)
	make -C $(NM)

$(NM): $(AR) bwa_conc_flags.patch
	tar jxvf $<
	cd $@ && patch -p1 < ../bwa_conc_flags.patch

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -f bwa
	rm -rf bwa-*.tar.bz2
	rm -rf bwa-*
