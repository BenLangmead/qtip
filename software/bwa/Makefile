#!gmake

VER=0.7.15
NM=bwa-$(VER)
AR=$(NM).tar.bz2
URL=https://github.com/lh3/bwa/releases/download/v$(VER)/$(AR)

bwa: $(NM)/bwa
	ln -s -f $< $@

$(NM)/bwa: $(NM)
	make -C $(NM)

$(NM): $(AR) bwa_conc_flags_$(VER).patch
	tar jxvf $<
	cd $@ && patch -p1 < ../bwa_conc_flags_$(VER).patch

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -f bwa
	rm -rf bwa-*.tar.bz2
	rm -rf bwa-*
