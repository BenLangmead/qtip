NM=htslib
VER=1.2.1
AR=$(NM)-$(VER).tar.bz2
URL=https://github.com/samtools/htslib/releases/download/$(VER)/$(AR)

$(NM)/install-root: $(NM)/Makefile
	cd $(NM) && ./configure --prefix=`pwd`/install-root
	make -C $(NM)
	make -C $(NM) install

$(NM)/Makefile: $(AR)
	tar xvfj $(AR)
	mv $(NM)-$(VER) $(NM)
	touch $@

$(AR):
	wget -O $(AR) $(URL)

.PHONY: clean
clean:
	rm -rf htslib*
