#!gmake

VER=2014-12-06
NM=gsnap
NM_VER=$(NM)-$(VER)
AR=gmap-gsnap-$(VER).tar.gz
URL=http://research-pub.gene.com/gmap/src/$(AR)

$(NM): $(NM)-src/install-root/bin/gsnap
	ln -s -f $< $@

$(NM)-src/install-root/bin/gsnap: $(NM)-src/Makefile
	$(MAKE) -C $(NM)-src install

$(NM)-src/Makefile: $(NM)-src
	mkdir -p $(NM)-src/install-root
	cd $(NM)-src && ./configure --prefix=`pwd`/install-root

$(NM)-src: $(AR)
	tar xvf $<
	mv gmap-$(VER) $@

$(AR):
	wget --no-check $(URL)

.PHONY: clean
clean:
	rm -f $(AR) $(NM)
	rm -rf $(NM)-src
