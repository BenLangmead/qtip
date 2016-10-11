#
# bowtie2
#

USER=BenLangmead
VER=2.2.9
NM=bowtie2
NMVER=$(NM)-$(VER)
AR=$(NMVER)-source.zip
URL=https://github.com/$(USER)/$(NM)/releases/download/v$(VER)/$(AR)

.PHONY: all
all: $(NM) $(NM)-build

$(NM)-build: $(NM)-dir/bowtie2-build-s
	ln -s -f $(NM)-dir/$(NM)-build $@

$(NM): $(NM)-dir/bowtie2-align-s
	ln -s -f $(NM)-dir/$(NM) $@

$(NM)-dir/bowtie2-build-s: $(NM)-dir
	make -C $(NM)-dir bowtie2-build-s

$(NM)-dir/bowtie2-align-s: $(NM)-dir
	make -C $(NM)-dir bowtie2-align-s

$(NM)-dir: $(AR)
	unzip $<
	mv $(NMVER) $@
	cd $@ && patch -p1 < ../bowtie2_$(VER)_mapq.patch

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -f $(NM) $(NM)-build
	rm -rf $(AR) $(NM)-dir
