#
# bowtie2
#

VER=2.3.3.1
NM=bowtie2
NMVER=$(NM)-$(VER)
AR=$(NMVER)-source.zip
SF_PROJ=bowtie-bio
URL=http://downloads.sourceforge.net/project/$(SF_PROJ)/$(NM)/$(VER)/$(AR)

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

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -f $(NM) $(NM)-build
	rm -rf $(AR) $(NM)-dir
