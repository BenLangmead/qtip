#!gmake

VER=1.0beta.18
AR=v$(VER).tar.gz
URL=https://github.com/amplab/snap/archive/$(AR)
NM=snap
NM_VER=snap-$(VER)

$(NM)/snap-aligner: $(NM)
	make -C $(NM) # CXX=g++-5

$(NM): $(NM_VER)
	mv $< $@
#	cd $@ && patch -p1 < ../snap_patch.diff
	cd $@ && patch -p1 < ../snap_features_patch.diff

$(NM_VER): $(AR)
	tar xvfz $<

$(AR):
	wget -O $(AR) $(URL)

.PHONY: clean
clean:
	rm -f bwa
	rm -f *.tar.gz

