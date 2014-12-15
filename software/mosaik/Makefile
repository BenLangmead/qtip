#!gmake

# might have to export BLD_PLATFORM=macosx64 first on Mac

VER=2.2.3
NM=MOSAIK
NM_VER=$(NM)-$(VER)
AR=$(NM)-source.tar
URL=https://mosaik-aligner.googlecode.com/files/$(AR)

.PHONY: all
all: MosaikBuild MosaikAligner

MosaikBuild: $(NM)/bin/MosaikBuild
	ln -s -f $< $@

MosaikAligner: $(NM)/bin/MosaikAligner
	ln -s -f $< $@

$(NM)/bin/MosaikBuild $(NM)/bin/MosaikAligner: $(NM)
	$(MAKE) -C $(NM)
	mv bin $(NM)
	mv obj $(NM)

$(NM): $(AR)
	tar xvf $<
	mv $(NM_VER)-source $@

$(AR):
	wget --no-check $(URL)

.PHONY: clean
clean:
	rm -f mosaik
	rm -rf MOSAIK-*
	rm -rf bin obj
	rm -f MosaikAligner MosaikBuild
