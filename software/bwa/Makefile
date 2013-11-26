#!gmake

# BTL 11/22/2013: Can't get grinder to install.  Perhaps it just takes
# a very long time to download the BioPerl packages required?  Stalls
# at 'Installing Bio::Root::Version...'

VER=0.7.5a
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
