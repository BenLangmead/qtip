#!gmake

# BTL 11/22/2013: Can't get grinder to install.  Perhaps it just takes
# a very long time to download the BioPerl packages required?  Stalls
# at 'Installing Bio::Root::Version...'

VER=0.5.3
NM=Grinder-$(VER)
AR=$(NM).tar.gz
URL=http://downloads.sourceforge.net/project/biogrinder/biogrinder/Grinder-$(VER)/$(AR)

$(NM)/???: $(NM)
	cd $(NM) && perl Makefile.PL --defaultdeps
	cd $(NM) && make

$(NM): $(AR)
	tar zxvf $<

$(AR):
	wget $(URL)
