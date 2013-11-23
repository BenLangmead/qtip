#!gmake

# BTL 11/22/2013: Could not download files from BGI servers

VER=111
NM=pIRS
NM_VER=$(NM)_$(VER)
AR=$(NM_VER).tgz
URL=ftp://ftp.genomics.org.cn/pub/$(NM)/$(AR)

#$(NM)/???: $(NM)
#	cd $(NM) && perl Makefile.PL --defaultdeps
#	cd $(NM) && make

$(NM): $(AR)
	tar zxvf $<

$(AR):
	wget $(URL)
