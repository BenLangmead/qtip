#!gmake

VER=1.0.3
NM=pbsim-$(VER)
AR=$(NM).tar.gz
URL=http://pbsim.googlecode.com/files/$(AR)

pbsim: $(NM)/src/pbsim
	ln -s -f $< $@

$(NM)/src/pbsim: $(NM)
	cd $< && ./configure --prefix=`pwd`/install-root
	cd $< && make

$(NM): $(AR)
	tar zxvf $<

$(AR):
	wget $(URL)
