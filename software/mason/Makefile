MASON_SRC=seqan-trunk/build/Release/bin/mason
RELEASE=seqan-v1.4.2

mason: $(MASON_SRC)
	ln -s -f seqan-trunk/build/Release/bin/mason

$(MASON_SRC): seqan-trunk/build
	cd seqan-trunk/build && mkdir -p Release && cd Release && cmake -DBoost_NO_BOOST_CMAKE=BOOL:ON -DCMAKE_BUILD_TYPE=Release ../.. && $(MAKE) mason

seqan-trunk/build: $(SVN)
	git clone --branch $(RELEASE) https://github.com/seqan/seqan.git seqan-trunk
	mkdir -p $@
