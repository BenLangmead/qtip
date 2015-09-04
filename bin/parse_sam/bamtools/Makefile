#!gmake

NCPUS=8

test: test.cpp bamtools/include/api/BamReader.h
	g++ -I bamtools/include -L bamtools/lib -lbamtools -o test test.cpp

bamtools/include/api/BamReader.h: bamtools/build/Makefile
	cd bamtools/build && make -j8

bamtools/build/Makefile: bamtools
	cd bamtools && mkdir build && cd build && cmake ..

bamtools:
	git clone git@github.com:pezmaster31/bamtools.git

clean:
	rm -f test
	rm -rf bamtools
