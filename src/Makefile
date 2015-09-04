#!gmake

.PHONY: all
all: qsim qsim-debug

qsim: qsim.cpp
	g++ -o $@ $^

qsim-debug: qsim.cpp
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -f qsim qsim-debug *.csv
